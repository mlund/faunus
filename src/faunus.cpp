#define DOCTEST_CONFIG_IMPLEMENT
#include <doctest/doctest.h>
#include <doctest/trompeloeil.hpp>
#include "core.h"
#include "mpicontroller.h"
#include "move.h"
#include "montecarlo.h"
#include "analysis.h"
#include "multipole.h"
#include "docopt.h"
#include "progress_tracker.h"
#include <cstdlib>
#include "spdlog/spdlog.h"
#include <spdlog/sinks/null_sink.h>
#include <spdlog/sinks/stdout_color_sinks.h>
#include <spdlog/sinks/basic_file_sink.h>
#include <iomanip>
#include <unistd.h>
#include <chrono>

#ifdef ENABLE_SID
#include "cppsid.h"
#include <thread>
using namespace std::this_thread;     // sleep_for, sleep_until
using namespace std::chrono_literals; // ns, us, ms, s, h, etc.
using std::chrono::system_clock;
// some forward declarations
std::pair<std::string, int> findSIDsong();
std::shared_ptr<CPPSID::Player> createLoadedSIDplayer();
#endif

using namespace Faunus;
using namespace std;

#define Q(x) #x
#define QUOTE(x) Q(x)

#ifndef FAUNUS_TIPSFILE
#define FAUNUS_TIPSFILE ""
#endif

static const char USAGE[] =
    R"(Faunus - the Monte Carlo code you're looking for!

    https://faunus.readthedocs.io

    Usage:
      faunus [-q] [--verbosity <N>] [--nobar] [--nopfx] [--notips] [--nofun] [--state=<file>] [--input=<file>] [--output=<file>]
      faunus (-h | --help)
      faunus --version
      faunus test <doctest-options>...

    Options:
      -i <file> --input <file>   Input file [default: /dev/stdin].
      -o <file> --output <file>  Output file [default: out.json].
      -s <file> --state <file>   State file to start from (.json/.ubj).
      -v <N> --verbosity <N>     Log verbosity level (0 = off, 1 = critical, ..., 6 = trace) [default: 4]
      -q --quiet                 Less verbose output. It implicates -v0 --nobar --notips --nofun.
      -h --help                  Show this screen.
      --nobar                    No progress bar.
      --nopfx                    Do not prefix input file with MPI rank.
      --notips                   Do not give input assistance
      --nofun                    No fun
      --version                  Show version.

    Multiple processes using MPI:

    1. input and output files are prefixed with "mpi{rank}."
    2. standard output is redirected to "mpi{rank}.stdout"
    3. Input prefixing can be suppressed with --nopfx
)";

using ProgressIndicator::ProgressTracker;

// forward declarations
std::shared_ptr<ProgressTracker> createProgressTracker(bool, unsigned int);

int main(int argc, const char **argv) {
    if (argc > 1) { // run unittests if the first argument equals "test"
        if (std::string(argv[1]) == "test") {
#ifdef DOCTEST_CONFIG_DISABLE
            std::cout << "error: this version of faunus does not include unittests\n";
            return EXIT_FAILURE;
#else
            faunus_logger = spdlog::basic_logger_mt("faunus", "unittests.log", true);
            faunus_logger->set_pattern("%L: %v");
            faunus_logger->set_level(spdlog::level::debug);
            return doctest::Context(argc, argv).run();
#endif
        }
    }

    using namespace Faunus::MPI;
    bool quiet = false, nofun = true; // conservative defaults
    try {

        auto starting_time = std::chrono::steady_clock::now(); // used to time the simulation

        std::string version = "Faunus";
#ifdef GIT_LATEST_TAG
        version += " "s + QUOTE(GIT_LATEST_TAG);
#endif
#ifdef GIT_COMMIT_HASH
        version += " git " + std::string(GIT_COMMIT_HASH);
#endif
#ifdef ENABLE_SID
        version += " [sid]";
#endif
#ifdef ENABLE_MPI
        version += " [mpi]";
#endif
#ifdef _OPENMP
        version += " [openmp]";
#endif
        auto args = docopt::docopt(USAGE, {argv + 1, argv + argc}, true, version);

        // --quiet; set quiet first as it also affects exception handling
        quiet = args["--quiet"].asBool();
        if (quiet) {
            cout.setstate(std::ios_base::failbit); // hold kæft
        }

        mpi.init(); // initialize MPI, if available

        // prepare loggers
        // TODO refactor to a standalone function and use cmd line options for different sinks, etc.
        faunus_logger = spdlog::stderr_color_mt("faunus");
        faunus_logger->set_pattern("[%n %P] %^%L: %v%$");
        mcloop_logger = spdlog::stderr_color_mt("mcloop");
        mcloop_logger->set_pattern("[%n %P] [%E.%f] %L: %v");

        // --verbosity (log level)
        long log_level =
            spdlog::level::off - (quiet ? 0 : args["--verbosity"].asLong()); // reverse sequence 0 → 6 to 6 → 0
        spdlog::set_level(static_cast<spdlog::level::level_enum>(log_level));

        // --notips
        if (!quiet && !args["--notips"].asBool()) {
            usageTip.load({FAUNUS_TIPSFILE});
        }

        // --nobar
        bool show_progress = !quiet && !args["--nobar"].asBool();

        // --nofun
        nofun = args["--nofun"].asBool();
        usageTip.asciiart = !quiet && !nofun;
#ifdef ENABLE_SID
        usageTip.asciiart = false; // if SID is enabled, disable ascii
#endif

        // --nopfx
        bool prefix = !args["--nopfx"].asBool();

        // --input
        json json_in;
        try {
            if (auto input = args["--input"].asString(); input == "/dev/stdin") {
                    std::cin >> json_in;
            } else {
                if (prefix) {
                    input = Faunus::MPI::prefix + input;
                }
                json_in = openjson(input);
            }
        } catch(json::parse_error& e) {
            faunus_logger->debug(e.what());
            throw ConfigurationError("empty or invalid input JSON");
        }

        {
            pc::temperature = json_in.at("temperature").get<double>() * 1.0_K;
            MetropolisMonteCarlo sim(json_in, mpi);

            // --state
            if (args["--state"]) {
                std::ifstream f;
                std::string state = Faunus::MPI::prefix + args["--state"].asString();
                std::string suffix = state.substr(state.find_last_of(".") + 1);
                bool binary = (suffix == "ubj");
                auto mode = std::ios::in;
                if (binary) {
                    mode = std::ifstream::ate | std::ios::binary; // ate = open at end
                }
                f.open(state, mode);
                if (f) {
                    json j;
                    faunus_logger->info("loading state file {}", state);
                    if (binary) {
                        size_t size = f.tellg(); // get file size
                        std::vector<std::uint8_t> v(size / sizeof(std::uint8_t));
                        f.seekg(0, f.beg); // go back to start
                        f.read((char *)v.data(), size);
                        j = json::from_ubjson(v);
                    } else {
                        f >> j;
                    }
                    sim.restore(j);
                } else {
                    throw std::runtime_error("state file error: " + state);
                }
            }

            // warn if initial system has a net charge
            {
                auto p = sim.getSpace().activeParticles();
                if (double system_charge = Faunus::monopoleMoment(p.begin(), p.end()); std::fabs(system_charge) > 0) {
                    faunus_logger->warn("non-zero system charge of {}e", system_charge);
                }
            }

            Analysis::CombinedAnalysis analysis(json_in.at("analysis"), sim.getSpace(), sim.getHamiltonian());

            auto &loop = json_in.at("mcloop");
            int macro = loop.at("macro");
            int micro = loop.at("micro");

            auto progress_tracker = createProgressTracker(show_progress, macro * micro);
            for (int i = 0; i < macro; i++) {
                for (int j = 0; j < micro; j++) {
                    if (progress_tracker && mpi.isMaster()) {
                        if(++(*progress_tracker) % 10 == 0) {
                            progress_tracker->display();
                        }
                    }
                    sim.move();
                    analysis.sample();
                }                   // end of micro steps
                analysis.to_disk(); // save analysis to disk
            }                       // end of macro steps
            if (progress_tracker && mpi.isMaster()) {
                progress_tracker->done();
            }

            faunus_logger->log((sim.relativeEnergyDrift() < 1E-9) ? spdlog::level::info : spdlog::level::warn,
                               "relative energy drift = {}", sim.relativeEnergyDrift());

            // --output
            if (std::ofstream file(Faunus::MPI::prefix + args["--output"].asString()); file) {
                json j;
                Faunus::to_json(j, sim);
                j["relative drift"] = sim.relativeEnergyDrift();
                j["analysis"] = analysis;
                if (mpi.nproc() > 1) {
                    j["mpi"] = mpi;
                }
#ifdef GIT_COMMIT_HASH
                j["git revision"] = GIT_COMMIT_HASH;
#endif
#ifdef __VERSION__
                j["compiler"] = __VERSION__;
#endif

                { // report on total simulation time
                    using namespace std::chrono;
                    auto ending_time = steady_clock::now();
                    auto secs = duration_cast<seconds>(ending_time - starting_time).count();
                    j["simulation time"] = {{"in minutes", secs / 60.0}, {"in seconds", secs}};
                }

                file << std::setw(4) << j << std::endl;
            }
        }

        mpi.finalize();

    } catch (std::exception &e) {
        faunus_logger->error(e.what());
        std::cerr << e.what() << std::endl;

        // ConfigurationError can carry a JSON snippet which should be shown for debugging.
        if (auto config_error = dynamic_cast<ConfigurationError*>(&e);
            config_error != nullptr && !config_error->attachedJson().empty()) {
            faunus_logger->debug("JSON snippet:\n{}", config_error->attachedJson().dump(4));
        }

        if (!usageTip.buffer.empty()) {
            // Use the srderr stream directly for more elaborated output of usage tip, optionally containing an ASCII
            // art. Level info seems appropriate.
            if (faunus_logger->should_log(spdlog::level::info)) {
                std::cerr << usageTip.buffer << std::endl;
            }
        }
#ifdef ENABLE_SID
        // easter egg...
        if (!quiet && !nofun && mpi.isMaster()) { // -> fun

            auto player = createLoadedSIDplayer(); // create C64 SID emulation and load a random tune

            if (player) {
                faunus_logger->info("error message music '{}' by {}, {} (6502/SID emulation)", player->title(),
                                    player->author(), player->info());
                faunus_logger->info("\033[1mpress ctrl-c to quit\033[0m");
                player->start();                         // start music
                sleep_for(10ns);                         // short delay
                sleep_until(system_clock::now() + 240s); // play for 4 minutes, then exit
                player->stop();
                std::cout << std::endl;
            }
        } // end of easter egg
#endif
        return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
}

#ifdef ENABLE_SID
/*
 * finds a random SID music file and picks a random sub-song if available
 *
 * will look for `music.json` as well as sid files in the following
 * directories and order:
 *
 * - BINARY_DIR/sids/
 * - INSTALL_PREFIX/share/faunus/sids
 */
std::pair<std::string, int> findSIDsong() {
    std::string filename;
    int subsong = -1;
    try {
        // look for json file with hvsc sid tune names
        std::string pfx;
        json json_music;
        for (std::string dir :
             {FAUNUS_BINARY_DIR, FAUNUS_INSTALL_PREFIX "/share/faunus/"}) { // installed and uninstalled cmake builds
            json_music = Faunus::openjson(dir + "/sids/music.json", false);
            if (!json_music.empty()) {
                pfx = dir + "/";
                break;
            }
        }
        if (not json_music.empty()) {
            json_music = json_music.at("songs"); // load playlist

            std::vector<size_t> weight; // weight for each tune (= number of subsongs)
            for (auto &i : json_music)
                weight.push_back(i.at("subsongs").size());
            std::discrete_distribution<size_t> dist(weight.begin(), weight.end());

            Faunus::random.seed();                                        // give global random a hardware seed
            auto it = json_music.begin() + dist(Faunus::random.engine);   // pick random tune weighted by subsongs
            auto subsongs = (*it).at("subsongs").get<std::vector<int>>(); // all subsongs
            subsong = *(Faunus::random.sample(subsongs.begin(), subsongs.end())) - 1; // random subsong
            filename = pfx + it->at("file").get<std::string>();
        }
    } catch (const std::exception &) {
        // silently ignore if something fails; it's just for fun!
    }
    return {filename, subsong};
}

std::shared_ptr<CPPSID::Player> createLoadedSIDplayer() {
    std::shared_ptr<CPPSID::Player> player;
    if (isatty(fileno(stdout))) {                        // only play music if on console
        if (not std::getenv("SSH_CLIENT")) {             // and not through a ssh connection
            player = std::make_shared<CPPSID::Player>(); // let's emulate a Commodore 64...
            auto tune = findSIDsong();                   // pick a song from our pre-defined library
            player->load(tune.first, tune.second);
        }
    }
    return player;
}
#endif

std::shared_ptr<ProgressTracker> createProgressTracker(bool show_progress, unsigned int steps) {
    using namespace ProgressIndicator;
    using namespace std::chrono;
    std::shared_ptr<ProgressTracker> tracker = nullptr;
    if(show_progress) {
        if (isatty(fileno(stdout))) {
            // show a progress bar on the console
            tracker = std::make_shared<ProgressBar>(steps);
        } else {
            // not in a console
            tracker = std::make_shared<TaciturnDecorator>(
                // hence print a new line
                std::make_shared<ProgressLog>(steps),
                // at most every 10 minutes or after 0.5% of progress, whatever comes first
                duration_cast<milliseconds>(minutes(10)), 0.005);
        }
    }
    return tracker;
}
