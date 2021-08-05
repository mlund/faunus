#define DOCTEST_CONFIG_IMPLEMENT
#include <doctest/doctest.h>
#include "mpicontroller.h"
#include "montecarlo.h"
#include "analysis.h"
#include "multipole.h"
#include "docopt.h"
#include "progress_tracker.h"
#include "spdlog/spdlog.h"
#include <spdlog/sinks/null_sink.h>
#include <spdlog/sinks/stdout_color_sinks.h>
#include <spdlog/sinks/basic_file_sink.h>

#ifdef ENABLE_SID
#include "cppsid.h"
#include <thread>
#include <unistd.h>
#include <chrono>
std::pair<std::string, int> findSIDsong();
std::shared_ptr<CPPSID::Player> createLoadedSIDplayer();
#endif

using namespace Faunus;

#define Q(x) #x
#define QUOTE(x) Q(x)

#ifndef FAUNUS_TIPSFILE
#define FAUNUS_TIPSFILE ""
#endif

std::shared_ptr<ProgressIndicator::ProgressTracker> createProgressTracker(bool, unsigned int);
int runUnittests(int argc, const char* const* argv);
std::string versionString();
void setInformationLevelAndLoggers(bool quiet, docopt::Options& args);
json getUserInput(docopt::Options& args);
void loadState(docopt::Options& args, MetropolisMonteCarlo& simulation);
void checkElectroNeutrality(MetropolisMonteCarlo& simulation);
void showErrorMessage(std::exception& exception);
void playRetroMusic();
template <typename TimePoint>
void saveOutput(TimePoint& starting_time, docopt::Options& args, MetropolisMonteCarlo& simulation,
                const Analysis::CombinedAnalysis& analysis);

void mainLoop(bool show_progress, const json& json_in, MetropolisMonteCarlo& simulation,
              Analysis::CombinedAnalysis& analysis);

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

int main(int argc, const char** argv) {
    if (argc > 1) { // run unittests if the first argument equals "test"
        if (std::string(argv[1]) == "test") {
            return runUnittests(argc, argv);
        }
    }
    bool fun = false;   //!< enable utterly unnecessarily stuff?
    bool quiet = false; //!< hold kaje?
    try {
        Faunus::MPI::mpi.init(); // initialize MPI, if available
        const auto starting_time = std::chrono::steady_clock::now();
        const auto version = versionString();
        auto args = docopt::docopt(USAGE, {argv + 1, argv + argc}, true, version);
        fun = !args["--nofun"].asBool();
        quiet = args["--quiet"].asBool();
        setInformationLevelAndLoggers(quiet, args);

        const auto input = getUserInput(args);
        MetropolisMonteCarlo simulation(input, Faunus::MPI::mpi);
        loadState(args, simulation);
        checkElectroNeutrality(simulation);
        Analysis::CombinedAnalysis analysis(input.at("analysis"), simulation.getSpace(), simulation.getHamiltonian());

        const bool show_progress = !quiet && !args["--nobar"].asBool();
        mainLoop(show_progress, input, simulation, analysis); // run simulation!

        saveOutput(starting_time, args, simulation, analysis);
        Faunus::MPI::mpi.finalize();
        return EXIT_SUCCESS;

    } catch (std::exception& e) {
        showErrorMessage(e);
        if (!quiet && fun) {
            playRetroMusic();
        }
        return EXIT_FAILURE;
    }
}

void showErrorMessage(std::exception& exception) {
    faunus_logger->error(exception.what());
    std::cerr << exception.what() << std::endl;

    // ConfigurationError can carry a JSON snippet which should be shown for debugging.
    if (auto* config_error = dynamic_cast<ConfigurationError*>(&exception);
        config_error != nullptr && !config_error->attachedJson().empty()) {
        faunus_logger->debug("json snippet:\n{}", config_error->attachedJson().dump(4));
    }
    if (!usageTip.buffer.empty()) {
        // Use the srderr stream directly for more elaborated output of usage tip, optionally containing an ASCII
        // art. Level info seems appropriate.
        if (faunus_logger->should_log(spdlog::level::info)) {
            std::cerr << usageTip.buffer << std::endl;
        }
    }
}

void playRetroMusic() {
#ifdef ENABLE_SID
    using std::chrono_literals::operator""s;
    using std::chrono_literals::operator""ns;
    if (MPI::mpi.isMaster()) {
        if (auto player = createLoadedSIDplayer()) { // create C64 SID emulation and load a random tune
            faunus_logger->info("error message music '{}' by {}, {} (6502/SID emulation)", player->title(),
                                player->author(), player->info());
            faunus_logger->info("\033[1mpress ctrl-c to quit\033[0m");
            player->start();                                                        // start music
            std::this_thread::sleep_for(10ns);                                      // short delay
            std::this_thread::sleep_until(std::chrono::system_clock::now() + 240s); // play for 4 minutes, then exit
            player->stop();
            std::cout << std::endl;
        }
    }
#endif
}

std::string versionString() {
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
    return version;
}

int runUnittests(int argc, const char* const* argv) {
#ifdef DOCTEST_CONFIG_DISABLE
    std::cerr << "this version of faunus does not include unittests" << std::endl;
    return EXIT_FAILURE;
#else
    faunus_logger = spdlog::basic_logger_mt("faunus", "unittests.log", true);
    faunus_logger->set_pattern("%L: %v");
    faunus_logger->set_level(spdlog::level::debug);
    return doctest::Context(argc, argv).run();
#endif
}

void mainLoop(bool show_progress, const json& json_in, MetropolisMonteCarlo& simulation,
              Analysis::CombinedAnalysis& analysis) {
    const auto& loop = json_in.at("mcloop");
    const auto macro = loop.at("macro").get<int>();
    const auto micro = loop.at("micro").get<int>();
    auto progress_tracker = createProgressTracker(show_progress, macro * micro);
    for (int i = 0; i < macro; i++) {
        for (int j = 0; j < micro; j++) {
            if (progress_tracker && MPI::mpi.isMaster()) {
                if (++(*progress_tracker) % 10 == 0) {
                    progress_tracker->display();
                }
            }
            simulation.move();
            analysis.sample();
        }                   // end of micro steps
        analysis.to_disk(); // save analysis to disk
    }                       // end of macro steps
    if (progress_tracker && MPI::mpi.isMaster()) {
        progress_tracker->done();
    }
    const double drift_tolerance = 1e-9;
    const auto level = simulation.relativeEnergyDrift() < drift_tolerance ? spdlog::level::info : spdlog::level::warn;
    faunus_logger->log(level, "relative energy drift = {:.3E}", simulation.relativeEnergyDrift());
}

void checkElectroNeutrality(MetropolisMonteCarlo& simulation) {
    auto particles = simulation.getSpace().activeParticles();
    const auto system_charge = Faunus::monopoleMoment(particles.begin(), particles.end());
    const double max_allowed_charge = 1e-6;
    if (std::fabs(system_charge) > max_allowed_charge) {
        faunus_logger->warn("non-zero system charge of {:.3E}e", system_charge);
    }
}

void setInformationLevelAndLoggers(bool quiet, docopt::Options& args) {
    // @todo refactor to use cmd line options for different sinks, etc.
    faunus_logger = spdlog::stderr_color_mt("faunus");
    faunus_logger->set_pattern("[%n %P] %^%L: %v%$");
    mcloop_logger = spdlog::stderr_color_mt("mcloop");
    mcloop_logger->set_pattern("[%n %P] [%E.%f] %L: %v");

    const long log_level = spdlog::level::off - (quiet ? 0 : args["--verbosity"].asLong()); // reverse: 0 → 6 to 6 → 0
    spdlog::set_level(static_cast<spdlog::level::level_enum>(log_level));
    if (quiet) {
        std::cout.setstate(std::ios_base::failbit); // muffle stdout
    }
    if (!quiet && !args["--notips"].asBool()) {
        usageTip.load({FAUNUS_TIPSFILE});
    }
    const bool fun = !args["--nofun"].asBool();
    usageTip.asciiart = !quiet && fun;
#ifdef ENABLE_SID
    usageTip.asciiart = false; // if SID is enabled, disable ascii
#endif
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
        for (std::string dir : {FAUNUS_BINARY_DIR, FAUNUS_INSTALL_PREFIX "/share/faunus/"}) {
            try {
                // look at installed and uninstalled cmake builds
                json_music = Faunus::openjson(dir + "/sids/music.json");
                if (!json_music.empty()) {
                    pfx = dir + "/";
                    break;
                }
            } catch (...) {
                // ignore any error
            }
        }
        if (not json_music.empty()) {
            json_music = json_music.at("songs"); // load playlist
            std::vector<size_t> weight;          // weight for each tune (= number of subsongs)
            for (const auto& song : json_music) {
                weight.push_back(song.at("subsongs").size());
            }
            std::discrete_distribution<size_t> dist(weight.begin(), weight.end());
            Faunus::random.seed(); // give global random a hardware seed
            const auto random_song =
                json_music.begin() + dist(Faunus::random.engine); // pick random tune weighted by subsongs
            const auto subsongs = random_song->at("subsongs").get<std::vector<int>>(); // all subsongs
            subsong = *(Faunus::random.sample(subsongs.begin(), subsongs.end())) - 1;  // random subsong
            filename = pfx + random_song->at("file").get<std::string>();
        }
    } catch (const std::exception&) {
        // silently ignore if something fails; it's just for fun!
    }
    return {filename, subsong};
}

std::shared_ptr<CPPSID::Player> createLoadedSIDplayer() {
    std::shared_ptr<CPPSID::Player> player;
    if (static_cast<bool>(isatty(fileno(stdout)))) {         // only play music if on console
        if (!static_cast<bool>(std::getenv("SSH_CLIENT"))) { // and not through a ssh connection
            player = std::make_shared<CPPSID::Player>();     // let's emulate a Commodore 64...
            const auto [filename, subsong] = findSIDsong();  // pick a song from our pre-defined library
            player->load(filename, subsong);
        }
    }
    return player;
}
#endif

std::shared_ptr<ProgressIndicator::ProgressTracker> createProgressTracker(bool show_progress, unsigned int steps) {
    using std::chrono::milliseconds;
    using std::chrono::minutes;
    std::shared_ptr<ProgressIndicator::ProgressTracker> tracker = nullptr;
    if (show_progress) {
        if (static_cast<bool>(isatty(fileno(stdout)))) { // show a progress bar on the console
            tracker = std::make_shared<ProgressIndicator::ProgressBar>(steps);
        } else { // not in a console
            tracker = std::make_shared<ProgressIndicator::TaciturnDecorator>(
                // hence print a new line
                std::make_shared<ProgressIndicator::ProgressLog>(steps),
                // at most every 10 minutes or after 0.5% of progress, whatever comes first
                std::chrono::duration_cast<milliseconds>(minutes(10)), 0.005);
        }
    }
    return tracker;
}

json getUserInput(docopt::Options& args) {
    try {
        json j;
        if (auto filename = args["--input"].asString(); filename == "/dev/stdin") {
            std::cin >> j;
        } else {
            if (!args["--nopfx"].asBool()) {
                filename = MPI::prefix + filename;
            }
            j = openjson(filename);
        }
        pc::temperature = j.at("temperature").get<double>() * 1.0_K;
        return j;
    } catch (json::parse_error& e) {
        faunus_logger->error(e.what());
        throw ConfigurationError("invalid input -> {}", e.what());
    }
}

void loadState(docopt::Options& args, MetropolisMonteCarlo& simulation) {
    if (args["--state"]) {
        const auto statefile = MPI::prefix + args["--state"].asString();
        const auto suffix = statefile.substr(statefile.find_last_of('.') + 1);
        const bool binary = (suffix == "ubj");
        auto mode = std::ios_base::in;
        if (binary) {
            mode = std::ios_base::ate | std::ios_base::binary; // ate = open at end
        }
        if (auto stream = std::ifstream(statefile, mode)) {
            json j;
            faunus_logger->info("loading state file {}", statefile);
            if (binary) {
                const auto size = stream.tellg(); // get file size
                std::vector<std::uint8_t> buffer(size / sizeof(std::uint8_t));
                stream.seekg(0, stream.beg);             // go back to start...
                stream.read((char*)buffer.data(), size); // ...and read into buffer
                j = json::from_ubjson(buffer);
            } else {
                stream >> j;
            }
            simulation.restore(j);
        } else {
            throw std::runtime_error("state file error -> "s + statefile);
        }
    }
}

template <typename TimePoint>
void saveOutput(TimePoint& starting_time, docopt::Options& args, MetropolisMonteCarlo& simulation,
                const Analysis::CombinedAnalysis& analysis) {

    if (std::ofstream stream(MPI::prefix + args["--output"].asString()); stream) {
        json j;
        to_json(j, simulation);
        j["relative drift"] = simulation.relativeEnergyDrift();
        j["analysis"] = analysis;
        if (MPI::mpi.nproc() > 1) {
            j["mpi"] = MPI::mpi;
        }
#ifdef GIT_COMMIT_HASH
        j["git revision"] = GIT_COMMIT_HASH;
#endif
#ifdef __VERSION__
        j["compiler"] = __VERSION__;
#endif

        { // report on total simulation time
            const auto ending_time = std::chrono::steady_clock::now();
            const auto elapsed_seconds =
                std::chrono::duration_cast<std::chrono::seconds>(ending_time - starting_time).count();
            j["simulation time"] = {{"in minutes", elapsed_seconds / 60.0}, {"in seconds", elapsed_seconds}};
        }

        stream << std::setw(2) << j << std::endl;
    }
}
