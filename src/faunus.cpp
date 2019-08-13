#include <iomanip>

#include "core.h"
#include "mpi.h"
#include "move.h"
#include "montecarlo.h"
#include "analysis.h"
//#include "multipole.h"
#include "docopt.h"
#include <cstdlib>
#include "ProgressBar.hpp"
#include "spdlog/spdlog.h"
#include <spdlog/sinks/null_sink.h>
#include <spdlog/sinks/stdout_color_sinks.h>

#ifdef ENABLE_SID
#include "cppsid.h"
#include <thread>
using namespace std::this_thread;     // sleep_for, sleep_until
using namespace std::chrono_literals; // ns, us, ms, s, h, etc.
using std::chrono::system_clock;
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

    http://github.com/mlund/faunus

    Usage:
      faunus [-q] [--nobar] [--nopfx] [--notips] [--nofun] [--state=<file>] [--input=<file>] [--output=<file>]
      faunus (-h | --help)
      faunus --version

    Options:
      -i <file> --input <file>   Input file [default: /dev/stdin].
      -o <file> --output <file>  Output file [default: out.json].
      -s <file> --state <file>   State file to start from (.json/.ubj).
      -q --quiet                 Less verbose output.
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

int main(int argc, char **argv) {
    using namespace Faunus::MPI;
    bool nofun;
    try {
        std::string version = "Faunus";
#ifdef GIT_LATEST_TAG
        version += " "s + QUOTE(GIT_LATEST_TAG);
#endif
#ifdef GIT_COMMIT_HASH
        version += " git " + std::string(GIT_COMMIT_HASH);
#endif
        auto args = docopt::docopt(USAGE, {argv + 1, argv + argc}, true, version);

        mpi.init(); // initialize MPI, if available

        // prepare loggers
        // TODO refactor to a standalone function and use cmd line options for different sinks, etc.
        faunus_logger = spdlog::stderr_color_mt("faunus");
        faunus_logger->set_pattern("[%n %P] %^%L: %v%$");
        faunus_logger->set_level(spdlog::level::info);
        mcloop_logger = spdlog::stderr_color_mt("mcloop");
        mcloop_logger->set_pattern("[%n %P] [%E.%f] %L: %v");
        mcloop_logger->set_level(spdlog::level::warn);

        // --notips
        if (not args["--notips"].asBool())
            usageTip.load({FAUNUS_TIPSFILE});

        // --nobar
        bool showProgress = !args["--nobar"].asBool();

        // --quiet
        bool quiet = args["--quiet"].asBool();
        if (quiet) {
            cout.setstate(std::ios_base::failbit); // hold kæft
            showProgress = false;
        }

        // --nofun
        nofun = args["--nofun"].asBool();
        if (nofun or quiet)
            usageTip.asciiart = false;
#ifdef ENABLE_SID
        usageTip.asciiart = false; // if SID is enabled, disable ascii
#endif

        // --nopfx
        bool prefix = !args["--nopfx"].asBool();

        // --input
        json j;
        auto input = args["--input"].asString();
        if (input == "/dev/stdin")
            std::cin >> j;
        else {
            if (prefix)
                input = Faunus::MPI::prefix + input;
            j = openjson(input);
        }

        {
            pc::temperature = j.at("temperature").get<double>() * 1.0_K;
            MCSimulation sim(j, mpi);

            // --state
            if (args["--state"]) {
                std::ifstream f;
                std::string state = Faunus::MPI::prefix + args["--state"].asString();
                std::string suffix = state.substr(state.find_last_of(".") + 1);
                bool binary = (suffix == "ubj");
                auto mode = std::ios::in;
                if (binary)
                    mode = std::ifstream::ate | std::ios::binary; // ate = open at end
                f.open(state, mode);
                if (f) {
                    json j;
                    if (not quiet)
                        faunus_logger->info("loading state file {}", state);
                    if (binary) {
                        size_t size = f.tellg(); // get file size
                        std::vector<std::uint8_t> v(size / sizeof(std::uint8_t));
                        f.seekg(0, f.beg); // go back to start
                        f.read((char *)v.data(), size);
                        j = json::from_ubjson(v);
                    } else
                        f >> j;
                    sim.restore(j);
                } else
                    throw std::runtime_error("state file error: " + state);
            }

            Analysis::CombinedAnalysis analysis(j.at("analysis"), sim.space(), sim.pot());

            auto &loop = j.at("mcloop");
            int macro = loop.at("macro");
            int micro = loop.at("micro");

            ProgressBar progressBar(macro * micro, 70);
            for (int i = 0; i < macro; i++) {
                for (int j = 0; j < micro; j++) {

                    if (showProgress and mpi.isMaster()) {
                        ++progressBar;
                        if (j % 10 == 0)
                            progressBar.display();
                    }

                    sim.move();
                    analysis.sample();
                }
            }
            if (showProgress and mpi.isMaster())
                progressBar.done();

            if (not quiet) {
                faunus_logger->log((sim.drift() < 1E-9) ? spdlog::level::info : spdlog::level::warn,
                                   "relative drift = {}", sim.drift());
            }

            // --output
            std::ofstream f(Faunus::MPI::prefix + args["--output"].asString());
            if (f) {
                json j;
                Faunus::to_json(j, sim);
                j["relative drift"] = sim.drift();
                j["analysis"] = analysis;
                if (mpi.nproc() > 1)
                    j["mpi"] = mpi;
#ifdef GIT_COMMIT_HASH
                j["git revision"] = GIT_COMMIT_HASH;
#endif
#ifdef __VERSION__
                j["compiler"] = __VERSION__;
#endif
                f << std::setw(4) << j << endl;
            }
        }

        mpi.finalize();

    } catch (std::exception &e) {
        faunus_logger->error(e.what());

        if (not usageTip.buffer.empty())
            faunus_logger->info(usageTip.buffer);
#ifdef ENABLE_SID
        // easter egg...
        if (not nofun) { // -> fun
            try {
                // look for json file with hvsc sid tune names
                std::string pfx;
                json j;
                for (std::string dir : {FAUNUS_INSTALL_PREFIX "/share/faunus/",
                                        FAUNUS_BINARY_DIR}) { // installed and uninstalled cmake builds
                    j = Faunus::openjson(dir + "/sids/music.json", false);
                    if (not j.empty()) {
                        pfx = dir + "/";
                        break;
                    }
                }
                if (not j.empty()) {
                    j = j.at("songs");                                            // load playlist
                    Faunus::random.seed();                                        // give global random a hardware seed
                    auto it = Faunus::random.sample(j.begin(), j.end());          // pick a random song
                    auto subsongs = (*it).at("subsongs").get<std::vector<int>>(); // all subsongs
                    int subsong = *(Faunus::random.sample(subsongs.begin(), subsongs.end())) - 1; // random subsong

                    CPPSID::Player player; // let's emulate a Commodore 64...

                    if (player.load(pfx + it->at("file").get<std::string>(), subsong)) {
                        std::cout << "You're comforted by C64 music '" << player.title() << "' by " << player.author()
                                  << ", " << player.info() << "\n\n"
                                  << "Press ctrl-c to quit." << std::flush;
                        player.start();                          // start music
                        sleep_for(10ns);                         // short delay
                        sleep_until(system_clock::now() + 240s); // play for 4 minutes, then exit
                        player.stop();
                        std::cout << std::endl;
                    }
                }
            } catch (const std::exception &) {
                // silently ignore if something fails; it's just for fun
            }
        } // end of fun
#endif
        return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
}
