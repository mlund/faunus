#include <iomanip>

#include "core.h"
#include "mpi.h"
#include "move.h"
#include "analysis.h"
#include "multipole.h"
#include "docopt.h"
#include <cstdlib>
#include "ProgressBar.hpp"

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
      faunus [-q] [--nobar] [--nopfx] [--notips] [--state=<file>] [--input=<file>] [--output=<file>]
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
      --version                  Show version.

    Multiple processes using MPI:

    1. input and output files are prefixed with "mpi{rank}."
    2. standard output is redirected to "mpi{rank}.stdout"
    3. Input prefixing can be suppressed with --nopfx
)";

int main( int argc, char **argv )
{
    using namespace Faunus::MPI;
    try {
        std::string version="Faunus";
#ifdef GIT_LATEST_TAG
        version += " "s + QUOTE(GIT_LATEST_TAG);
#endif
#ifdef GIT_COMMIT_HASH
        version += " git " + std::string(GIT_COMMIT_HASH);
#endif
        auto args = docopt::docopt( USAGE,
                { argv + 1, argv + argc }, true, version);

        mpi.init(); // initialize MPI, if available

        // --notips
        if (not args["--notips"].asBool())
            usageTip.load( {FAUNUS_TIPSFILE} ); 

        // --nobar
        bool showProgress = !args["--nobar"].asBool();

        // --quiet
        bool quiet = args["--quiet"].asBool();
        if (quiet) {
            cout.setstate( std::ios_base::failbit ); // hold kæft
            showProgress = false;
        }

        // --nopfx
        bool prefix = !args["--nopfx"].asBool();

        // --input
        json j;
        auto input = args["--input"].asString();
        if (input=="/dev/stdin")
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
                bool binary = (suffix=="ubj");
                auto mode = std::ios::in;
                if (binary)
                    mode = std::ifstream::ate | std::ios::binary; // ate = open at end
                f.open(state, mode);
                if (f) {
                    json j;
                    if (not quiet)
                        mpi.cout() << "Loading state file '" << state << "'" << endl;
                    if (binary) {
                        size_t size = f.tellg(); // get file size
                        std::vector<std::uint8_t> v(size/sizeof(std::uint8_t));
                        f.seekg(0, f.beg);       // go back to start
                        f.read((char*)v.data(), size);
                        j = json::from_ubjson(v);
                    } else
                        f >> j;
                    sim.restore(j);
                } else
                    throw std::runtime_error("Error loading state file '" + state + "'");
            }

            Analysis::CombinedAnalysis analysis(j.at("analysis"), sim.space(), sim.pot());

            auto& loop = j.at("mcloop");
            int macro = loop.at("macro");
            int micro = loop.at("micro");

            ProgressBar progressBar(macro*micro, 70);
            for (int i=0; i<macro; i++) {
                for (int j=0; j<micro; j++) {

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

            if (not quiet)
                mpi.cout() << "relative drift = " << sim.drift() << endl;

            // --output
            std::ofstream f(Faunus::MPI::prefix + args["--output"].asString());
            if (f) {
                json j;
                Faunus::to_json(j, sim);
                j["relative drift"] = sim.drift();
                j["analysis"] = analysis;
                if (mpi.nproc()>1)
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
        std::cerr << e.what() << endl;
#ifdef ENABLE_SID
        // easter egg...
        CPPSID::Player player;
        player.load("Beyond_the_Zero.sid", 0);
        player.start();
        sleep_for(10ns);
        sleep_until(system_clock::now() + 120s);
        player.stop();
#endif
        return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
}
