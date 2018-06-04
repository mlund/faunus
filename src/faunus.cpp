#include "core.h"
#include "move.h"
#include "mpi.h"
#include "multipole.h"
#include "docopt.h"
#include <cstdlib>
#include "ProgressBar.hpp"

using namespace Faunus;
using namespace std;

#ifdef CUBOID_NOPBC
typedef Geometry::CuboidNoPBC Tgeometry;
#elif SPHERE
typedef Geometry::Sphere Tgeometry;
#else
typedef Geometry::Cuboid Tgeometry;
#endif
typedef Particle<Charge> Tparticle;

static const char USAGE[] =
R"(Faunus - the Monte Carlo code you're looking for!

    http://github.com/mlund/neofaunus

    Usage:
      faunus [-q] [--nobar] [--nopfx] [--state=<file>] [--input=<file>] [--output=<file>] [--rerun=<traj>]
      faunus (-h | --help)
      faunus --version

    Options:
      -i <file> --input <file>   Input file [default: /dev/stdin].
      -o <file> --output <file>  Output file [default: out.json].
      -s <file> --state <file>   State file to start from.
      --rerun <traj>             Rerun from trajectory (.xtc)
      -q --quiet                 Less verbose output.
      -h --help                  Show this screen.
      --nobar                    No progress bar.
      --nopfx                    Do not prefix input file with MPI rank.
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
        std::string version="Faunus Mk2";
#ifdef GIT_COMMIT_HASH
        version += " git " + std::string(GIT_COMMIT_HASH);
#endif
        auto args = docopt::docopt( USAGE,
                { argv + 1, argv + argc }, true, version);

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

        pc::temperature = j.at("temperature").get<double>() * 1.0_K;
        MCSimulation<Tgeometry,Tparticle> sim(j, mpi);

        // --state
        if (args["--state"]) {
            std::string state = Faunus::MPI::prefix + args["--state"].asString();
            std::ifstream f(state);
            if (f) {
                if (!quiet)
                    mpi.cout() << "Loading state file '" << state << "'" << endl;
                json j;
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

                if (showProgress && mpi.isMaster()) {
                    ++progressBar;
                    if (j % 10 == 0)
                        progressBar.display();
                }

                sim.move();
                analysis.sample();
            }
        }
        if (showProgress && mpi.isMaster())
            progressBar.done();

        if (!quiet)
            mpi.cout() << "relative drift = " << sim.drift() << endl;

        // --output
        std::ofstream f(Faunus::MPI::prefix + args["--output"].asString());
        if (f) {
            json j = sim;
            j["analysis"] = analysis;
            if (mpi.nproc()>1)
                j["mpi"] = mpi;
#ifdef GIT_COMMIT_HASH
            j["git revision"] = GIT_COMMIT_HASH;
#endif
            f << std::setw(4) << j << endl;
        }

    } catch (std::exception &e) {
        std::cerr << e.what() << endl;
        return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
}
