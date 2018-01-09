#include "core.h"
#include "move.h"
#include "multipole.h"
#include "mpi.h"
#include "docopt.h"
#include <cstdlib>
#include "ProgressBar.hpp"

using namespace Faunus;
using namespace std;

typedef Geometry::Cuboid Tgeometry;
typedef Particle<Charge> Tparticle;

static const char USAGE[] =
R"(Hoth - the Monte Carlo game you're looking for!

    http://github.com/mlund/faunus

    Usage:
      faunus [-q] [--nobar] [--input=<file>] [--output=<file>]
      faunus (-h | --help)
      faunus --version

    Options:
      
      -i <file> --input <file>   Input file [default: /dev/stdin].
      -o <file> --output <file>  Output file [default: out.json].
      --nobar                    No progress bar.
      -q --quiet                 Less verbose output.
      -h --help                  Show this screen.
      --version                  Show version.
)";

int main( int argc, char **argv )
{
    Faunus::MPI::MPIController mpi; // OK also if MPI is unavailable

    try {
        auto args = docopt::docopt( USAGE,
                { argv + 1, argv + argc }, true, "Faunus 2.0.0");

        bool showProgress = !args["--nobar"].asBool();

        bool quiet = args["--quiet"].asBool();
        if (quiet) {
            cout.setstate( std::ios_base::failbit ); // hold kæft
            showProgress = false;
        }

        json j;
        auto input = args["--input"].asString();
        if (input=="/dev/stdin")
            std::cin >> j;
        else
            j = openjson(mpi.prefix + input);

        pc::temperature = j.at("temperature").get<double>() * 1.0_K;
        MCSimulation<Tgeometry,Tparticle> sim(j);
        Analysis::CombinedAnalysis analysis(j, sim.space(), sim.pot());

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
        progressBar.done();

        if (!quiet)
            mpi.cout << "relative drift = " << sim.drift() << endl;

        std::ofstream f(mpi.prefix + args["--output"].asString());
        if (f) {
            json j = sim;
            j["analysis"] = analysis;
            if (mpi.nproc()>1)
                j["mpi"] = mpi;
            f << std::setw(4) << j << endl;
        }

    } catch (std::exception &e) {
        std::cerr << e.what() << endl;
        return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
}
