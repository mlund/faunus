#define DOCTEST_CONFIG_IMPLEMENT
#include <doctest/doctest.h>
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
#include <range/v3/view/common.hpp>

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
    R"(PDB Scatter

    https://faunus.readthedocs.io

    Usage:
      pdb_scatter [-q] [--multiplications <N>] [--verbosity <N>] [--output=<file>] <files>...
      pdb_scatter (-h | --help)
      pdb_scatter --version

    Options:
      -o <file> --output <file>         Output file [default: scatter.dat].
      -p <N> --multiplications          Determines upper q-value [default: 15]
      -v <N> --verbosity <N>            Log verbosity level (0 = off, 1 = critical, ..., 6 = trace) [default: 4]
      -q --quiet                        Less verbose output. It implicates -v0 --nobar --notips --nofun.
      -h --help                         Show this screen.
      --version                         Show version.
)";

using ProgressIndicator::ProgressTracker;

// forward declarations
std::shared_ptr<ProgressTracker> createProgressTracker(bool, unsigned int);

int main(int argc, const char **argv) {
    bool quiet = false;
    try {

        auto starting_time = std::chrono::steady_clock::now(); // used to time the simulation
        std::string version = "PDB scatter";
        version += "0.1";
        auto args = docopt::docopt(USAGE, {argv + 1, argv + argc}, true, version);

        // --quiet; set quiet first as it also affects exception handling
        quiet = args["--quiet"].asBool();
        if (quiet) {
            cout.setstate(std::ios_base::failbit); // hold kæft
        }

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

        // Hard coded atom types required for PQR loader
        Faunus::atoms = R"([
             { "Na": { "sigma": 3.166, "eps": 0.65, "mw": 1.0 } },
             { "Cl": { "sigma": 2.0, "eps": 0.0, "mw": 1.0 } }
             ])"_json.get<decltype(atoms)>();

        // Scattering analysis of multiple input PQR files
        const auto output_file = args["--output"].asString();
        const auto multiplications = args["--multiplications"].asLong();
        faunus_logger->info("q multiplications set to {}", multiplications);
        Scatter::StructureFactorPBC scatter(multiplications);
        ParticleVector particles;
        for (const auto filename : args["<files>"].asStringList()) {
            faunus_logger->info("analysing {s}", filename);
            particles.clear();
            auto box_dimensions = FormatPQR::load(filename, particles, true);
            auto positions = particles | ranges::cpp20::views::transform([&](const auto &p) { return p.pos; });
            scatter.sample(positions, box_dimensions);
        }
        IO::write(output_file, scatter.getSampling());

    } catch (std::exception &e) {
        faunus_logger->error(e.what());
        if (!usageTip.buffer.empty()) {
            faunus_logger->error(usageTip.buffer);
        }
        return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
}