#undef DOCTEST_CONFIG_DISABLE
#define DOCTEST_CONFIG_IMPLEMENT
#define DOCTEST_CONFIG_SUPER_FAST_ASSERTS

#include <doctest/doctest.h>
#include <spdlog/spdlog.h>
#include <spdlog/sinks/null_sink.h>
#include <spdlog/sinks/basic_file_sink.h>
#include <spdlog/sinks/stdout_color_sinks.h>

#include "aux/eigensupport.h"
#include "atomdata_test.h"
#include "random.h"
#include "core_test.h"
#include "tensor_test.h"
#include "mpicontroller.h"
#include "auxiliary.h"
#include "molecule.h"
#include "group_test.h"
#include "geometry.h"
#include "regions.h"
#include "space.h"
#include "average.h"
#include "potentials_test.h"
#include "tabulate.h"
#include "move.h"
#include "penalty.h"
#include "celllist.h"
#include "functionparser.h"
#include "multipole.h"
#include "energy_test.h"

int main(int argc, char** argv) {
    Faunus::faunus_logger = spdlog::basic_logger_mt("faunus", "unittests.log", true);
    Faunus::faunus_logger->set_pattern("%L: %v");
    Faunus::faunus_logger->set_level(spdlog::level::debug);

    doctest::Context context;
    context.applyCommandLine(argc, argv);
    int res = context.run();

    if(context.shouldExit()) // important - query flags (and --exit) rely on the user doing this
        return res;          // propagate the result of the tests
    return res;              // the result from doctest is propagated here as well
}
