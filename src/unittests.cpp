#undef DOCTEST_CONFIG_DISABLE

#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#define DOCTEST_CONFIG_SUPER_FAST_ASSERTS

#include <doctest/doctest.h>
#include "random.h"
#include "core.h"
#include "mpi.h"
#include "auxiliary.h"
#include "molecule.h"
#include "group.h"
#include "geometry.h"
#include "space.h"
#include "average.h"
#include "tabulate.h"
#include "move.h"
#include "penalty.h"
#include "celllist.h"

