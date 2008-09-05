cmake_minimum_required(VERSION 2.6)
project ("faunus")
# ------------------------------
#   Default build is "Release"
# ------------------------------
IF(NOT CMAKE_BUILD_TYPE)
  SET(CMAKE_BUILD_TYPE Release CACHE STRING
    "Build type: None Debug Release RelWithDebInfo MinSizeRel."
    FORCE)
ENDIF(NOT CMAKE_BUILD_TYPE)

# -----------------
#   User options
# -----------------
option (OPENMP "Enable OpenMP parallel execution on multicore CPU's")
option (GROMACS "Enable Gromacs xtc trajectory output" off)
mark_as_advanced(CLEAR CMAKE_CXX_COMPILER CMAKE_CXX_FLAGS_RELEASE)
mark_as_advanced(EXECUTABLE_OUTPUT_PATH LIBRARY_OUTPUT_PATH CMAKE_OSX_SYSROOT)

# ----------------------------
#   Locate external packages
# ----------------------------
set (CMAKE_MODULE_PATH cmake/modules)
include(TestCXXAcceptsFlag)
find_package(OpenBabel2)
if (GROMACS)
  find_package(GMX)
  if (GMX_MOTIF)
    find_package(X11)
    find_package(Motif)
  endif (GMX_MOTIF)
endif(GROMACS)

# ------------------
#   Compiler flags
# ------------------
function( fau_setcmp cid flags fomp )
if (CMAKE_CXX_COMPILER_ID MATCHES "${cid}")
  set (DUMMY "${flags}" PARENT_SCOPE)
  CHECK_CXX_ACCEPTS_FLAG( ${fomp} OMP_CHECK)
  if (OPENMP AND OMP_CHECK)
    set(DUMMY "${flags} ${fomp}" PARENT_SCOPE)
  endif(OPENMP AND OMP_CHECK)
endif (CMAKE_CXX_COMPILER_ID MATCHES "${cid}")
set (CMAKE_CXX_FLAGS_RELEASE "${DUMMY}" CACHE STRING "C++ compiler flags (automatic)" FORCE)
endfunction( fau_setcmp )

fau_setcmp("GNU" "-O3 -w -funroll-loops -DNDEBUG" "-fopenmp" ) # GNU
fau_setcmp("Intel" "-O3 -w -DNDEBUG" "-openmp")  # Intel
fau_setcmp("Sun" "-fast -w -DNDEBUG" "-xopenmp") # Sun
fau_setcmp("pgCC" "-O3 -w -DNDEBUG" "-mp")       # Portland
fau_setcmp("pathCC" "-Ofast -DNDEBUG" " ")       # Pathscale

# ---------------------------
#   Add source directories
# ---------------------------
add_subdirectory (classes)
add_subdirectory (examples)
add_subdirectory (legacy/fortran)
add_subdirectory (doc EXCLUDE_FROM_ALL)

