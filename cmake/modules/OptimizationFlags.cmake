# Compiler specific optimization flags.
# - Set FAUNUS_OPENMP to test and enable OpenMP support
# $Mikael Lund, 2008

include(TestCXXAcceptsFlag)

function( fau_setcmp cxxid flags fomp )
if (NOT SKIP_COMPILER_CHECK)
  if (CMAKE_CXX_COMPILER_ID MATCHES "${cxxid}")
    set (newflags "${flags}")
    CHECK_CXX_ACCEPTS_FLAG( ${fomp} OMP_CHECK)
    if (FAUNUS_OPENMP AND OMP_CHECK)
      set(newflags "${newflags} ${fomp}")
    endif(FAUNUS_OPENMP AND OMP_CHECK)
    set (CMAKE_CXX_FLAGS_RELEASE "${newflags}" CACHE STRING "C++ compiler flags (automatic)" FORCE)
    set (SKIP_COMPILER_CHECK ON CACHE INTERNAL "Skip Compiler Test")
  endif (CMAKE_CXX_COMPILER_ID MATCHES "${cxxid}")
endif (NOT SKIP_COMPILER_CHECK)
endfunction( fau_setcmp )

fau_setcmp("GNU" "-O3 -funroll-loops -W -Wno-sign-compare -Wconversion -DNDEBUG" "-fopenmp" ) # GNU
fau_setcmp("Intel" "-O3 -w1 -Wno-unknown-pragmas -DNDEBUG" "-openmp") # Intel
fau_setcmp("Sun" "-fast -DNDEBUG" "-xopenmp") # Sun
fau_setcmp("pgCC" "-O3 -DNDEBUG" "-mp")       # Portland
fau_setcmp("pathCC" "-Ofast -DNDEBUG" " ")    # Pathscale

