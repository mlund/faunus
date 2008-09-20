# Compiler specific optimization flags.
# - Set FAUNUS_OPENMP to test and enable OpenMP support
# $Mikael Lund, 2008

function( fau_setcmp cxxid flags fomp )
if (CMAKE_CXX_COMPILER_ID MATCHES "${cxxid}")
  set (newflags "${flags}")
  CHECK_CXX_ACCEPTS_FLAG( ${fomp} OMP_CHECK)
  if (FAUNUS_OPENMP AND OMP_CHECK)
    set(newflags "${newflags} ${fomp}")
  endif(FAUNUS_OPENMP AND OMP_CHECK)
  set (CMAKE_CXX_FLAGS_RELEASE "${newflags}" CACHE STRING "C++ compiler flags (automatic)" FORCE)
endif (CMAKE_CXX_COMPILER_ID MATCHES "${cxxid}")
endfunction( fau_setcmp )

fau_setcmp("GNU" "-O3 -w -funroll-loops -DNDEBUG" "-fopenmp" ) # GNU
fau_setcmp("Intel" "-O3 -w -DNDEBUG" "-openmp")  # Intel
fau_setcmp("Sun" "-fast -w -DNDEBUG" "-xopenmp") # Sun
fau_setcmp("pgCC" "-O3 -w -DNDEBUG" "-mp")       # Portland
fau_setcmp("pathCC" "-Ofast -DNDEBUG" " ")       # Pathscale

