# Compiler specific optimization flags.
# - Set FAUNUS_OPENMP to test and enable OpenMP support
# $Mikael Lund, 2008

function( fau_setcmp cid flags fomp )
if (CMAKE_CXX_COMPILER_ID MATCHES "${cid}")
  set (DUMMY "${flags}" PARENT_SCOPE)
  CHECK_CXX_ACCEPTS_FLAG( ${fomp} OMP_CHECK)
  if (FAUNUS_OPENMP AND OMP_CHECK)
    set(DUMMY "${flags} ${fomp}" PARENT_SCOPE)
  endif(FAUNUS_OPENMP AND OMP_CHECK)
endif (CMAKE_CXX_COMPILER_ID MATCHES "${cid}")
set (CMAKE_CXX_FLAGS_RELEASE "${DUMMY}" CACHE STRING "C++ compiler flags (automatic)" FORCE)
endfunction( fau_setcmp )

fau_setcmp("GNU" "-O3 -w -funroll-loops -DNDEBUG" "-fopenmp" ) # GNU
fau_setcmp("Intel" "-O3 -w -DNDEBUG" "-openmp")  # Intel
fau_setcmp("Sun" "-fast -w -DNDEBUG" "-xopenmp") # Sun
fau_setcmp("pgCC" "-O3 -w -DNDEBUG" "-mp")       # Portland
fau_setcmp("pathCC" "-Ofast -DNDEBUG" " ")       # Pathscale

