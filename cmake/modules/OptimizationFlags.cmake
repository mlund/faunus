# Compiler specific optimization flags.
# - Set ENABLE_OPENMP to enable OpenMP support
# $Mikael Lund, 2008
# See http://fedetft.wordpress.com/2009/12/21/cmake-part-2-compiler-flags/

unset(CMAKE_CXX_FLAGS)

if (CMAKE_CXX_COMPILER_ID MATCHES "GNU")
  set(CMAKE_CXX_FLAGS "-funroll-loops -W -Wno-sign-compare -Wno-conversion")
elseif (CMAKE_CXX_COMPILER_ID MATCHES "Intel")
  set(CMAKE_CXX_FLAGS "-w1 -Wno-unknown-pragmas")
elseif (CMAKE_CXX_COMPILER_ID MATCHES "CLang")
  set(CMAKE_CXX_FLAGS "-Wno-unused-parameter")
elseif (CMAKE_CXX_COMPILER_ID MATCHES "Sun")
  set(CMAKE_CXX_FLAGS "")
elseif (CMAKE_CXX_COMPILER_ID MATCHES "PGI")
  set(CMAKE_CXX_FLAGS "")
elseif (CMAKE_CXX_COMPILER_ID MATCHES "pathCC")
  set(CMAKE_CXX_FLAGS "")
endif()

if (ENABLE_OPENMP AND OPENMP_FOUND)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}" CACHE STRING "Additional compilation flags" FORCE)
