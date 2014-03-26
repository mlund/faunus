# Compiler specific optimization flags.
# - Set ENABLE_OPENMP to enable OpenMP support
# $Mikael Lund, 2008
# See http://fedetft.wordpress.com/2009/12/21/cmake-part-2-compiler-flags/
unset(CMAKE_CXX_FLAGS)

if (CMAKE_CXX_COMPILER_ID MATCHES "GNU")
  set(CMAKE_CXX_FLAGS "-std=c++0x -funroll-loops -Wall -Wno-unknown-pragmas -Wextra -Wno-unused-parameter -Wno-unused-local-typedefs")
elseif (CMAKE_CXX_COMPILER_ID MATCHES "Intel")
  set(CMAKE_CXX_FLAGS "-std=c++11 -Wall -Wcheck -wd2259,981,869,383 -Wno-unknown-pragmas")
elseif (CMAKE_CXX_COMPILER_ID MATCHES "Clang")
  set(CMAKE_CXX_FLAGS "-Wextra -pedantic -std=c++11 -Wno-long-long -Wno-unused-parameter -Wno-unknown-pragmas -Wno-c++0x-extensions -Wno-variadic-macros")
elseif (CMAKE_CXX_COMPILER_ID MATCHES "Sun")
  set(CMAKE_CXX_FLAGS "-erroff=wvarhidemem,hidevf,badstring,largeshift -errtags=yes")
elseif (CMAKE_CXX_COMPILER_ID MATCHES "PGI")
  set(CMAKE_CXX_FLAGS "")
elseif (CMAKE_CXX_COMPILER_ID MATCHES "pathCC")
  set(CMAKE_CXX_FLAGS "")
endif()

if (ENABLE_OPENMP)
  find_package(OpenMP)
  if (OPENMP_FOUND)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
  endif()
endif()

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}" CACHE STRING "Additional compilation flags" FORCE)
