# Compiler specific optimization flags.
# - Set ENABLE_OPENMP to enable OpenMP support
# $Mikael Lund, 2008
# See http://fedetft.wordpress.com/2009/12/21/cmake-part-2-compiler-flags/
#set(CMAKE_CXX_FLAGS "-pedantic -Wall -stdlib=libc++ -Wno-long-long -Wno-unknown-pragmas -Wc++0x-extensions")
unset(CMAKE_CXX_FLAGS)

if (CMAKE_CXX_COMPILER_ID MATCHES "GNU")
  set(CMAKE_CXX_FLAGS "-std=c++0x -funroll-loops -Wall -Wno-unknown-pragmas")
elseif (CMAKE_CXX_COMPILER_ID MATCHES "Intel")
  set(CMAKE_CXX_FLAGS "-std=c++0x -w1 -Wno-unknown-pragmas")
elseif (CMAKE_CXX_COMPILER_ID MATCHES "Clang")
  set(CMAKE_CXX_FLAGS "-W -pedantic -stdlib=libc++ -Wno-long-long -Wno-unused-parameter -Wno-unknown-pragmas -Wno-c++0x-extensions")
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
