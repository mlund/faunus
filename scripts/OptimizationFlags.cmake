# Compiler specific optimization flags
unset(CMAKE_CXX_FLAGS)

# GNU
if (CMAKE_CXX_COMPILER_ID MATCHES "GNU")
  set(CMAKE_CXX_FLAGS "-std=c++11 -funroll-loops -Wall -Wno-unknown-pragmas -Wextra -Wno-unused-parameter -Wno-reorder")

# Intel
# (list of warning codes: icpc -diag-dump a.cpp)
elseif (CMAKE_CXX_COMPILER_ID MATCHES "Intel")
  set(CMAKE_CXX_FLAGS "-std=c++11 -Wall -Wcheck -wd2259,3180 -Wno-unknown-pragmas")

# Clang
elseif (CMAKE_CXX_COMPILER_ID MATCHES "Clang")
  set(CMAKE_CXX_FLAGS "-Wextra -pedantic -std=c++11 -Wno-unused-parameter -Wno-unknown-pragmas")
  if(APPLE)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -stdlib=libc++")
  endif()

# The following are untested -- C++11 support unknown
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
