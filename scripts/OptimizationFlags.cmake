# Compiler specific optimization flags
unset(CMAKE_CXX_FLAGS)

# GNU
if (CMAKE_CXX_COMPILER_ID MATCHES "GNU")
  set(CMAKE_CXX_FLAGS "-funroll-loops -Wall -Wno-unknown-pragmas -Wextra -Wno-unused-parameter -Wno-reorder -Wno-misleading-indentation")

# Intel
elseif (CMAKE_CXX_COMPILER_ID MATCHES "Intel")
  set(CMAKE_CXX_FLAGS "-Wall -Wcheck -wd2259,3180 -Wno-unknown-pragmas")

# Clang
elseif (CMAKE_CXX_COMPILER_ID MATCHES "Clang")
  set(CMAKE_CXX_FLAGS "-Wextra -pedantic -Wno-unused-parameter -Wno-unknown-pragmas")
  if(APPLE)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -stdlib=libc++")
  endif()
endif()

if (ENABLE_OPENMP)
  find_package(OpenMP)
  if (OPENMP_FOUND)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
  endif()
endif()

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}" CACHE STRING "Additional compilation flags" FORCE)
