#  Compiler specific flags
## GCC
if (CMAKE_CXX_COMPILER_ID MATCHES "GNU")
    add_compile_options(-Wall -Wextra -Wpedantic -Wunreachable-code -Wstrict-aliasing
        -Wno-sign-compare -Wno-unused-local-typedefs -Wno-unknown-pragmas)
## Clang
elseif (CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    add_compile_options(-Wall -Wextra -Wpedantic -Wunreachable-code -fsized-deallocation
        -Wstrict-aliasing -Wno-sign-compare -Wno-unused-local-typedef -Wno-unknown-pragmas)
endif()

# in Debug mode, all warnings are treated as errors and we want no optimisations
#add_compile_options($<$<CONFIG:Debug>:-Werror>)
add_compile_options($<$<CONFIG:Debug>:-O0>)

# in Release mode, add aggressive optimizations
add_compile_options($<$<CONFIG:Release>:-march=native>)
add_compile_options($<$<CONFIG:Release>:-Ofast>)

option(ENABLE_APPROXMATH "Use approximate math" off)
if (ENABLE_APPROXMATH)
    add_definitions(-DFAU_APPROXMATH)
endif ()


