# in Debug mode, all warnings are treated as errors and we want no optimisations
#add_compile_options($<$<CONFIG:Debug>:-Werror>)
add_compile_options($<$<CONFIG:Debug>:-O0>)

# in Release mode, add aggressive optimizations
add_compile_options($<$<CONFIG:Release>:-Ofast>)
add_compile_options($<$<CONFIG:Release>:-tbb>) # icpx specific
add_compile_options($<$<CONFIG:RelWithDebInfo>:-tbb>) # icpx specific
add_compile_options($<$<CONFIG:Release>:-g>)
#add_compile_options($<$<CONFIG:RelWithDebInfo>:-fiopenmp>)
#add_compile_options($<$<CONFIG:Release>:-fiopenmp>)
add_compile_options($<$<CONFIG:RelWithDebInfo>:-fno-finite-math-only>)
add_compile_options($<$<CONFIG:Release>:-fno-finite-math-only>)

#if (CMAKE_SYSTEM_PROCESSOR MATCHES "(x86)|(X86)|(amd64)|(AMD64)")
    add_compile_options($<$<CONFIG:Release>:-march=native>)
#endif()

option(ENABLE_APPROXMATH "Use approximate math" off)
if (ENABLE_APPROXMATH)
    add_definitions(-DFAU_APPROXMATH)
endif ()

option(ENABLE_OPENMP "Try to use OpenMP parallisation" on)
if (ENABLE_OPENMP)
    find_package(OpenMP)
    if (OPENMP_FOUND)
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    endif ()
endif ()
