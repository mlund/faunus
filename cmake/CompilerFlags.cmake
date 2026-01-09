# in Debug mode, all warnings are treated as errors and we want no optimisations
#add_compile_options($<$<CONFIG:Debug>:-Werror>)
add_compile_options($<$<CONFIG:Debug>:-O0>)

# in Release mode, add aggressive optimizations
add_compile_options($<$<CONFIG:Release>:-O3>)
add_compile_options($<$<CONFIG:Release>:-ffast-math>)
add_compile_options($<$<CONFIG:Release>:-fno-finite-math-only>)

if (CMAKE_SYSTEM_PROCESSOR MATCHES "(x86)|(X86)|(amd64)|(AMD64)")
    add_compile_options($<$<CONFIG:Release>:-march=native>)
endif()

if (APPLE AND CMAKE_SYSTEM_PROCESSOR MATCHES "arm64")
    if (CMAKE_CXX_COMPILER_ID MATCHES "GNU")
        add_compile_options($<$<CONFIG:Release>:-mcpu=native>)
    endif()
    if (CMAKE_CXX_COMPILER_ID MATCHES "Clang")
        add_compile_options($<$<CONFIG:Release>:-mcpu=apple-a14>)
    endif()
endif()


option(ENABLE_APPROXMATH "Use approximate math" off)
if (ENABLE_APPROXMATH)
    target_compile_definitions(project_options INTERFACE FAU_APPROXMATH)
endif ()

option(ENABLE_OPENMP "Try to use OpenMP parallelisation" on)
if (ENABLE_OPENMP)
    find_package(OpenMP)
    if (OpenMP_CXX_FOUND)
        target_link_libraries(project_options INTERFACE OpenMP::OpenMP_CXX)
    endif ()
endif ()
