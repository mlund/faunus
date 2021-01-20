cmake_minimum_required(VERSION 3.16 FATAL_ERROR)

project(faunus)
include(cmake/StandardProjectSettings.cmake)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED YES)
set(CMAKE_CXX_EXTENSIONS YES)

if(EXISTS ${CMAKE_BINARY_DIR}/conanbuildinfo.cmake)
   include(${CMAKE_BINARY_DIR}/conanbuildinfo.cmake)
   conan_basic_setup(TARGETS)
else()
    message(WARNING "First run `conan install -b missing -if ${CMAKE_BINARY_DIR}`")
endif()

# Interface library for all Faunus targets
add_library(project_options INTERFACE)
target_compile_features(project_options INTERFACE cxx_std_17)

# Link this 'library' to use the warnings specified in CompilerWarnings.cmake
add_library(project_warnings INTERFACE)

# Enable cache system
include(cmake/Cache.cmake)

include(cmake/CompilerWarnings.cmake)
set_project_warnings(project_warnings)

set(THREADS_PREFER_PTHREAD_FLAG ON)
find_package(PythonInterp 3 REQUIRED)
find_package(Threads REQUIRED)
find_package(ZLIB REQUIRED)
enable_testing()

include(cmake/CompilerFlags.cmake)
include(cmake/ExternalTools.cmake)
include(cmake/GITversion.cmake)

add_subdirectory(src)
add_subdirectory(docs)
add_subdirectory(examples)

