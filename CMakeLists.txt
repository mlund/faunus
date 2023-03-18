cmake_minimum_required(VERSION 3.24)

project(faunus)
include(cmake/StandardProjectSettings.cmake)
include(cmake/CPM.cmake)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED YES)
set(CMAKE_CXX_EXTENSIONS YES)

# Interface library for all Faunus targets
add_library(project_options INTERFACE)
target_compile_features(project_options INTERFACE cxx_std_20)

# Link this 'library' to use the warnings specified in CompilerWarnings.cmake
add_library(project_warnings INTERFACE)

# Enable cache system
include(cmake/Cache.cmake)

include(cmake/CompilerWarnings.cmake)
set_project_warnings(project_warnings)

set(THREADS_PREFER_PTHREAD_FLAG ON)
find_package(Python 3 COMPONENTS Interpreter REQUIRED)
find_package(Threads REQUIRED)
find_package(ZLIB REQUIRED)
enable_testing()

include(cmake/CompilerFlags.cmake)
include(cmake/ExternalTools.cmake)
include(cmake/GITversion.cmake)

add_subdirectory(src)
add_subdirectory(docs)
add_subdirectory(examples)

