cmake_minimum_required(VERSION 3.8.2 FATAL_ERROR)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

project(neofaunus)

include(ExternalProject)

ExternalProject_Add( cpp11range
    PREFIX "${CMAKE_SOURCE_DIR}/3rdparty"
    GIT_REPOSITORY "https://github.com/klmr/cpp11-range.git"
    GIT_TAG "master"
    LOG_DOWNLOAD ON CONFIGURE_COMMAND "" BUILD_COMMAND "" INSTALL_COMMAND "" )
ExternalProject_Get_Property(cpp11range source_dir)
set(Cpp11RangeIncludeDir ${source_dir})

ExternalProject_Add( modernjson
    PREFIX "${CMAKE_SOURCE_DIR}/3rdparty"
    URL "https://github.com/nlohmann/json/releases/download/v2.1.1/json.hpp"
    DOWNLOAD_NO_EXTRACT 0 DOWNLOAD_NO_PROGRESS ON
    CONFIGURE_COMMAND "" BUILD_COMMAND "" INSTALL_COMMAND "" )
ExternalProject_Get_Property(modernjson source_dir)
set(JSONIncludeDir ${source_dir}/../../src/)
message(${JSONIncludeDir})

ExternalProject_Add( range
    PREFIX "${CMAKE_SOURCE_DIR}/3rdparty"
    URL "https://bitbucket.org/AraK/range/get/tip.tar.gz"
    URL_MD5 622daaa483f540aced28c6f7434ff083
    LOG_DOWNLOAD ON CONFIGURE_COMMAND "" BUILD_COMMAND "" INSTALL_COMMAND "" )
ExternalProject_Get_Property(range source_dir)
set(RangeIncludeDir ${source_dir})

ExternalProject_Add( eigen
    PREFIX "${CMAKE_SOURCE_DIR}/3rdparty"
    URL "http://bitbucket.org/eigen/eigen/get/3.3.4.tar.gz"
    URL_MD5 1a47e78efe365a97de0c022d127607c3
    LOG_DOWNLOAD ON CONFIGURE_COMMAND "" BUILD_COMMAND "" INSTALL_COMMAND "" )
ExternalProject_Get_Property(eigen source_dir)
set(EigenIncludeDir ${source_dir})

add_executable(point src/point.cpp)
add_dependencies(point modernjson range eigen)
target_include_directories(point PRIVATE
    ${EigenIncludeDir} ${RangeIncludeDir} ${JSONIncludeDir})

