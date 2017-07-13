cmake_minimum_required(VERSION 3.8)
project(neofaunus)

set (CMAKE_CXX_STANDARD 14)

include(ExternalProject)

ExternalProject_Add( cpp11range
    GIT_REPOSITORY "https://github.com/klmr/cpp11-range.git"
    GIT_TAG "master"
    CONFIGURE_COMMAND "" BUILD_COMMAND "" INSTALL_COMMAND "" )
ExternalProject_Get_Property(cpp11range source_dir)
set(Cpp11RangeIncludeDir ${source_dir})

ExternalProject_Add( modernjson
    GIT_REPOSITORY "https://github.com/nlohmann/json.git"
    GIT_TAG "master"
    CONFIGURE_COMMAND "" BUILD_COMMAND "" INSTALL_COMMAND "" )
ExternalProject_Get_Property(modernjson source_dir)
set(JSONIncludeDir ${source_dir}/src)

ExternalProject_Add( range
    HG_REPOSITORY "https://bitbucket.org/AraK/range"
    CONFIGURE_COMMAND "" BUILD_COMMAND "" INSTALL_COMMAND "" )
ExternalProject_Get_Property(range source_dir)
set(RangeIncludeDir ${source_dir})


ExternalProject_Add( eigen
    HG_REPOSITORY "https://bitbucket.org/eigen/eigen"
    HG_TAG "3.3.4"
    CONFIGURE_COMMAND "" BUILD_COMMAND "" INSTALL_COMMAND "" )
ExternalProject_Get_Property(eigen source_dir)
set(EigenIncludeDir ${source_dir})


message(${RangeIncludeDir})
message(${Cpp11RangeIncludeDir})
message(${EigenIncludeDir})
message(${JSONIncludeDir})

add_executable(point point.cpp)
target_include_directories(point PRIVATE ${EigenIncludeDir} ${RangeIncludeDir} ${JSONIncludeDir})


