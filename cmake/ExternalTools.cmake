include(ExternalProject)
include(FetchContent)

###############
# CPM Packages
###############

find_package(docopt CONFIG REQUIRED)
find_package(cereal CONFIG REQUIRED)
find_package(doctest CONFIG REQUIRED)
find_package(Eigen3 CONFIG REQUIRED)

CPMAddPackage("gh:gabime/spdlog@1.11.0")
CPMAddPackage("gh:ericniebler/range-v3#0.12.0")
#CPMAddPackage("gh:docopt/docopt.cpp#v0.6.3")
#CPMAddPackage("gh:doctest/doctest#v2.4.9")
CPMAddPackage("gh:mateidavid/zstr#v1.0.6")
CPMAddPackage("gh:martinus/nanobench#v4.3.9")
CPMAddPackage("gh:pybind/pybind11#v2.10.1")
CPMAddPackage("gh:imneme/pcg-cpp#ffd522e7188bef30a00c74dc7eb9de5faff90092")
CPMAddPackage("gh:ArashPartow/exprtk#93a9f44f99b910bfe07cd1e933371e83cea3841c")

CPMAddPackage(
    NAME mpl GITHUB_REPOSITORY rabauke/mpl DOWNLOAD_ONLY YES
    GIT_TAG ff9512fc61195b6c7e643e234789b0b937d28ee3
)

CPMAddPackage(
    NAME nlohmann_json VERSION 3.11.2
    URL https://github.com/nlohmann/json/releases/download/v3.11.2/include.zip
    OPTIONS "JSON_BuildTests OFF"
)

#CPMAddPackage(
#    NAME Eigen VERSION 3.4.0 DOWNLOAD_ONLY YES
#    URL https://gitlab.com/libeigen/eigen/-/archive/3.4.0/eigen-3.4.0.tar.gz
#)

#CPMAddPackage(
#    NAME cereal VERSION 1.3.2 GITHUB_REPOSITORY USCiLab/cereal
#    OPTIONS "SKIP_PORTABILITY_TEST ON" "JUST_INSTALL_CEREAL ON"
#)

CPMAddPackage("gh:pybind/pybind11_json#0.2.13")


###################################
# Configure CPM packages if needed
###################################

set_property(TARGET spdlog PROPERTY POSITION_INDEPENDENT_CODE ON)
#set_property(TARGET docopt PROPERTY POSITION_INDEPENDENT_CODE ON)

if (nlohmann_json_ADDED)
    add_library(nlohmann_json INTERFACE IMPORTED)
    target_include_directories(nlohmann_json INTERFACE ${nlohmann_json_SOURCE_DIR}/include)
endif()

add_compile_definitions("NLOHMANN_JSON_HPP") # older versions used this macro. Now it's suffixed with "_"

if(mpl_ADDED)
    add_library(mpl INTERFACE IMPORTED)
    target_include_directories(mpl INTERFACE ${mpl_SOURCE_DIR})
endif()

#if(Eigen_ADDED)
#    add_library(Eigen INTERFACE IMPORTED)
#    target_include_directories(Eigen INTERFACE ${Eigen_SOURCE_DIR})
#endif()

if (pcg-cpp_ADDED)
    add_library(pcg-cpp INTERFACE)
    target_include_directories(pcg-cpp INTERFACE "${pcg-cpp_SOURCE_DIR}/include")
endif()
option(ENABLE_PCG "Enable PCG random number generator" off)
if (ENABLE_PCG)
    add_definitions("-DENABLE_PCG")
endif()

if (exprtk_ADDED)
    add_definitions("-Dexprtk_disable_string_capabilities")
    add_definitions("-Dexprtk_disable_rtl_io_file")
    add_library(exprtk INTERFACE)
    target_include_directories(exprtk INTERFACE "${exprtk_SOURCE_DIR}")
endif()

###################
# PROGRESS TRACKER
###################

ExternalProject_Add(
    project_progresstracker
    PREFIX "${CMAKE_CURRENT_BINARY_DIR}/_deps"
    CMAKE_ARGS -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER} -DCMAKE_POSITION_INDEPENDENT_CODE=on
    BUILD_COMMAND ${CMAKE_MAKE_PROGRAM} progresstracker
    INSTALL_COMMAND ""
    LOG_DOWNLOAD ON
    DOWNLOAD_EXTRACT_TIMESTAMP true
    URL https://github.com/mlund/progress-cpp/archive/74c33b1eb21417fef9e5fc2b02c7dbe1d533010c.zip
    URL_HASH SHA256=45e2e83a351d44fc1723aecdf1fbf7cee1afc5d44b7190128d8fd6b4437d15b4
)
ExternalProject_Get_Property(project_progresstracker binary_dir)
ExternalProject_Get_Property(project_progresstracker source_dir)
set(ProgressTrackerIncludeDir ${source_dir})
add_library(progresstracker STATIC IMPORTED GLOBAL)
add_dependencies(progresstracker project_progresstracker)
set_property(TARGET progresstracker PROPERTY IMPORTED_LOCATION ${binary_dir}/libprogresstracker.a)


#########
# CPPSID
#########

option(ENABLE_SID "Enable SID emulation" off)
if(ENABLE_SID)
    find_package(SDL2 CONFIG REQUIRED)
    ExternalProject_Add(project_cppsid
        PREFIX "${CMAKE_CURRENT_BINARY_DIR}/_deps"
        CMAKE_ARGS -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER} -DCMAKE_POSITION_INDEPENDENT_CODE=on
        BUILD_COMMAND ${CMAKE_MAKE_PROGRAM} "cppsid"
        INSTALL_COMMAND "" LOG_DOWNLOAD ON
        UPDATE_DISCONNECTED ON
        URL_MD5 b420c4c114e00a147c2c9a974249f0d4
        DOWNLOAD_EXTRACT_TIMESTAMP true
        URL "https://github.com/mlund/cppsid/archive/v0.2.1.tar.gz")
    ExternalProject_Get_Property(project_cppsid binary_dir)
    ExternalProject_Get_Property(project_cppsid source_dir)
    set(CppsidIncludeDir ${source_dir}/include)
    add_library(cppsid STATIC IMPORTED GLOBAL)
    add_dependencies(cppsid project_cppsid)
    set_property(TARGET cppsid PROPERTY IMPORTED_LOCATION ${binary_dir}/libcppsid.a)
endif()

############
# INTEL TBB
############

option(ENABLE_TBB "Enable Intel TBB" off)
if (ENABLE_TBB)
    find_package(TBB REQUIRED COMPONENTS tbb)
    target_link_libraries(project_options INTERFACE TBB::tbb)
endif ()

##########
# XRDFILE
##########

ExternalProject_Add(
    project_xdrfile
    PREFIX "${CMAKE_CURRENT_BINARY_DIR}/_deps"
    URL "https://github.com/chemfiles/xdrfile/archive/8935d749e1f43a87221089588d1cc3f37a0354b0.tar.gz"
    URL_HASH SHA256=a5530703fd07a5baadc9ba75d806fe0844d7b3da0e16f5adbb966660a1cd6828
    PATCH_COMMAND patch -p1 < ${CMAKE_SOURCE_DIR}/cmake/patches/xdrfile-01.patch
    BUILD_COMMAND ${CMAKE_MAKE_PROGRAM} xdrfile-static
    DOWNLOAD_EXTRACT_TIMESTAMP true
    UPDATE_DISCONNECTED ON
    CMAKE_ARGS -Wno-dev -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER} -DCMAKE_POSITION_INDEPENDENT_CODE=on
    LOG_DOWNLOAD ON INSTALL_COMMAND "")

ExternalProject_Get_Property(project_xdrfile source_dir)
ExternalProject_Get_Property(project_xdrfile binary_dir)
set(XdrfileIncludeDir ${source_dir}/include)
add_library(xdrfile STATIC IMPORTED)
set_property(TARGET xdrfile PROPERTY IMPORTED_LOCATION ${binary_dir}/libxdrfile-static.a)
add_dependencies(xdrfile project_xdrfile)
set_target_properties(xdrfile PROPERTIES POSITION_INDEPENDENT_CODE TRUE)


##############
# TROMPELOEIL
##############

FetchContent_Declare(
    trompeloeil
    URL "https://github.com/rollbear/trompeloeil/archive/v41.tar.gz"
    DOWNLOAD_EXTRACT_TIMESTAMP true
    URL_HASH SHA256=48986b507497f027e4fa1144a08c2d0b6d81fb476fad024956f8104448ca9ad8 DOWNLOAD_EXTRACT_TIMESTAMP true)
FetchContent_GetProperties(trompeloeil)
if(NOT trompeloeil_POPULATED)
    FetchContent_Populate(trompeloeil)
endif()

###########
# FREESASA
###########

option(ENABLE_FREESASA "Fetch 3rd-party SASA calculation software" on)
if (ENABLE_FREESASA)
    ExternalProject_Add(
            project_freesasa
            PREFIX "${CMAKE_CURRENT_BINARY_DIR}/_deps"
            LOG_DOWNLOAD ON
            DOWNLOAD_EXTRACT_TIMESTAMP true
            URL https://github.com/mittinatten/freesasa/releases/download/2.0.3/freesasa-2.0.3.tar.gz
            URL_HASH SHA256=ba1d4f7e9dd51ae2452b5c3a80ac34039d51da4826dae1dbe173cd7a1d6aca94
            # -fPIC flag is needed to link with pyfaunus
            CONFIGURE_COMMAND CFLAGS=-fPIC <SOURCE_DIR>/configure --disable-xml --disable-json
            BUILD_COMMAND ${CMAKE_MAKE_PROGRAM}
            INSTALL_COMMAND ""
    )
    ExternalProject_Get_Property(project_freesasa source_dir)
    ExternalProject_Get_Property(project_freesasa binary_dir)
    add_library(freesasa STATIC IMPORTED)
    set_property(TARGET freesasa PROPERTY IMPORTED_LOCATION ${binary_dir}/src/libfreesasa.a)
    add_definitions("-DENABLE_FREESASA")
    include_directories(SYSTEM "${source_dir}/src")
    add_dependencies(freesasa project_freesasa)
endif ()

################
# COULOMBGALORE
################

FetchContent_Declare(
    coulombgalore
    URL https://github.com/mlund/coulombgalore/archive/4055f58538d781acccb2937ab4580855fcba31f8.tar.gz
    URL_HASH MD5=922f0c5988c0f70c887d65b7cf2762ac
    DOWNLOAD_EXTRACT_TIMESTAMP true)
FetchContent_GetProperties(coulombgalore)
if(NOT coulombgalore_POPULATED)
    FetchContent_Populate(coulombgalore)
endif()
include_directories(SYSTEM ${coulombgalore_SOURCE_DIR})

# Add third-party headers to include path. Note this is done with SYSTEM
# to disable potential compiler warnings

include_directories(SYSTEM ${trompeloeil_SOURCE_DIR}/include ${nanobench_SOURCE_DIR}/src/include
    ${Pybind11IncludeDir} ${CppsidIncludeDir} ${XdrfileIncludeDir} ${ProgressTrackerIncludeDir})
