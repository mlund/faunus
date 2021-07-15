include(ExternalProject)
include(FetchContent)

###############
# CPM Packages
###############

CPMAddPackage("gh:gabime/spdlog@1.8.5")
CPMAddPackage("gh:ericniebler/range-v3#0.11.0")
CPMAddPackage("gh:onqtam/doctest#2.4.6")
CPMAddPackage("gh:nlohmann/json@3.9.1")
add_compile_definitions("NLOHMANN_JSON_HPP") # older versions used this macro. Now it's suffixed with "_"

CPMAddPackage(
    NAME cereal VERSION 1.3.0 GITHUB_REPOSITORY USCiLab/cereal
    OPTIONS "SKIP_PORTABILITY_TEST ON" "JUST_INSTALL_CEREAL ON"
)


#############
# PCG RANDOM
#############

option(ENABLE_PCG "Enable PCG random number generator" off)
if (ENABLE_PCG)
    FetchContent_Declare(
        pcg-cpp
        URL https://github.com/imneme/pcg-cpp/archive/ffd522e7188bef30a00c74dc7eb9de5faff90092.tar.gz
        URL_HASH MD5=051b969bbaf924f35f2159813f93e341)
    FetchContent_GetProperties(pcg-cpp)
    if(NOT pcg-cpp_POPULATED)
        FetchContent_Populate(pcg-cpp)
    endif()
    include_directories(SYSTEM ${pcg-cpp_SOURCE_DIR}/include)
    add_definitions("-DENABLE_PCG")
endif()

#########
# EXPRTK 
#########

FetchContent_Declare(
    exprtk
    URL https://github.com/ArashPartow/exprtk/archive/e0e880c3797ea363d24782ba63fe362f7d94f89c.zip
    URL_HASH MD5=772293e80f8353961fcc8a2b337e8dec)
FetchContent_GetProperties(exprtk)
if(NOT exprtk_POPULATED)
    FetchContent_Populate(exprtk)
endif()
add_definitions("-Dexprtk_disable_string_capabilities")
add_definitions("-Dexprtk_disable_rtl_io_file")

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
    URL https://github.com/mlund/progress-cpp/archive/74c33b1eb21417fef9e5fc2b02c7dbe1d533010c.zip
    URL_HASH SHA256=45e2e83a351d44fc1723aecdf1fbf7cee1afc5d44b7190128d8fd6b4437d15b4
)
ExternalProject_Get_Property(project_progresstracker binary_dir)
ExternalProject_Get_Property(project_progresstracker source_dir)
set(ProgressTrackerIncludeDir ${source_dir})
add_library(progresstracker STATIC IMPORTED GLOBAL)
add_dependencies(progresstracker project_progresstracker)
set_property(TARGET progresstracker PROPERTY IMPORTED_LOCATION ${binary_dir}/libprogresstracker.a)

#######
# ZSTR
#######

FetchContent_Declare(
    zstr
    URL "https://github.com/mateidavid/zstr/archive/v1.0.1.tar.gz"
    URL_HASH MD5=42de51b1c6adac0ec957a24088ef7523)
FetchContent_GetProperties(zstr)
if(NOT zstr_POPULATED)
    FetchContent_Populate(zstr)
endif()

#############
# DOCOPT.CPP
#############

ExternalProject_Add(project_docopt
    PREFIX "${CMAKE_CURRENT_BINARY_DIR}/_deps"
    CMAKE_ARGS -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER} -DCMAKE_POSITION_INDEPENDENT_CODE=on
    BUILD_COMMAND ${CMAKE_MAKE_PROGRAM} docopt_s
    INSTALL_COMMAND ""
    LOG_DOWNLOAD ON
    UPDATE_DISCONNECTED ON
    URL_MD5 c6290672c8dae49a01774297a51046fe
    URL "https://github.com/docopt/docopt.cpp/archive/v0.6.3.tar.gz")

ExternalProject_Get_Property(project_docopt binary_dir)
ExternalProject_Get_Property(project_docopt source_dir)
set(DocoptIncludeDir ${source_dir})
add_library(docopt STATIC IMPORTED GLOBAL)
add_dependencies(docopt project_docopt)
set_property(TARGET docopt PROPERTY IMPORTED_LOCATION ${binary_dir}/libdocopt.a)

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
        URL "https://github.com/mlund/cppsid/archive/v0.2.1.tar.gz")
    ExternalProject_Get_Property(project_cppsid binary_dir)
    ExternalProject_Get_Property(project_cppsid source_dir)
    set(CppsidIncludeDir ${source_dir}/include)
    add_library(cppsid STATIC IMPORTED GLOBAL)
    add_dependencies(cppsid project_cppsid)
    set_property(TARGET cppsid PROPERTY IMPORTED_LOCATION ${binary_dir}/libcppsid.a)
endif()

###########
# PYBIND11
###########

FetchContent_Declare(
    pybind11
    URL https://github.com/pybind/pybind11/archive/v2.5.0.tar.gz
    URL_HASH MD5=1ad2c611378fb440e8550a7eb6b31b89)

############
# INTEL TBB
############

option(ENABLE_TBB "Enable Intel TBB" off)
if (ENABLE_TBB)
    if (DEFINED TBB_DIR)
        find_package(TBB REQUIRED COMPONENTS tbb)
        target_link_libraries(project_options INTERFACE TBB::tbb)
    else ()
        FetchContent_Declare(
            tbb URL https://github.com/wjakob/tbb/archive/806df70ee69fc7b332fcf90a48651f6dbf0663ba.tar.gz
            URL_HASH MD5=63fda89e88d34da63ddcef472e7725ef
            )
        if (NOT tbb_POPULATED)
            FetchContent_Populate(tbb)
            add_subdirectory(${tbb_SOURCE_DIR} ${tbb_BINARY_DIR})
            set_target_properties(tbb_static PROPERTIES COMPILE_FLAGS "-w")
            include_directories(SYSTEM ${tbb_SOURCE_DIR}/include)
            target_link_libraries(project_options INTERFACE tbb_static)
        endif ()
    endif ()
endif ()

############
# NANOBENCH
############

FetchContent_Declare(
    nanobench
    URL "https://github.com/martinus/nanobench/archive/v3.1.0.tar.gz"
    URL_HASH MD5=e646fb61164a60921c1a1834fbca24bc)
FetchContent_GetProperties(nanobench)
if(NOT nanobench_POPULATED)
    FetchContent_Populate(nanobench)
endif()

########
# EIGEN
########

#CPMAddPackage("gl:nlohmann/json@3.9.1")


CPMAddPackage(
    NAME Eigen
    VERSION 3.3.9
    URL https://gitlab.com/libeigen/eigen/-/archive/3.3.9/eigen-3.3.9.tar.gz
    # Eigen's CMakelists are not intended for library use
    DOWNLOAD_ONLY YES 
)
if(Eigen_ADDED)
    add_library(Eigen INTERFACE IMPORTED)
    target_include_directories(Eigen INTERFACE ${Eigen_SOURCE_DIR})
endif()

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
    URL "https://github.com/rollbear/trompeloeil/archive/v39.tar.gz"
    URL_HASH SHA256=10506e48abd605740bc9ed43e34059f5068bc80af14476bd129a3ed3b54d522f)
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
    URL_HASH MD5=922f0c5988c0f70c887d65b7cf2762ac)
FetchContent_GetProperties(coulombgalore)
if(NOT coulombgalore_POPULATED)
    FetchContent_Populate(coulombgalore)
endif()
include_directories(SYSTEM ${coulombgalore_SOURCE_DIR})

# Add third-party headers to include path. Note this is done with SYSTEM
# to disable potential compiler warnings

include_directories(SYSTEM ${eigen_SOURCE_DIR} ${modernjson_SOURCE_DIR}/include ${rangev3_SOURCE_DIR}/include
    ${doctest_SOURCE_DIR} ${trompeloeil_SOURCE_DIR}/include ${nanobench_SOURCE_DIR}/src/include
    ${Pybind11IncludeDir} ${DocoptIncludeDir} ${CppsidIncludeDir} ${XdrfileIncludeDir} ${SpdlogIncludeDir}
    ${ProgressTrackerIncludeDir} ${exprtk_SOURCE_DIR} ${cereal_SOURCE_DIR}/include ${zstr_SOURCE_DIR}/src)
