option(ENABLE_FREESASA "Fetch 3rd-party SASA calculation software" on)

include(ExternalProject)
include(FetchContent)

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


#########
# SPDLOG
#########

ExternalProject_Add(
    project_spdlog
    PREFIX "${CMAKE_CURRENT_BINARY_DIR}/_deps"
    LOG_DOWNLOAD ON
    URL https://github.com/gabime/spdlog/archive/v1.6.1.tar.gz
    URL_HASH SHA256=378a040d91f787aec96d269b0c39189f58a6b852e4cbf9150ccfacbe85ebbbfc
    CMAKE_ARGS -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER} -DSPDLOG_INSTALL=off -DCMAKE_POSITION_INDEPENDENT_CODE=on
    BUILD_COMMAND ${CMAKE_MAKE_PROGRAM} spdlog
    INSTALL_COMMAND ""
)
ExternalProject_Get_Property(project_spdlog source_dir)
ExternalProject_Get_Property(project_spdlog binary_dir)
set(SpdlogIncludeDir ${source_dir}/include)
add_library(spdlog STATIC IMPORTED)
set_property(TARGET spdlog PROPERTY IMPORTED_LOCATION ${binary_dir}/libspdlog.a)
add_dependencies(spdlog project_spdlog)

##############
# MODERN JSON
##############

FetchContent_Declare(
    modernjson
    URL "https://github.com/nlohmann/json/releases/download/v3.8.0/include.zip"
    URL_HASH SHA256=8590fbcc2346a3eefc341935765dd57598022ada1081b425678f0da9a939a3c0)
FetchContent_GetProperties(modernjson)
if(NOT modernjson_POPULATED)
    FetchContent_Populate(modernjson)
endif()
add_compile_definitions("NLOHMANN_JSON_HPP") # older versions used this macro. Now it's suffixed with "_"

#########
# CEREAL
#########

FetchContent_Declare(
    cereal
    URL "https://github.com/USCiLab/cereal/archive/v1.3.0.tar.gz"
    URL_HASH MD5=4342e811f245403646c4175258f413f1)
FetchContent_GetProperties(cereal)
if(NOT cereal_POPULATED)
    FetchContent_Populate(cereal)
endif()

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

###########
# RANGE-V3
###########

FetchContent_Declare(
    rangev3
    URL "https://github.com/ericniebler/range-v3/archive/0.10.0.tar.gz"
    URL_HASH MD5=2a68ca385f70d62398365cc2923d0573)

FetchContent_GetProperties(rangev3)
if(NOT rangev3_POPULATED)
    FetchContent_Populate(rangev3)
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

option(ENABLE_SID "Enable SID emulation" off)
if(ENABLE_SID)
    find_package(SDL2 CONFIG)
endif()


###########
# PYBIND11
###########

FetchContent_Declare(
    pybind11
    URL https://github.com/pybind/pybind11/archive/v2.5.0.tar.gz
    URL_HASH MD5=1ad2c611378fb440e8550a7eb6b31b89)

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

FetchContent_Declare(
    eigen
    URL "http://bitbucket.org/eigen/eigen/get/3.3.7.tar.gz"
    URL_HASH MD5=f2a417d083fe8ca4b8ed2bc613d20f07)
FetchContent_GetProperties(eigen)
if(NOT eigen_POPULATED)
    FetchContent_Populate(eigen)
endif()
 
##########
# XRDFILE
##########

ExternalProject_Add(
    project_xdrfile
    PREFIX "${CMAKE_CURRENT_BINARY_DIR}/_deps"
    URL "https://github.com/wesbarnett/libxdrfile/archive/2.1.2.tar.gz"
    URL_MD5 ee114404b4a01613b2f0167a2ad92536
    PATCH_COMMAND echo "add_library(xdrfile-static STATIC \${SRCS})" >> CMakeLists.txt
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

##########
# DOCTEST
##########

FetchContent_Declare(
    doctest
    URL "https://github.com/onqtam/doctest/archive/2.4.0.tar.gz"
    URL_MD5 0b679d6294f97be4fb11cb26e801fda6)
FetchContent_GetProperties(doctest)
if(NOT doctest_POPULATED)
    FetchContent_Populate(doctest)
    add_definitions(-DDOCTEST_CONFIG_DISABLE)
endif()

###########
# FREESASA
###########

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

include_directories(SYSTEM ${eigen_SOURCE_DIR} ${doctest_SOURCE_DIR} ${modernjson_SOURCE_DIR}/include ${rangev3_SOURCE_DIR}/include
    ${nanobench_SOURCE_DIR}/src/include
    ${Pybind11IncludeDir} ${DocoptIncludeDir} ${CppsidIncludeDir} ${XdrfileIncludeDir} ${SpdlogIncludeDir}
    ${ProgressTrackerIncludeDir} ${exprtk_SOURCE_DIR} ${cereal_SOURCE_DIR}/include ${zstr_SOURCE_DIR}/src)
