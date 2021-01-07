include(ExternalProject)
include(FetchContent)

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


############
# INTEL TBB
############

option(ENABLE_TBB "Enable Intel TBB" off)
if (ENABLE_TBB)
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

include_directories(SYSTEM ${nanobench_SOURCE_DIR}/src/include
    ${CppsidIncludeDir} ${XdrfileIncludeDir} ${ProgressTrackerIncludeDir} ${cereal_SOURCE_DIR}/include ${zstr_SOURCE_DIR}/src)
