include_directories(${CMAKE_SOURCE_DIR}/include)

# -------------------------------
#   Faunus library object files
# -------------------------------
set(objs slump physconst inputfile point io titrate container
    hardsphere group particles analysis species widom histogram
    mcloop ensemble faunus messagepool
    moves/base moves/translational moves/rotational moves/charge
    moves/miscmove moves/saltbath moves/clustermove #moves/rosenbluth 
    moves/crankshaft moves/membrane moves/reptation moves/switch
    moves/replicaexchange moves/eqtitrate notification) 
set_source_files_properties(${objs} PROPERTIES LANGUAGE CXX)

#----- xdrfile static library for gromacs output -----
add_subdirectory(xdrfile)

# ---=---------------------------------
#   Mersenne Twister requires C++ TR1
# ----=--------------------------------
check_include_file_cxx("tr1/random" HAVETR1)
if(HAVETR1)
  add_definitions(-DHAVETR1)
  if(ENABLE_TWISTER)
    add_definitions(-DMERSENNETWISTER)
  endif()
endif()

# -------------------------------------
#   Get subversion revision of source
# -------------------------------------
if(Subversion_FOUND)
  Subversion_WC_INFO(${CMAKE_HOME_DIRECTORY} Faunus ERROR_QUIET)
  if(Faunus_WC_FOUND)
    add_definitions("-D__SVN_REV__=${Faunus_WC_REVISION}")
  endif()
endif()

# -------------------------------------
#   Growl notification (MacOS X)
# -------------------------------------
find_program( GROWLNOTIFYEXE NAMES "growlnotify" HINTS "/usr/local/bin/" )
if(GROWLNOTIFYEXE)
  add_definitions(-DGROWLNOTIFYEXE="${GROWLNOTIFYEXE}")
endif()

# -----------------------
#   Link with openbabel
# -----------------------
if(ENABLE_BABEL AND OPENBABEL_FOUND)
  set(objs ${objs} iobabel)
  include_directories(${OPENBABEL_INCLUDE_DIR})
  set(LINKLIBS ${LINKLIBS} ${OPENBABEL_LIBRARIES})
  add_definitions(-DBABEL)
endif()

# -----------------------------
#   Link with fortran library
# -----------------------------
if(ENABLE_FORTRAN AND CMAKE_Fortran_COMPILER_WORKS)
  set(LINKLIBS ${LINKLIBS} fortranfunc)
endif()

# --------------------
#   Faunus libraries
# --------------------
# static
add_library(libfaunus.a STATIC EXCLUDE_FROM_ALL ${objs})
set_target_properties(libfaunus.a PROPERTIES OUTPUT_NAME faunus)
target_link_libraries(libfaunus.a xdrfile ${LINKLIBS})

# dynamic
add_library(libfaunus SHARED ${objs})
set_target_properties(libfaunus PROPERTIES
                      OUTPUT_NAME faunus
                      )
target_link_libraries(libfaunus xdrfile ${LINKLIBS})
install(TARGETS libfaunus LIBRARY DESTINATION lib)

# dynamic for hypersphere
add_library(libhyperfaunus SHARED ${objs})
set_target_properties(libhyperfaunus PROPERTIES
                      OUTPUT_NAME hyperfaunus
                      COMPILE_DEFINITIONS "HYPERSPHERE"
                      EXCLUDE_FROM_ALL TRUE
                      )
target_link_libraries(libhyperfaunus xdrfile ${LINKLIBS})

#----- header install target -----
install(DIRECTORY ${CMAKE_SOURCE_DIR}/include/faunus
        DESTINATION include
        FILES_MATCHING PATTERN "*.h"
        PATTERN ".svn" EXCLUDE
        PATTERN "CMakeFiles" EXCLUDE)

#----- Add example programs -----
add_subdirectory(examples)

