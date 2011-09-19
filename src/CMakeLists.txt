include_directories(${CMAKE_SOURCE_DIR}/include)

# -------------------------------
#   Faunus library object files
# -------------------------------
set(objs
  titrate
  drift energy faunus geometry group histogram inputfile io mcloop move
  physconst point potentials slump space species)
set_source_files_properties(${objs} PROPERTIES LANGUAGE CXX)

#----- xdrfile static library for gromacs output -----
add_subdirectory(xdrfile)

# -------------------------------------
#   Check for C++11
# -------------------------------------
check_include_file_cxx("tuple" TUPLE)
if (NOT TUPLE)
  message(FATAL_ERROR "Compiler does support C++11")
endif()

# -------------------------------------
#   Mersenne Twister requires
# ----=--------------------------------
if(ENABLE_TWISTER)
  check_include_file_cxx("random" HAVERANDOM)
  if(HAVERANDOM)
    add_definitions(-DMERSENNETWISTER)
  endif()
endif()

# -------------------------------------
#   Use approximate exp() funtion?
# -------------------------------------
if(ENABLE_FAST_EXP)
  add_definitions(-DFAU_FAST_EXP)
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

# dynamic for cigarparticles
add_library(libcigarfaunus SHARED ${objs})
set_target_properties(libcigarfaunus PROPERTIES
                      OUTPUT_NAME cigarfaunus
                      COMPILE_DEFINITIONS "CIGARPARTICLE"
                      EXCLUDE_FROM_ALL TRUE
                      )
target_link_libraries(libcigarfaunus xdrfile ${LINKLIBS})


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

# -----------------------------------------
#   Add "cppcheck" static analyzer target
# -----------------------------------------
find_program(CPPCHECK_EXE "cppcheck")
if (CPPCHECK_EXE)
  add_custom_target(
    cppcheck
    COMMAND ${CPPCHECK_EXE}
    -I ${CMAKE_SOURCE_DIR}/include -j2 --quiet --enable=style
    ${CMAKE_SOURCE_DIR}/src/*.cpp ${CMAKE_SOURCE_DIR}/src/moves/*.cpp
    COMMENT "Running C++ static analyzer..."
    )
  set_target_properties(cppcheck PROPERTIES EXCLUDE_FROM_ALL TRUE )
endif()

#----- Add example programs -----
add_subdirectory(programs)

