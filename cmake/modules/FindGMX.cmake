# - Find gromacs gmx lib
# Find the native gromacs includes and library
#
#  GMX_INCLUDE_DIR - where to find include files
#  GMX_LIBRARIES   - List of libraries
#  GMX_FOUND       - True if gromacs is found.
#  GMX_MOTIF       - True if Motif/X is required for linking.

SET(GMX_MOTIF FALSE)
IF (GMX_INCLUDE_DIR)
  # Already in cache, be silent
  SET(GMX_FIND_QUIETLY TRUE)
ENDIF (GMX_INCLUDE_DIR)

FIND_PATH( GMX_INCLUDE_DIR xtcio.h
           HINTS "$ENV{GMXLDLIB}/../" /opt/local /usr/local /sw /usr
           PATH_SUFFIXES include/gromacs )

SET(GMX_NAMES gmx)
FIND_LIBRARY( GMX_LIBRARIES NAMES ${GMX_NAMES}
              HINTS "$ENV{GMXLDLIB}/../" /opt/local /usr/local /sw /usr
              PATH_SUFFIXES lib lib64 )

if(GMX_INCLUDE_DIR AND GMX_LIBRARIES)
  execute_process(
    COMMAND ${CMAKE_AR} -t ${GMX_LIBRARIES} mgmx.o
    RESULT_VARIABLE ARRC OUTPUT_QUIET ERROR_QUIET)
  if (ARRC EQUAL 0)
    set(GMX_MOTIF TRUE)
  endif (ARRC EQUAL 0)
else(GMX_INCLUDE_DIR AND GMX_LIBRARIES)
  set(GMX_LIBRARIES)
  set(GMX_INCLUDE_DIR)
endif(GMX_INCLUDE_DIR AND GMX_LIBRARIES)

if (GMX_MOTIF) 
  find_package(X11)
  find_package(Motif)
  if (X11_FOUND AND MOTIF_FOUND)
    set(GMX_LIBRARIES ${GMX_LIBRARIES} ${X11_Xt_LIB} ${MOTIF_LIBRARIES})
    set(GMX_INCLUDE_DIR ${GMX_INCLUDE_DIR} ${X11_INCLUDE_DIR} ${MOTIF_INCLUDE_DIR})
  else (X11_FOUND AND MOTIF_FOUND)
    set (GMX_FOUND FALSE)
    set (GMX_LIBRARIES)
    set (GMX_INCLUDE_DIR)
  endif (X11_FOUND AND MOTIF_FOUND)
endif (GMX_MOTIF)

include(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(GMX DEFAULT_MSG GMX_LIBRARIES GMX_INCLUDE_DIR)


