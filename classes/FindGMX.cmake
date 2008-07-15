# - Find gromacs gmx lib
# Find the native gromacs includes and library
#
#  GMX_INCLUDE_DIR - where to find include files
#  GMX_LIBRARIES   - List of libraries when using zlib.
#  GMX_FOUND       - True if gromacs is found.

IF (GMX_INCLUDE_DIR)
  # Already in cache, be silent
  SET(GMX_FIND_QUIETLY TRUE)
ENDIF (GMX_INCLUDE_DIR)

FIND_PATH(GMX_INCLUDE_DIR xtcio.h
  PATHS
  "$ENV{GMXLDLIB}/../"
  /usr/local /opt/local /sw /usr
  PATH_SUFFIXES include/gromacs
  NO_DEFAULT_PATH
)

SET(GMX_NAMES gmx_mpi gmx)
FIND_LIBRARY(GMX_LIBRARY NAMES ${GMX_NAMES}
  PATHS
  "$ENV{GMXLDLIB}/../"
  /usr/local /opt/local /sw /usr
  PATH_SUFFIXES lib lib64
  NO_DEFAULT_PATH
)

# handle the QUIETLY and REQUIRED arguments and set GMX_FOUND to TRUE if 
# all listed variables are TRUE
INCLUDE(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(GMX DEFAULT_MSG GMX_LIBRARY GMX_INCLUDE_DIR)

IF(GMX_FOUND)
  SET( GMX_LIBRARIES ${GMX_LIBRARY} )
ELSE(GMX_FOUND)
  SET( GMX_LIBRARIES )
ENDIF(GMX_FOUND)

MARK_AS_ADVANCED( GMX_LIBRARY GMX_INCLUDE_DIR )
