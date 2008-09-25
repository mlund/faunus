# FindGMX.cmake - Copyright (C) 2008 Mikael Lund 
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or 
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along
# with this program; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.

# Try to find Gromacs includes and libgmx library
#   GMX_INCLUDE_DIR - where to find include files
#   GMX_LIBRARIES   - List of libraries
#   GMX_FOUND       - True if gromacs is found.
#   GMX_MOTIF       - True if Motif/X11 is needed.

if (GMX_INCLUDE_DIR)
  # Already in cache, be silent
  set(GMX_FIND_QUIETLY TRUE)
else (GMX_INCLUDE_DIR)
  set(GMX_MOTIF FALSE)

  find_path( GMX_INCLUDE_DIR xtcio.h
    HINTS "$ENV{GMXLDLIB}/../" /opt/local /usr/local /sw /usr
    PATH_SUFFIXES include/gromacs gromacs/include )

  set(GMX_NAMES gmx gmx_d gmx_mpi)
  find_library( GMX_LIBRARIES NAMES ${GMX_NAMES}
    HINTS "$ENV{GMXLDLIB}/../" /opt/local /usr/local /sw /usr
    PATH_SUFFIXES lib lib64 gromacs/lib )

  # Does libgmx need X11 and Motif?
  if(GMX_INCLUDE_DIR AND GMX_LIBRARIES)
    execute_process(
      COMMAND ${CMAKE_AR} -t ${GMX_LIBRARIES} mgmx.o
      RESULT_VARIABLE ARRC
      OUTPUT_QUIET ERROR_QUIET )
    if (ARRC EQUAL 0)
      set(GMX_MOTIF TRUE)
    endif (ARRC EQUAL 0)
  else(GMX_INCLUDE_DIR AND GMX_LIBRARIES)
    set(GMX_LIBRARIES)
    set(GMX_INCLUDE_DIR)
  endif(GMX_INCLUDE_DIR AND GMX_LIBRARIES)

  # Add X11 and Motif if needed
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
endif (GMX_INCLUDE_DIR)

