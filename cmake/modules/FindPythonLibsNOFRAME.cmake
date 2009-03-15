
# FindPythonLibsNOFRAME.cmake - Copyright (C) 2008 Mikael Lund 
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
#
#  PYTHONLIBS_FOUND    = found?
#  PYTHON_LIBRARIES    = path to the python library
#  PYTHON_INCLUDE_PATH = path to where Python.h is found

set(CMAKE_FIND_FRAMEWORK "NEVER")

find_path(
  PYTHON_INCLUDE_PATH Python.h
  HINTS /opt/local/ /usr/local /sw /usr
  PATH_SUFFIXES include/python2.6 include/python2.5 include/python2.4
  )

find_library(
  PYTHON_LIBRARIES
  NAMES python2.6 python2.5 python2.4
  HINTS /opt/local /usr/local /sw /usr
  PATH_SUFFIXES lib lib64
  )

include(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(PythonLibs DEFAULT_MSG PYTHON_LIBRARIES PYTHON_INCLUDE_PATH)

