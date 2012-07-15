cmake_minimum_required(VERSION 2.8 FATAL_ERROR)
project("Faunus")
set(CMAKE_LEGACY_CYGWIN_WIN32 0) # Remove when CMake >= 2.8.4 is required

# ----------------------------
#   Locate external packages
# ----------------------------
set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/modules)
find_package(OpenBabel2)
#find_package(PythonLibs)
find_package(Subversion-patched)
find_package(MPI)
include(CTest)
include(CheckIncludeFileCXX)
include(CheckCXXSourceCompiles)
enable_testing()

# ------------------------------
#   Default build is "Release"
# ------------------------------
include(OptimizationFlags)
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Debug
    CACHE STRING "Build type: None Debug Release RelWithDebInfo MinSizeRel." FORCE)
endif()

# -----------------
#   User options
# -----------------
option(ENABLE_BABEL    "Try to use OpenBabel for file I/O" off)
option(ENABLE_OPENMP   "Try to use OpenMP parallization" off)
option(ENABLE_FORTRAN  "Try to compile legacy fortran code" off)
option(ENABLE_MPI      "Enable MPI parallel code (experimental!)" off)
option(ENABLE_TWISTER  "Enable Mersenne Twister random number generator" off)
#option(ENABLE_SWIG     "Try to create SWIG modules for python, tcl, ruby etc. (experimental!)" off)
option(ENABLE_FAST_EXP "Approximate e^x function for certain pair potentials." off)
option(ENABLE_UNICODE  "Use unicode characters in output" on)
mark_as_advanced( CLEAR CMAKE_VERBOSE_MAKEFILE CMAKE_CXX_COMPILER CMAKE_CXX_FLAGS )
mark_as_advanced( EXECUTABLE_OUTPUT_PATH LIBRARY_OUTPUT_PATH
                  CMAKE_OSX_ARCHITECTURES CMAKE_OSX_SYSROOT GCCXML DART_TESTING_TIMEOUT)

# -------------------------------
#   Check C/C++ code generation
# -------------------------------
CHECK_CXX_SOURCE_COMPILES("int main() { unsigned long long int i; }" unsigned_long_long_int)
CHECK_CXX_SOURCE_COMPILES("int main() { long double f; }" long_double)
CHECK_CXX_SOURCE_COMPILES("int main() { int a[0]; for (int &i:a) i++; }" range_based_for_loop)
CHECK_CXX_SOURCE_COMPILES("int main() { auto i=0; }" auto_type_inference)
CHECK_CXX_SOURCE_COMPILES("int main() { nullptr; }" nullptr)
CHECK_CXX_SOURCE_COMPILES("#include<memory>\nint main() {std::shared_ptr<int> i;}" smart_pointers)
CHECK_CXX_SOURCE_COMPILES("int main() { struct a { virtual void f() final; }; }" final_keyword)

if(NOT (auto_type_inference AND unsigned_long_long_int AND long_double AND range_based_for_loop AND nullptr AND smart_pointers) )
  message(FATAL_ERROR "Required C++11 features unsupported by compiler.")
endif()
if (NOT (final_keyword))
  add_definitions(-DNO_EXPLICIT_OVERRIDE)
endif()

# --------------------------------
#   Add Revision control targets
# --------------------------------
if(Subversion_FOUND)
  add_custom_target(rccommit svn commit)

endif()
 
# ---------------------------
#   Add source directories
# ---------------------------
add_subdirectory(tests)
#add_subdirectory(legacy/fortran)
add_subdirectory(src)
add_subdirectory(doc EXCLUDE_FROM_ALL)

#if(ENABLE_SWIG)
#  add_subdirectory(misc/swig)
#endif()

if (CMAKE_BUILD_TYPE MATCHES "Debug")
  message(
    "...............................................................\n"
    "  Faunus is compiled in Debug mode to enable checks and\n"
    "  assertions highly recommended when developing and testing.\n"
    "  For production, switch to Release for optimal performance:\n\n"
    "  -DCMAKE_BUILD_TYPE=Release\n"
    "***************************************************************\n")
endif()
if (CMAKE_BUILD_TYPE MATCHES "Release")
  message(
    "......................................................................\n"
    "  Faunus is running in Release mode recommended for production runs.\n"
    "  For setting up simulations it is HIGHLY recommended to use Debug\n"
    "  mode to enable run-time checks and assertions:\n\n"
    "  -DCMAKE_BUILD_TYPE=Debug\n"
    "**********************************************************************\n")
endif()

