find_package(SWIG 2.0)

if (SWIG_FOUND)
  INCLUDE(${SWIG_USE_FILE})

  INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/include)
  INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

  SET_SOURCE_FILES_PROPERTIES("interface.i" PROPERTIES CPLUSPLUS ON)
  SET_SOURCE_FILES_PROPERTIES("interface.i" PROPERTIES SWIG_FLAGS "-w315,325,401,503 -includeall")

  SET(CMAKE_SWIG_FLAGS ${CMAKE_SWIG_FLAGS} "-DNO_EXPLICIT_OVERRIDE")

  if(ENABLE_BABEL AND OPENBABEL_FOUND)
    include_directories(${OPENBABEL_INCLUDE_DIR})
    add_definitions(-DBABEL)
    SET(CMAKE_SWIG_FLAGS ${CMAKE_SWIG_FLAGS} "-DBABEL")
  endif()

  find_package(PythonLibs)
  find_package(PythonInterp)

  if(PYTHONLIBS_FOUND)

    INCLUDE_DIRECTORIES(${PYTHON_INCLUDE_PATH})
    SWIG_ADD_MODULE(pyfaunus python "interface.i")
    SWIG_LINK_LIBRARIES(pyfaunus libfaunus ${PYTHON_LIBRARIES})

    if (PYTHONINTERP_FOUND)
      add_test(SWIG_python_module ${PYTHON_EXECUTABLE} "-c" "import pyfaunus")
      set_tests_properties(SWIG_python_module PROPERTIES WORKING_DIRECTORY "${CMAKE_BINARY_DIR}/swig")
    endif()

  endif()

  #find_package(Ruby)
  if(RUBY_FOUND)
    INCLUDE_DIRECTORIES(${RUBY_INCLUDE_DIRS})
    SWIG_ADD_MODULE(rubyfaunus ruby "interface.i")
    SWIG_LINK_LIBRARIES(rubyfaunus libfaunus ${RUBY_LIBRARY})
    set_target_properties(rubyfaunus PROPERTIES EXCLUDE_FROM_ALL TRUE)
  endif()

  #find_package(TCL)
  if(TCL_FOUND)
    INCLUDE_DIRECTORIES(${TCL_INCLUDE_PATH})
    SWIG_ADD_MODULE(tclfaunus tcl "interface.i")
    SWIG_LINK_LIBRARIES(tclfaunus libfaunus ${TCL_LIBRARY})
    set_target_properties(tclfaunus PROPERTIES EXCLUDE_FROM_ALL TRUE)
  endif()
endif()

