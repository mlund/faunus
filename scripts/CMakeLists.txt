if (ENABLE_PYTHON)
    INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/include)
    INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

    find_package(PythonLibs)
    find_package(PythonInterp)

    if (PYTHONLIBS_FOUND)
        add_library(pyfaunus SHARED pyfaunus.cpp)
        set_target_properties(pyfaunus PROPERTIES PREFIX "")
        set_target_properties(pyfaunus PROPERTIES SUFFIX ".so")
        target_link_libraries(pyfaunus libfaunus ${PYTHON_LIBRARY})
        INCLUDE_DIRECTORIES(${PYTHON_INCLUDE_PATH})
    endif ()

endif ()
