include_directories(${CMAKE_SOURCE_DIR}/include)

add_executable(pythontest pythontest.cpp)
set_target_properties(pythontest PROPERTIES
                      RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}")
target_link_libraries(pythontest libfaunus)

