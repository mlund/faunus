include_directories(${CMAKE_SOURCE_DIR}/include)

#----- add examples to build process -----
fau_example(example_minimal "./" minimal.cpp)
set_target_properties(example_minimal PROPERTIES OUTPUT_NAME "minimal")

fau_example(example_polymers "./" polymers.cpp)
set_target_properties(example_polymers PROPERTIES OUTPUT_NAME "polymers")
add_test( example_polymers ${CMAKE_CURRENT_SOURCE_DIR}/polymers.run )


