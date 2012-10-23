include_directories(${CMAKE_SOURCE_DIR}/include)
if(ENABLE_MPI AND MPI_CXX_FOUND)
  fau_example(example_temper "./" temper.cpp)
  set_target_properties(example_temper PROPERTIES OUTPUT_NAME "temper")
  add_test( example_temper ${CMAKE_CURRENT_SOURCE_DIR}/temper.run ${MPIEXEC})
  set_tests_properties(example_temper PROPERTIES WORKING_DIRECTORY "${CMAKE_BINARY_DIR}/src/examples/")
endif()

#----- add examples to build process -----
fau_example(unittests "./" unittests.cpp)
add_test( unittests ${CMAKE_BINARY_DIR}/src/examples/unittests )

fau_example(example_minimal "./" minimal.cpp)
set_target_properties(example_minimal PROPERTIES OUTPUT_NAME "minimal")

fau_example(example_polymers "./" polymers.cpp)
set_target_properties(example_polymers PROPERTIES OUTPUT_NAME "polymers")
add_test( example_polymers ${CMAKE_CURRENT_SOURCE_DIR}/polymers.run )

fau_example(example_slitpolymer "./" slitpolymer.cpp)
set_target_properties(example_slitpolymer PROPERTIES OUTPUT_NAME "slitpolymer")
add_test( example_slitpolymer ${CMAKE_CURRENT_SOURCE_DIR}/slitpolymer.run )


