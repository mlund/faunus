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

fau_dipoleexample(example_stockmayer "./" stockmayer.cpp)
set_target_properties(example_stockmayer PROPERTIES OUTPUT_NAME "stockmayer")

fau_example(example_titrate "./" titrate.cpp)
set_target_properties(example_titrate PROPERTIES OUTPUT_NAME "titrate")

fau_example(example_polymers "./" polymers.cpp)
set_target_properties(example_polymers PROPERTIES OUTPUT_NAME "polymers")
add_test( example_polymers ${CMAKE_CURRENT_SOURCE_DIR}/polymers.run )

fau_example(example_slitpolymer "./" slitpolymer.cpp)
set_target_properties(example_slitpolymer PROPERTIES OUTPUT_NAME "slitpolymer")
add_test( example_slitpolymer ${CMAKE_CURRENT_SOURCE_DIR}/slitpolymer.run )

fau_example(example_membrane "./" membrane.cpp)
set_target_properties(example_membrane PROPERTIES OUTPUT_NAME "membrane")
add_test( example_membrane ${CMAKE_CURRENT_SOURCE_DIR}/membrane.run )

fau_example(example_bulk "./" bulk.cpp)
set_target_properties(example_bulk PROPERTIES OUTPUT_NAME "bulk")
add_test( example_bulk ${CMAKE_CURRENT_SOURCE_DIR}/bulk.run )

fau_example(example_HS "./" HS.cpp)
set_target_properties(example_HS PROPERTIES OUTPUT_NAME "HS")

