include_directories(${CMAKE_SOURCE_DIR}/include)
if(ENABLE_MPI AND MPI_CXX_FOUND)
  fau_example(example_temper "./" temper.cpp)
  set_target_properties(example_temper PROPERTIES OUTPUT_NAME "temper")
  add_test( example_temper ${CMAKE_CURRENT_SOURCE_DIR}/temper.run ${MPIEXEC})
  set_tests_properties(example_temper PROPERTIES WORKING_DIRECTORY "${CMAKE_BINARY_DIR}/src/examples/")

  fau_example(example_manybody "./" manybody.cpp)
  set_target_properties(example_manybody PROPERTIES OUTPUT_NAME "manybody")
  #add_test( example_manybody ${CMAKE_CURRENT_SOURCE_DIR}/membrane.run )
endif()

#----- add examples to build process -----
fau_dipoleexample(unittests "./" unittests.cpp)
add_test( unittests ${CMAKE_BINARY_DIR}/src/examples/unittests )
set_target_properties(unittests PROPERTIES EXCLUDE_FROM_ALL false)

fau_example(example_minimal "./" minimal.cpp)
set_target_properties(example_minimal PROPERTIES OUTPUT_NAME "minimal")

fau_example(example_grand "./" grand.cpp)
set_target_properties(example_grand PROPERTIES OUTPUT_NAME "grand")
add_test( example_grand python ${CMAKE_CURRENT_SOURCE_DIR}/grand.py )

fau_dipoleexample(example_stockmayer "./" stockmayer.cpp)
fau_dipoleexample(example_stockmayer_polarizable "./" stockmayer.cpp)
set_target_properties(example_stockmayer PROPERTIES OUTPUT_NAME "stockmayer" EXCLUDE_FROM_ALL FALSE)
set_target_properties(example_stockmayer_polarizable PROPERTIES OUTPUT_NAME "stockmayer_polarizable")
set_target_properties(example_stockmayer_polarizable PROPERTIES COMPILE_DEFINITIONS "POLARIZE;DIPOLEPARTICLE")
add_test( example_stockmayer ${CMAKE_CURRENT_SOURCE_DIR}/stockmayer.run )

fau_example(example_titrate "./" titrate.cpp)
set_target_properties(example_titrate PROPERTIES OUTPUT_NAME "titrate")

fau_example(example_stripes "./" stripes.cpp)
set_target_properties(example_stripes PROPERTIES OUTPUT_NAME "stripes")

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

fau_example(example_bulk_coulomb "./" bulk.cpp)
set_target_properties(example_bulk_coulomb PROPERTIES OUTPUT_NAME "bulk_coulomb" EXCLUDE_FROM_ALL TRUE)
set_target_properties(example_bulk_coulomb PROPERTIES COMPILE_DEFINITIONS "COULOMB")

fau_example(example_bulk_dh "./" bulk.cpp)
set_target_properties(example_bulk_dh PROPERTIES OUTPUT_NAME "bulk_dh" EXCLUDE_FROM_ALL TRUE)
set_target_properties(example_bulk_dh PROPERTIES COMPILE_DEFINITIONS "DEBYEHUCKEL")

fau_example(example_water "./" water.cpp)
set_target_properties(example_water PROPERTIES OUTPUT_NAME "water")
add_test( example_water ${CMAKE_CURRENT_SOURCE_DIR}/water.run )

fau_example(example_water2 "./" water2.cpp)
set_target_properties(example_water2 PROPERTIES OUTPUT_NAME "water2")
add_test( example_water2 ${CMAKE_CURRENT_SOURCE_DIR}/water2.run )


fau_cigarexample(example_cigarsenergy "./" cigarsenergy.cpp)
set_target_properties(example_cigarsenergy PROPERTIES OUTPUT_NAME "cigarsenergy" EXCLUDE_FROM_ALL FALSE)
add_test( example_cigarsenergy ${CMAKE_CURRENT_SOURCE_DIR}/cigarsenergy.run )

fau_cigarexample(example_cigars2fibrils "./" cigars2fibrils.cpp)
set_target_properties(example_cigars2fibrils PROPERTIES OUTPUT_NAME "cigars2fibrils" EXCLUDE_FROM_ALL FALSE)

