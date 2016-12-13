include_directories(${CMAKE_SOURCE_DIR}/include)
if (ENABLE_MPI AND MPI_CXX_FOUND)
    fau_example(example_temper "./" temper.cpp)
    set_target_properties(example_temper PROPERTIES OUTPUT_NAME "temper")
    add_test(example_temper ${CMAKE_CURRENT_SOURCE_DIR}/temper.run ${MPIEXEC})
    set_tests_properties(example_temper PROPERTIES WORKING_DIRECTORY "${CMAKE_BINARY_DIR}/src/examples/")

    fau_example(example_manybody "./" manybody.cpp)
    set_target_properties(example_manybody PROPERTIES OUTPUT_NAME "manybody")
    #add_test( example_manybody ${CMAKE_CURRENT_SOURCE_DIR}/membrane.run )
endif ()

#----- add examples to build process -----
fau_example(unittests "./" unittests.cpp)
add_test(unittests ${CMAKE_BINARY_DIR}/src/examples/unittests
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
set_target_properties(unittests PROPERTIES EXCLUDE_FROM_ALL false)

fau_example(example_minimal "./" minimal.cpp)
set_target_properties(example_minimal PROPERTIES OUTPUT_NAME "minimal")

fau_example(example_penalty "./" penalty.cpp)
set_target_properties(example_penalty PROPERTIES OUTPUT_NAME "penalty")
add_test(example_penalty python ${CMAKE_CURRENT_SOURCE_DIR}/penalty.py)

fau_example(example_grand "./" grand.cpp)
set_target_properties(example_grand PROPERTIES OUTPUT_NAME "grand")
add_test(example_grand python ${CMAKE_CURRENT_SOURCE_DIR}/grand.py)

fau_example(example_gctit "./" gctit.cpp)
set_target_properties(example_gctit PROPERTIES OUTPUT_NAME "gctit")
#add_test( example_gctit python ${CMAKE_CURRENT_SOURCE_DIR}/gctit.py )

fau_example(example_stockmayer "./" stockmayer.cpp)
fau_example(example_stockmayer_polarizable "./" stockmayer.cpp)
set_target_properties(example_stockmayer PROPERTIES OUTPUT_NAME "stockmayer" EXCLUDE_FROM_ALL FALSE)
set_target_properties(example_stockmayer_polarizable PROPERTIES OUTPUT_NAME "stockmayer_polarizable")
set_target_properties(example_stockmayer_polarizable PROPERTIES COMPILE_DEFINITIONS "POLARIZE;DIPOLEPARTICLE")
add_test(example_stockmayer ${CMAKE_CURRENT_SOURCE_DIR}/stockmayer.run)

#fau_example(example_titrate "./" titrate.cpp)
#set_target_properties(example_titrate PROPERTIES OUTPUT_NAME "titrate")

fau_example(example_stripes "./" stripes.cpp)
set_target_properties(example_stripes PROPERTIES OUTPUT_NAME "stripes")

fau_example(example_polymers "./" polymers.cpp)
set_target_properties(example_polymers PROPERTIES OUTPUT_NAME "polymers")
add_test(example_polymers python ${CMAKE_CURRENT_SOURCE_DIR}/polymers.py)

# gcmol
fau_example(example_gcmol "./" gcmol.cpp)
set_target_properties(example_gcmol PROPERTIES OUTPUT_NAME "gcmol")
add_test(example_gcmol python ${CMAKE_CURRENT_SOURCE_DIR}/gcmol.py)
add_custom_command(TARGET example_gcmol PRE_LINK COMMAND python gcmol.py --norun
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})

fau_example(example_slitpolymer "./" slitpolymer.cpp)
set_target_properties(example_slitpolymer PROPERTIES OUTPUT_NAME "slitpolymer")
add_test(example_slitpolymer python ${CMAKE_CURRENT_SOURCE_DIR}/slitpolymer.py)

fau_example(example_membrane "./" membrane.cpp)
set_target_properties(example_membrane PROPERTIES OUTPUT_NAME "membrane")
add_test(example_membrane ${CMAKE_CURRENT_SOURCE_DIR}/membrane.run)

fau_example(example_bulk "./" bulk.cpp)
set_target_properties(example_bulk PROPERTIES OUTPUT_NAME "bulk")
add_test(example_bulk python ${CMAKE_CURRENT_SOURCE_DIR}/bulk.py)

fau_example(example_seawater "./" seawater.cpp)
set_target_properties(example_seawater PROPERTIES OUTPUT_NAME "seawater")
add_test(example_seawater ${CMAKE_CURRENT_SOURCE_DIR}/seawater.sh)

fau_example(example_twobody "./" twobody.cpp)
set_target_properties(example_twobody PROPERTIES OUTPUT_NAME "twobody")

fau_example(example_keesom "./" keesom.cpp)
set_target_properties(example_keesom PROPERTIES OUTPUT_NAME "keesom")

fau_example(example_bulk_coulomb "./" bulk.cpp)
set_target_properties(example_bulk_coulomb PROPERTIES OUTPUT_NAME "bulk_coulomb" EXCLUDE_FROM_ALL TRUE)
set_target_properties(example_bulk_coulomb PROPERTIES COMPILE_DEFINITIONS "COULOMB")

fau_example(example_bulk_dh "./" bulk.cpp)
set_target_properties(example_bulk_dh PROPERTIES OUTPUT_NAME "bulk_dh" EXCLUDE_FROM_ALL TRUE)
set_target_properties(example_bulk_dh PROPERTIES COMPILE_DEFINITIONS "DEBYEHUCKEL")

# water2
fau_example(example_water2 "./" water2.cpp)
set_target_properties(example_water2 PROPERTIES OUTPUT_NAME "water2")
add_test(example_water2 python ${CMAKE_CURRENT_SOURCE_DIR}/water2.py)
add_custom_command(TARGET example_water2 PRE_LINK COMMAND python water2.py --norun
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})

fau_example(example_water2_ewald "./" water2.cpp)
set_target_properties(example_water2_ewald PROPERTIES
        OUTPUT_NAME "water2" COMPILE_DEFINITIONS "EWALD" EXCLUDE_FROM_ALL TRUE)

fau_example(example_cigarsenergy "./" cigarsenergy.cpp)
set_target_properties(example_cigarsenergy PROPERTIES OUTPUT_NAME "cigarsenergy" EXCLUDE_FROM_ALL FALSE)
add_test(example_cigarsenergy ${CMAKE_CURRENT_SOURCE_DIR}/cigarsenergy.run)

fau_example(example_cigars2fibrils "./" cigars2fibrils.cpp)
set_target_properties(example_cigars2fibrils PROPERTIES OUTPUT_NAME "cigars2fibrils" EXCLUDE_FROM_ALL FALSE)
add_test(example_cigars2fibrils ${CMAKE_CURRENT_SOURCE_DIR}/cigars2fibrils.run)
