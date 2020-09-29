#
# Most examples also act as tests
#
set(EXAMPLES_DIR ${CMAKE_SOURCE_DIR}/examples)
set(YASON ${CMAKE_SOURCE_DIR}/scripts/yason.py)
set(JSON_COMPARE ${CMAKE_SOURCE_DIR}/scripts/jsoncompare.py)
set(CSV_COMPARE ${CMAKE_SOURCE_DIR}/scripts/csvcompare.py)

add_test(
        NAME unittests
        COMMAND sh -c "$<TARGET_FILE:faunus> test --test-case-exclude=*Benchmark*")

if (TARGET pyfaunus)
    add_test(
            NAME pythontest
            COMMAND sh -c "${PYTHON_EXECUTABLE} pythontest.py -v"
            WORKING_DIRECTORY ${EXAMPLES_DIR})
    set_tests_properties(pythontest PROPERTIES ENVIRONMENT "PYTHONPATH=$<TARGET_FILE_DIR:pyfaunus>")
endif()

add_test(
        NAME swapconf
        COMMAND sh -c "${PYTHON_EXECUTABLE} ${YASON} swapconf.yml\
    | $<TARGET_FILE:faunus> --nobar --state swapconf.state.json\
    && ${PYTHON_EXECUTABLE} ${JSON_COMPARE} swapconf.out.json out.json --tol 0.05"
        WORKING_DIRECTORY ${EXAMPLES_DIR}/swapconf)

add_test(
        NAME doublelayer
        COMMAND sh -c "${PYTHON_EXECUTABLE} ${YASON} doublelayer.yml\
    | $<TARGET_FILE:faunus> --nobar --state doublelayer.state.json\
    && ${PYTHON_EXECUTABLE} ${JSON_COMPARE} doublelayer.out.json out.json --tol 0.05"
        WORKING_DIRECTORY ${EXAMPLES_DIR}/doublelayer)

add_test(
        NAME penalty
        COMMAND sh -c "jupyter nbconvert --execute penalty.ipynb --ExecutePreprocessor.timeout=None"
        WORKING_DIRECTORY ${EXAMPLES_DIR}/penalty)
set_property(TEST penalty APPEND PROPERTY ENVIRONMENT FAUNUS_EXECUTABLE=$<TARGET_FILE:faunus>)
set_property(TEST penalty APPEND PROPERTY ENVIRONMENT YASON_EXECUTABLE=${YASON})

add_test(
        NAME multipole
        COMMAND sh -c "jupyter nbconvert --execute multipole.ipynb --ExecutePreprocessor.timeout=None"
        WORKING_DIRECTORY ${EXAMPLES_DIR}/multipole)
set_property(TEST multipole APPEND PROPERTY ENVIRONMENT FAUNUS_EXECUTABLE=$<TARGET_FILE:faunus>)
set_property(TEST multipole APPEND PROPERTY ENVIRONMENT YASON_EXECUTABLE=${YASON})

add_test(
        NAME virial
        COMMAND sh -c "jupyter nbconvert --execute virial.ipynb --ExecutePreprocessor.timeout=None"
        WORKING_DIRECTORY ${EXAMPLES_DIR}/virial)
set_property(TEST virial APPEND PROPERTY ENVIRONMENT FAUNUS_EXECUTABLE=$<TARGET_FILE:faunus>)
set_property(TEST virial APPEND PROPERTY ENVIRONMENT YASON_EXECUTABLE=${YASON})

add_test(
        NAME minimal
        COMMAND sh -c "${PYTHON_EXECUTABLE} ${YASON} minimal.yml\
    | $<TARGET_FILE:faunus> --nobar --state minimal.state.json\
    && ${PYTHON_EXECUTABLE} ${JSON_COMPARE} minimal.out.json out.json --tol 0.02"
        WORKING_DIRECTORY ${EXAMPLES_DIR}/minimal)

add_test(
        NAME chain
        COMMAND sh -c "${PYTHON_EXECUTABLE} ${YASON} chain.yml\
    | $<TARGET_FILE:faunus> --nobar --state chain.state.json\
    && ${PYTHON_EXECUTABLE} ${JSON_COMPARE} chain.out.json out.json --tol 0.05"
        WORKING_DIRECTORY ${EXAMPLES_DIR}/chain)

add_test(
        NAME polymers
        COMMAND sh -c "${PYTHON_EXECUTABLE} ${YASON} polymers.yml\
    | $<TARGET_FILE:faunus> --nobar --state polymers.state.json\
    && ${PYTHON_EXECUTABLE} ${JSON_COMPARE} polymers.out.json out.json --tol 0.05"
        WORKING_DIRECTORY ${EXAMPLES_DIR}/polymers)

add_test(
        NAME speciation
        COMMAND sh -c "${PYTHON_EXECUTABLE} ${YASON} speciation.yml\
    | $<TARGET_FILE:faunus> --quiet --state speciation.state.json\
    && ${PYTHON_EXECUTABLE} ${JSON_COMPARE} speciation.out.json out.json --tol 0.02"
        WORKING_DIRECTORY ${EXAMPLES_DIR}/speciation)

add_test(
        NAME speciation_group
        COMMAND sh -c "${PYTHON_EXECUTABLE} ${YASON} speciation_group.yml\
    | $<TARGET_FILE:faunus> --quiet --state speciation_group.state.json\
    && ${PYTHON_EXECUTABLE} ${JSON_COMPARE} speciation_group.out.json out.json --tol 0.03"
        WORKING_DIRECTORY ${EXAMPLES_DIR}/speciation_group)

add_test(
        NAME speciation_swapconf
        COMMAND sh -c "${PYTHON_EXECUTABLE} ${YASON} speciation_swapconf.yml\
    | $<TARGET_FILE:faunus> --quiet --state speciation_swapconf.state.json\
    && ${PYTHON_EXECUTABLE} ${JSON_COMPARE} speciation_swapconf.out.json out.json --tol 0.015"
        WORKING_DIRECTORY ${EXAMPLES_DIR}/speciation_swapconf)

add_test(
        NAME titration
        COMMAND sh -c "${PYTHON_EXECUTABLE} ${YASON} titration.yml\
    | $<TARGET_FILE:faunus> --quiet --state titration.state.json\
    && ${PYTHON_EXECUTABLE} ${JSON_COMPARE} titration.out.json out.json --tol 0.05"
        WORKING_DIRECTORY ${EXAMPLES_DIR}/titration)

add_test(
        NAME calciumhydroxide
        COMMAND sh -c "${PYTHON_EXECUTABLE} ${YASON} calciumhydroxide.yml\
    | $<TARGET_FILE:faunus> --quiet --state calciumhydroxide.state.json\
    && ${PYTHON_EXECUTABLE} ${JSON_COMPARE} calciumhydroxide.out.json out.json --tol 0.02"
        WORKING_DIRECTORY ${EXAMPLES_DIR}/calciumhydroxide)

add_test(
        NAME calciumhydroxide_molecular
        COMMAND sh -c "${PYTHON_EXECUTABLE} ${YASON} calciumhydroxide_molecular.yml\
    | $<TARGET_FILE:faunus> --quiet --state calciumhydroxide_molecular.state.json\
    && ${PYTHON_EXECUTABLE} ${JSON_COMPARE} calciumhydroxide_molecular.out.json out.json --tol 0.02"
        WORKING_DIRECTORY ${EXAMPLES_DIR}/calciumhydroxide_molecular)

add_test(
        NAME gouychapman
        COMMAND sh -c "${PYTHON_EXECUTABLE} ${YASON} gouychapman.yml\
    | $<TARGET_FILE:faunus> --quiet --state gouychapman.state.json\
    && ${PYTHON_EXECUTABLE} ${JSON_COMPARE} gouychapman.out.json out.json --tol 0.015"
        WORKING_DIRECTORY ${EXAMPLES_DIR}/gouychapman)

add_test(
        NAME water-ewald-NOCHECKS
        COMMAND sh -c "${PYTHON_EXECUTABLE} ${YASON} ewald.yml | $<TARGET_FILE:faunus>"
        WORKING_DIRECTORY ${EXAMPLES_DIR}/water)

add_test(
        NAME phosphate
        COMMAND sh -c "${PYTHON_EXECUTABLE} ${YASON} phosphate.yml\
    | $<TARGET_FILE:faunus> --quiet --state phosphate.state.json\
    && ${PYTHON_EXECUTABLE} ${JSON_COMPARE} phosphate.out.json out.json --tol 0.02"
        WORKING_DIRECTORY ${EXAMPLES_DIR}/phosphate)

add_test(
        NAME sodiumphosphate
        COMMAND sh -c "${PYTHON_EXECUTABLE} ${YASON} sodiumphosphate.yml\
    | $<TARGET_FILE:faunus> --quiet --state sodiumphosphate.state.json\
    && ${PYTHON_EXECUTABLE} ${JSON_COMPARE} sodiumphosphate.out.json out.json --tol 0.05"
        WORKING_DIRECTORY ${EXAMPLES_DIR}/sodiumphosphate)

add_test(
        NAME cluster
        COMMAND sh -c "${PYTHON_EXECUTABLE} ${YASON} cluster.yml\
    | $<TARGET_FILE:faunus> --quiet --state cluster.state.json\
    && ${PYTHON_EXECUTABLE} ${JSON_COMPARE} cluster.out.json out.json --tol 0.1"
        WORKING_DIRECTORY ${EXAMPLES_DIR}/cluster)

add_test(
        NAME cluster-ideal
        COMMAND sh -c "${PYTHON_EXECUTABLE} ${YASON} cluster-ideal.yml\
    | $<TARGET_FILE:faunus> --quiet\
    && ${PYTHON_EXECUTABLE} ${JSON_COMPARE} cluster-ideal.out.json out.json --tol 0.01"
        WORKING_DIRECTORY ${EXAMPLES_DIR}/cluster)

add_test(
        NAME bulk
        COMMAND sh -c "${PYTHON_EXECUTABLE} ${YASON} bulk.yml\
    | $<TARGET_FILE:faunus> --quiet --state bulk.state.json\
    ; ${PYTHON_EXECUTABLE} ${JSON_COMPARE} bulk.out.json out.json --tol 0.05"
        WORKING_DIRECTORY ${EXAMPLES_DIR}/bulk)
set_property(TEST bulk PROPERTY FIXTURES_SETUP bulk)

add_test(
        NAME bulk-replay
        COMMAND sh -c "${PYTHON_EXECUTABLE} ${YASON} bulk-replay.yml\
    | $<TARGET_FILE:faunus> --quiet \
    ; diff -q ../bulk/traj.xtc traj.xtc \
    && ${CSV_COMPARE} -q --tol 0.04 --small 0.1 ../bulk/rdf.dat rdf.dat"
        WORKING_DIRECTORY ${EXAMPLES_DIR}/bulk-replay)
set_property(TEST bulk-replay PROPERTY FIXTURES_REQUIRED bulk)

add_test(
        NAME seawater
        COMMAND sh -c "${PYTHON_EXECUTABLE} ${YASON} seawater.yml\
    | $<TARGET_FILE:faunus> --quiet"
        WORKING_DIRECTORY ${EXAMPLES_DIR}/seawater)

add_test(
        NAME stockmayer
        COMMAND sh -c "${PYTHON_EXECUTABLE} ${YASON} stockmayer.yml\
    | $<TARGET_FILE:faunus> --quiet --state stockmayer.state.json\
    && ${PYTHON_EXECUTABLE} ${JSON_COMPARE} stockmayer.out.json out.json --tol 0.065"
        WORKING_DIRECTORY ${EXAMPLES_DIR}/stockmayer)

add_test(
        NAME membrane
        COMMAND sh -c "${PYTHON_EXECUTABLE} ${YASON} membrane.yml\
    | $<TARGET_FILE:faunus> --quiet --state membrane.state.json\
    && ${PYTHON_EXECUTABLE} ${JSON_COMPARE} membrane.out.json out.json --tol 0.065"
        WORKING_DIRECTORY ${EXAMPLES_DIR}/membrane)

install(FILES
        ${EXAMPLES_DIR}/README.md
        ${EXAMPLES_DIR}/bulk/bulk.yml
        ${EXAMPLES_DIR}/bulk-replay/bulk-replay.yml
        ${EXAMPLES_DIR}/cluster/cluster.yml
        ${EXAMPLES_DIR}/cluster/cluster.agr
        ${EXAMPLES_DIR}/cluster/cluster.out.json
        ${EXAMPLES_DIR}/cluster/cluster.state.json
        ${EXAMPLES_DIR}/cluster/water.xyz
        ${EXAMPLES_DIR}/doublelayer/doublelayer.yml
        ${EXAMPLES_DIR}/doublelayer/doublelayer.agr
        ${EXAMPLES_DIR}/doublelayer/doublelayer.out.json
        ${EXAMPLES_DIR}/doublelayer/doublelayer.state.json
        ${EXAMPLES_DIR}/doublelayer/doublelayer.mfcorr.dat
        ${EXAMPLES_DIR}/gouychapman/gouychapman.yml
        ${EXAMPLES_DIR}/gouychapman/gouychapman.out.json
        ${EXAMPLES_DIR}/gouychapman/gouychapman.state.json
        ${EXAMPLES_DIR}/isochoric.yml
        ${EXAMPLES_DIR}/isochoric.xyz
        ${EXAMPLES_DIR}/membrane/membrane.yml
        ${EXAMPLES_DIR}/membrane/membrane.out.json
        ${EXAMPLES_DIR}/membrane/membrane.state.json
        ${EXAMPLES_DIR}/lipid.xyz
        ${EXAMPLES_DIR}/minimal/minimal.yml
        ${EXAMPLES_DIR}/minimal/minimal.out.json
        ${EXAMPLES_DIR}/minimal/minimal.state.json
        ${EXAMPLES_DIR}/penalty/penalty.yml
        ${EXAMPLES_DIR}/penalty/penalty.ipynb
        ${EXAMPLES_DIR}/multipole/multipole.yml
        ${EXAMPLES_DIR}/multipole/multipole.ipynb
        ${EXAMPLES_DIR}/multipole/multipole.state.ubj
        ${EXAMPLES_DIR}/virial/virial.yml
        ${EXAMPLES_DIR}/virial/virial.ipynb
        ${EXAMPLES_DIR}/chain/chain.yml
        ${EXAMPLES_DIR}/chain/chain.xyz
        ${EXAMPLES_DIR}/chain/chain.state.json
        ${EXAMPLES_DIR}/chain/chain.out.json
        ${EXAMPLES_DIR}/polymers/polymers.yml
        ${EXAMPLES_DIR}/polymers/polymers.readme
        ${EXAMPLES_DIR}/polymers/polymers.pdf
        ${EXAMPLES_DIR}/polymers/polymers.state.json
        ${EXAMPLES_DIR}/polymers/polymers.out.json
        ${EXAMPLES_DIR}/pythontest.py
        ${EXAMPLES_DIR}/seawater/seawater.yml
        ${EXAMPLES_DIR}/speciation/speciation.yml
        ${EXAMPLES_DIR}/speciation/speciation.out.json
        ${EXAMPLES_DIR}/speciation/speciation.state.json
        ${EXAMPLES_DIR}/speciation_group/speciation_group.yml
        ${EXAMPLES_DIR}/speciation_group/speciation_group.out.json
        ${EXAMPLES_DIR}/speciation_group/speciation_group.state.json
        ${EXAMPLES_DIR}/speciation_swapconf/speciation_swapconf.yml
        ${EXAMPLES_DIR}/speciation_swapconf/speciation_swapconf.out.json
        ${EXAMPLES_DIR}/speciation_swapconf/speciation_swapconf.state.json
        ${EXAMPLES_DIR}/stockmayer/stockmayer.yml
        ${EXAMPLES_DIR}/stockmayer/stockmayer.out.json
        ${EXAMPLES_DIR}/stockmayer/stockmayer.state.json
        ${EXAMPLES_DIR}/swapconf/swapconf.yml
        ${EXAMPLES_DIR}/swapconf/swapconf.out.json
        ${EXAMPLES_DIR}/swapconf/swapconf.state.json
        ${EXAMPLES_DIR}/swapconf/swapconf.conformations.pqr
        ${EXAMPLES_DIR}/swapconf/swapconf.weights.dat
        ${EXAMPLES_DIR}/titration/titration.yml
        ${EXAMPLES_DIR}/titration/titration.out.json
        ${EXAMPLES_DIR}/titration/titration.state.json
        ${EXAMPLES_DIR}/calciumhydroxide/calciumhydroxide.yml
        ${EXAMPLES_DIR}/calciumhydroxide/calciumhydroxide.out.json
        ${EXAMPLES_DIR}/calciumhydroxide/calciumhydroxide.state.json
        ${EXAMPLES_DIR}/calciumhydroxide_molecular/calciumhydroxide_molecular.yml
        ${EXAMPLES_DIR}/calciumhydroxide_molecular/calciumhydroxide_molecular.out.json
        ${EXAMPLES_DIR}/calciumhydroxide_molecular/calciumhydroxide_molecular.state.json
        ${EXAMPLES_DIR}/water/ewald.yml
        DESTINATION share/faunus/examples)
