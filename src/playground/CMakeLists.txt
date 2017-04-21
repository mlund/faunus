# Bosse's playground
fau_example(boj-grand "./bosse" grand.cpp)

# Bearnix playground
fau_example(bearnix-cyl "./bearnix" cyl.cpp)

# Carolina's playground
fau_example(cc-polymers "./cragnell" polymers.cpp)

# Fernando's playground
fau_example(fer-cyl "./fernando" cyl.cpp)

# Mikael's playground
fau_example(mlund-volume "./mlund" volume.cpp)
fau_example(mlund-cosan "./mlund" cosan.cpp)
fau_example(mlund-angular "./mlund" angular.cpp)
fau_example(mlund-rotor "./mlund" twobodyrotor.cpp)
fau_example(mlund-test "./mlund" test.cpp)
fau_example(mlund-cluster "./mlund" cluster.cpp)
fau_example(mlund-manybody "./mlund" manybody.cpp)
fau_example(mlund-membrane "./mlund" membrane.cpp)
fau_example(mlund-membrane-penalty "./mlund" membrane.cpp)
set_target_properties( mlund-test PROPERTIES COMPILE_DEFINITIONS "SOMETHING")
set_target_properties( mlund-membrane-penalty PROPERTIES COMPILE_DEFINITIONS "PENALTY")
if(ENABLE_MPI AND MPI_CXX_FOUND)
  fau_example(mlund-membrane-temper "./mlund" membrane-temper.cpp)
  fau_example(mlund-rennes "./mlund/rennes" twobody.cpp)
endif()

# Anil's playground
fau_example(akurut-manybody "./akurut" manybody.cpp)
fau_example(akurut-gouychapman "./akurut" gouychapman.cpp)
fau_example(akurut-gouychapman-slit "./akurut" gouychapman.cpp)
set_target_properties( akurut-gouychapman-slit PROPERTIES COMPILE_DEFINITIONS "SLIT")

# Cedric's playground
fau_example(cedric-polymers "./cedric" polymers.cpp)

# Robert's playground
fau_example(vacha-gouychapman "./vacha" gouychapman.cpp)
fau_example(vacha-gouychapman-temper "./vacha" gouychapman.cpp)
set_target_properties( vacha-gouychapman-temper
  PROPERTIES COMPILE_DEFINITIONS "TEMPER")

fau_example(cigar "./vacha" cigar.cpp)
fau_example(sc-npt "./vacha" sc-npt.cpp)
fau_example(psc-nvt "./vacha" psc-nvt.cpp)
fau_example(test2 "./vacha" test2.cpp)

# Stenqvist's playground
fau_example(stenqvist-nemo "./stenqvist" nemo.cpp)
set_target_properties(stenqvist-nemo PROPERTIES OUTPUT_NAME "nemo")

# Axel's playground
fau_example(axel-clustertest "./axel/clustertest" cluster.cpp)
set_target_properties(axel-clustertest PROPERTIES OUTPUT_NAME "clustertest")

fau_example(axel-tempertest "./axel/tempertest" temper.cpp)
set_target_properties(axel-tempertest PROPERTIES OUTPUT_NAME "tempertest")

# Alexei's playground
fau_example(alexei-hs "./alexei" HS.cpp)
set_target_properties(alexei-hs PROPERTIES OUTPUT_NAME "HS")
fau_example(alexei_cyl "./alexei" cyl.cpp)
set_target_properties(alexei_cyl PROPERTIES OUTPUT_NAME "cyl")

# Jana's playground
fau_example(jana-wall "./jana" wall.cpp)

if(ENABLE_POWERSASA)
  fau_example(heyda-polymers "./heyda" polymers.cpp)
endif()

if (EXISTS ${MYPLAYGROUND})
  add_subdirectory(${MYPLAYGROUND} ${MYPLAYGROUND})
endif()
