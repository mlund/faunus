# Mikael's playground
fau_example(mlund-cosan "./mlund" cosan.cpp)
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

if(ENABLE_MPI AND MPI_CXX_FOUND)
#  fau_example(mpitest "./mlund" mpitest.cpp)
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

fau_cigarexample(cigar "./vacha" cigar.cpp)
set_target_properties( cigar PROPERTIES COMPILE_DEFINITIONS "CIGARPARTICLE")

fau_cigarexample(cigar-npt "./vacha" cigar-npt.cpp)
set_target_properties( cigar-npt PROPERTIES COMPILE_DEFINITIONS "CIGARPARTICLE")

