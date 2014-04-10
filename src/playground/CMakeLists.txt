#----- function to add a generic C++ example (CIGARPARTICLES) -----
function( fau_dipoleexample tname tdir tsrc )
  add_executable( ${tname} "${tdir}/${tsrc}" )
  set_source_files_properties( "${tdir}/${tsrc}" PROPERTIES LANGUAGE CXX)
  set_target_properties(${tname}
    PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/${tdir}/"
    COMPILE_DEFINITIONS "DIPOLEPARTICLE"
    EXCLUDE_FROM_ALL TRUE )
  target_link_libraries(${tname} libdipolefaunus)
  install (DIRECTORY "${tdir}"
           DESTINATION "share/faunus/examples"
           PATTERN ".svn" EXCLUDE)
endfunction( fau_dipoleexample )

# Bosse's playground
fau_example(boj-grand "./bosse" grand.cpp)

# Bearnix playground
fau_example(bearnix-cyl "./bearnix" cyl.cpp)

# Carolina's playground
fau_example(cc-polymers "./cragnell" polymers.cpp)

# Carolina's playground
fau_example(fer-cyl "./fernando" cyl.cpp)

# Mikael's playground
fau_example(mlund-volume "./mlund" volume.cpp)
fau_example(mlund-cosan "./mlund" cosan.cpp)
fau_example(mlund-angular "./mlund" angular.cpp)
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

fau_cigarexample(sc-npt "./vacha" sc-npt.cpp)
set_target_properties( sc-npt PROPERTIES COMPILE_DEFINITIONS "CIGARPARTICLE")

fau_cigarexample(psc-nvt "./vacha" psc-nvt.cpp)
set_target_properties( psc-nvt PROPERTIES COMPILE_DEFINITIONS "CIGARPARTICLE")

fau_cigarexample(test2 "./vacha" test2.cpp)
set_target_properties( psc-nvt PROPERTIES COMPILE_DEFINITIONS "CIGARPARTICLE")

# Stenqvist playground
fau_dipoleexample(stenqvist-nemo "./stenqvist" nemo.cpp)
set_target_properties(stenqvist-nemo PROPERTIES OUTPUT_NAME "nemo")

# Axels playground
fau_example(axel-clustertest "./axel/clustertest" cluster.cpp)
set_target_properties(axel-clustertest PROPERTIES OUTPUT_NAME "clustertest")

fau_example(axel-tempertest "./axel/tempertest" temper.cpp)
set_target_properties(axel-tempertest PROPERTIES OUTPUT_NAME "tempertest")

# Alexei's playground
fau_example(alexei-hs "./alexei" HS.cpp)
set_target_properties(alexei-hs PROPERTIES OUTPUT_NAME "HS")
#fau_dipoleexample(alexei_stockmayer "./alexei" stockmayer.cpp)
#fau_dipoleexample(alexei_stockmayer_polarizable "./alexei" stockmayer.cpp)
#set_target_properties(alexei_stockmayer PROPERTIES OUTPUT_NAME "stockmayer")
#set_target_properties(alexei_stockmayer_polarizable PROPERTIES OUTPUT_NAME "stockmayer_polarizable")
#set_target_properties(alexei_stockmayer_polarizable PROPERTIES COMPILE_DEFINITIONS "POLARIZE;DIPOLEPARTICLE")

# Jana's playground
fau_example(jana-wall "./jana" wall.cpp)

if (EXISTS ${MYPLAYGROUND})
  add_subdirectory(${MYPLAYGROUND} ${MYPLAYGROUND})
endif()
