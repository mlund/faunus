include_directories($faunus_source_dir/classes)
set(objs inputfile io titrate point physconst slump container
    potentials hardsphere group particles analysis species) 
set_source_files_properties(${objs} PROPERTIES LANGUAGE CXX)
add_library (faunus SHARED ${objs})

