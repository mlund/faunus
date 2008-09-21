set(objs xtc xdrf smalloc fatal gmxfio futil) 
add_library (libxtc STATIC ${objs})
set_target_properties(libxtc PROPERTIES OUTPUT_NAME xtc)
install(FILES "xtc.h" DESTINATION include/faunus)
