include_directories(${CMAKE_SOURCE_DIR}/include/xdrfile)
set(xdrobjs xdrfile xdrfile_xtc xdrfile_trr)
set_source_files_properties( ${xdrobjs} PROPERTIES LANGUAGE C)
add_library(xdrfile STATIC ${xdrobjs})
