include_directories(${CMAKE_SOURCE_DIR}/include/xdrfile)
set(xdrobjs xdrfile xdrfile_xtc xdrfile_trr)
set_source_files_properties( ${xdrobjs} PROPERTIES LANGUAGE C)
add_library( xdrfile STATIC ${xdrobjs} )
set_target_properties(xdrfile PROPERTIES COMPILE_FLAGS "-fPIC" )

# OpenSolaris dislikes the RPC/XDR implementation in xdrfile.
if (CMAKE_HOST_SYSTEM_NAME MATCHES "SunOS")
  add_library( rpcsvc SHARED IMPORTED )
  set_target_properties( rpcsvc PROPERTIES IMPORTED_LOCATION /usr/lib/librpcsvc.so )
  target_link_libraries( xdrfile rpcsvc )
  add_definitions( -DHAVE_RPC_XDR_H )
endif (CMAKE_HOST_SYSTEM_NAME MATCHES "SunOS")
