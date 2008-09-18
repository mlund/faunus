set(objs pyfaunus)
set_source_files_properties(${objs} PROPERTIES LANGUAGE CXX)

add_library(pyfaunus SHARED ${objs})

include_directories("/usr/include/python2.5")
include_directories("../../include")
target_link_libraries(pyfaunus boost_python libfaunus)
link_directories("/usr/lib/python2.5")              # Not needed?

set_target_properties(pyfaunus PROPERTIES PREFIX "")

