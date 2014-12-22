include_directories(${CMAKE_SOURCE_DIR}/include)

#----- examples where openbabel is needed -----
if (ENABLE_BABEL AND OPENBABEL_FOUND)
  include_directories(${OPENBABEL_INCLUDE_DIR})
  #fau_example(babel2aam tools babel2aam.cpp)
endif (ENABLE_BABEL AND OPENBABEL_FOUND)

