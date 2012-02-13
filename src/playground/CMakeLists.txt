# Mikael's playground
fau_example(mlund-cosan "./mlund" cosan.cpp)
fau_example(mlund-test "./mlund" test.cpp)
fau_example(mlund-cluster "./mlund" cluster.cpp)
fau_example(mlund-manybody "./mlund" manybody.cpp)
set_target_properties( mlund-test PROPERTIES COMPILE_DEFINITIONS "SOMETHING")

