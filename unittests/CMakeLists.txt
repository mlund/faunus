# ---------------------------
#   Add unittests
# ---------------------------
add_test( titration_implicitsalt titration_implicitsalt/test.run )
add_test( titration_grandcanonical titration_grandcanonical/test.run )
add_test( titration_simple titration_simple/test.run )
add_test( isobaric_minimg_rigidmolecules isobaric_minimg_rigidmolecules/test.run )
add_test( twobody twobody/test.run )
add_test( bulk_electrolyte_sphere bulk_electrolyte_sphere/test.run )
add_test( bulk_electrolyte_minimumimage bulk_electrolyte_minimumimage/test.run )
add_test( bulk_electrolyte_hypersphere bulk_electrolyte_hypersphere/1.run )

