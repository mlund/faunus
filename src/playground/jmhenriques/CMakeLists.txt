fau_example(peptide_bulk "." peptide.cpp)
fau_example(peptide_surface "." peptide.cpp)
set_target_properties(peptide_surface PROPERTIES COMPILE_DEFINITIONS "SLIT")
