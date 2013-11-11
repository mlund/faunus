fau_example(jmhenriques-peptide_bulk "." peptide.cpp)
fau_example(jmhenriques-peptide_surface "." peptide.cpp)
set_target_properties(jmhenriques-peptide_surface PROPERTIES COMPILE_DEFINITIONS "SLIT")
