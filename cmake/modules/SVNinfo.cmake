#
# Cmake script to extract the Subversion revision
# of the current directory. Zero is returned if the 
# dir is not under version control.
# Mikael Lund, 2009
#
find_package(Subversion)
if (Subversion_FOUND)
  set(SVNREVISION 0)
  execute_process(
    COMMAND ${Subversion_SVN_EXECUTABLE} info
    RESULT_VARIABLE SVNRC
    OUTPUT_VARIABLE SVNOUT
    ERROR_QUIET )
  if(SVNRC EQUAL 0)
    execute_process(
      COMMAND echo "${SVNOUT}"
      COMMAND grep "Revision: "
      COMMAND awk "{print $2}"
      OUTPUT_VARIABLE SVNREVISION
      ERROR_QUIET )
  endif(SVNRC EQUAL 0)
endif (Subversion_FOUND)
