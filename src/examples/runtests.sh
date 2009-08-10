#!/bin/bash

# Arguments
# $1 = directory with test files
# $2 = name of binary program to execute
function run {
if [ -x $1/../$2 ]
then
  echo -e "\033[1mTesting $2.\033[0m"
  tput sgr0
  cd $1
  rm -f confout.aam
  ../$2 test.conf #1>/dev/null
  cd - > /dev/null
fi
}

exdir="./"

run ${exdir}/titration/test1 pka
run ${exdir}/twobody/test1 twobody

exit 0
