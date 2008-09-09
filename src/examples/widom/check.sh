#!/bin/bash

error="1"  # in percent
excess_ref="-0.65823"
tmp=".tmp"

echo -n "Widom test. Running. "
rm -f $tmp
./widom > $tmp
excess_cal="`cat ${tmp} | grep "Excess ch" | gawk '{print $6}'`"
#echo  $excess_cal $excess_ref

echo "import sys
import math
if (abs(${excess_ref}-${excess_cal})/${excess_ref})*100 > ${error}:
  print 'FAILED!'
  sys.exit(1)
print 'Passed.'
sys.exit(0)
" > ${tmp}
python ${tmp}
rm -f ${tmp}
