#!/bin/bash

host=`hostname -s`
mac=`uname -s`
rel=`uname -r`
sys=`uname -p`
date=`date "+%d%m%y-%H%M"`
#date=`date "+%d.%m.%y"`
id="$HOSTNAME"
suffix="ERROR"
if [ -z $CC ]; then
  CC=cc
fi
if [ -z $CXX ]; then
  CXX=g++
fi

log="faunlog.${host}-${CXX}_${sys}-${mac}${rel}"

if [ -e Makefile ]; then
  make clean
fi
if [ -e CMakeCache.txt ]; then
  rm CMakeCache.txt
fi

if [ "${1}" = "force" ]; then
  svnold=0;
else
  svnold=`svn info | grep "Revision:" | gawk '{print $2}'`
fi
svn update > /dev/null
svnnew=`svn info | grep "Revision:" | gawk '{print $2}'`

if [ ${svnnew} -gt ${svnold} ]; then
  log="${log}-svn${svnnew}"
  echo "Subversion revision: ${svnnew}" > ${log} 2>&1
  cmake -DENABLE_BABEL=off . >> ${log} 2>&1
  if [ $? = 0 ]; then
    make libfaunus pka twobody widom_cube >> ${log} 2>&1
    if [ $? = 0 ]; then
      suffix="buildOK"
      tbeg=`date "+%s"`
      make test >> ${log} 2>&1
      if [ $? = 0 ]; then
        suffix="${suffix}.testOK"
      else
        suffix="${suffix}.testERROR"
      fi
      testtime=$[`date "+%s"`-tbeg]
    fi
  fi
  mv ${log} ${log}.${suffix}.${testtime}secs
fi
