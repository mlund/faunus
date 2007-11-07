#!/bin/bash

file=$1
res=$2

function getph() {
  ph=`cat ${file} | grep "pH=" | sed -e 's/pH=//g' | awk '{print $1}'`
}

function getcharges() {
  z=`cat ${file} | grep ${res} | sed -e 's/'${res}': //g'`
}

getph
getcharges

echo "$ph $z"
