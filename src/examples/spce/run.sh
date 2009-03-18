#!/bin/bash

#--- Write config file ---
function mkinput() {
echo "
macrosteps       ${macrosteps}
microsteps       ${microsteps}
bjerrum          ${bjerrum}
nion1            ${nion1}
nion2            ${nion2}
tion1            ${tion1}
tion2            ${tion2}
protein          ${protein}
atomfile         ${atomfile}
cellradius       ${cellradius}
boxlen           ${boxlen}
distmax          ${distmax}
mrtdpr           ${mrtdpr}
mrtdpt           ${mrtdpt}
mtdp             ${mtdp}
mrdp             ${mrdp}
multtdp          ${multtdp}
multrdp          ${multrdp}
multinum         ${multinum}
t                ${t}
rot              ${rot}
tr               ${tr}
multi            ${multi}
epso             ${epso}
rfield_cavity    ${rfield_cavity}
rfield_epso      ${rfield_epso}
rfield_epsi      ${rfield_epsi}
rfield_bjerrum   ${rfield_bjerrum}
rfield_steps     ${rfield_steps}
splash           ${splash}
" > spce.conf
}
function move(d) {
  mv systemenergy.dat ${d}sys.dat
  mv dipoledist.dat   ${d}dipdist.dat
  mv dipolexdist.dat  ${d}dipxdist.dat
  mv rdf-cell-OW.dat  ${d}rdf-c-OW.dat
  mv OW-OW.dat        ${d}OW-OW.dat
  mv kirkwood.dat     ${d}kirkwood.dat
  mv mu-OW.dat        ${d}mu-OW.dat
}

#--- Input parameters ---

macrosteps=10
microsteps=2000
bjerrum=560.2
nion1=0
nion2=0
tion1="NA"
tion2="CL"
protein="calbindin.aam"
atomfile="../../../misc/faunatoms.dat"
#geometry
cellradius=30
boxlen=60
distmax=30
#markov parameters
mrtdpr=1.0
mrtdpt=0.35
mtdp=0.7
mrdp=2.0
multinum=4
multtdp=0.5
multrdp=0.5
t=1
rot=4
tr=4
multi=1
#ionisation analysis
epso=80
rfield_cavity=5.0
rfield_epso=80.0
rfield_epsi=1.0
rfield_bjerrum=7.0025
rfield_steps=100

splash="no"
for epso in 80  
do
  suffix="test"
  microsteps=1000
  mkinput
  ./spce > ${suffix}
  microsteps=1000
#  mkinput
#  ./spce > RF/${suffix}
#  mv rdf-OW-OW.dat RF/${suffix}rdf-OW-OW.dat
#  mv rdf-cell-OW.dat RF/${suffix}rdf-cell-OW.dat
#  mv mu-OW.dat RF/${suffix}mu-OW.dat
#  mv kirkwood.dat RF/${suffix}kirkwood.dat
#  mv confout.gro RF/${suffix}confout.gro
#  mv systemenergy.dat RF/${suffix}syse.dat
#  mv dipoledist.dat RF/${suffix}dipoledist.dat
#  mv widom.dat RF/${suffix}widom.dat
#  cp confout.aam RF/${suffix}conf.aam
done

