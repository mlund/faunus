#!/bin/bash

# A simple configuration file.
# Lines containing [ or # or blank
# lines will be ignored# lines will be ignored
                       
# System               # System
function bottleInput() {
echo "macrosteps            $macrosteps        
microsteps            $microsteps        
boxlen                $boxlen            
polymer               $polymer           
N_polymer             $N_polymer         
iw2_21                $iw2_21            
iw2_31                $iw2_31            
iw2_23                $iw2_23            
wd_21                 $wd_21             
wd_31                 $wd_31             
wd_23                 $wd_23             
trans_f               $trans_f           
moltrans_dp           $moltrans_dp       
clt_f                 $clt_f             
rep_f                 $rep_f             
crank_f               $crank_f           
crankshaft_dp         $crankshaft_dp     
crankshaft_max        $crankshaft_max    
crankshaft_min        $crankshaft_min    
rot_f                 $rot_f             
molrot_dp             $molrot_dp         
tit_f                 $tit_f             
Type_1                $Type_1            
Type_2                $Type_2            
pc                    $pc                
pK                    $pK                
cluster_def           $cluster_def       
bjerrum               $bjerrum           
atomfile              $atomfile          
" > $id.conf
}

#harmonic_k            $harmonic_k       
#harmonic_req          $harmonic_req
     
function controllerInput {
echo "
macrosteps $macro
microsteps $micro
temper     $temper
" > replica_chains.conf
}

macrosteps=10  
microsteps=10        
boxlen=1000            
polymer="../sq_well/9-bead.mol2"           
N_polymer=50         
iw2_21=870.25            
iw2_31=870.25            
iw2_23=121            
wd_21=0             
wd_31=0             
wd_23=-1.5             
trans_f=1          
moltrans_dp=200       
clt_f=1             
rep_f=1             
crank_f=1           
crankshaft_dp=3     
crankshaft_max=1    
crankshaft_min=3    
rot_f=1             
molrot_dp=3         
tit_f=1             
Type_1="B"            
Type_2="C"            
pc=0                
pK=0                
cluster_def=3       
bjerrum=7.1          
atomfile=../sq_well/faunatoms.dat          

# 1. Generate replica input files
for id in a b c d e f
do
  if [ "$id" == "a" ]; then wd_23=-0.7; fi
  if [ "$id" == "b" ]; then wd_23=-0.9; fi
  if [ "$id" == "c" ]; then wd_23=-1.1; fi
  if [ "$id" == "d" ]; then wd_23=-1.3; fi
  if [ "$id" == "e" ]; then wd_23=-1.4; fi
  if [ "$id" == "f" ]; then wd_23=-1.5; fi
  bottleInput
done

export OMP_NUM_THREADS=2

# 2. Eq run
rm -f *.dump
macro=5
micro=10
temper="yes"
controllerInput
./replica_chains
exit

# 3. Production
micro=10000
temper="yes"
controllerInput
./replica

