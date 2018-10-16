#!/bin/bash

function mkinput() {
echo '
atomlist: # define all atom types below
    - Na: {q:  0.0, eps: 1.0, sigma: 4.0, dp: 40}
    - Ca: {q:  2.0, eps: 1.0, sigma: 4.0, dp: 40}
    - K: {q:  0.0, eps: 1.0, sigma: 4.0, dp: 40}
    - OH: {q: -1.0, eps: 1.0, sigma: 4.0, dp: 10}
    - H4SiO4: {q: 0.0, eps: 1.0, sigma: 4.0, dp: 40}
    - H3SiO4: {q: -1.0, eps: 1.0, sigma: 4.0, dp: 40}
    - H2SiO4: {q: -2.0, eps: 1.0, sigma: 4.0, dp: 40}
    - ID1 : {q: 0.0, eps: 1.0, sigma: 2.0, dp: 20}
    - ID2 : {q: 1.0, eps: 1.0, sigma: 2.0, dp: 20}
    - ID3 : {q: -2.0, eps: 1.0, sigma: 2.0, dp: 20}
    - ID4 : {q: -1.0, eps: 1.0, sigma: 2.0, dp: 20}

moleculelist:
    - OH_m: {atoms: [OH], atomic: true}
    - Ca_m: {atoms: [Ca], atomic: true}
    - H4SiO4_m: {atoms: [H4SiO4], atomic: true}
    - H3SiO4_m: {atoms: [H3SiO4], atomic: true}
    - H2SiO4_m: {atoms: [H2SiO4], atomic: true}
    - CaOH2: {atoms: [Ca, OH, OH], atomic: true}
#    - CSH1_m: {atoms: [ID1], atomic: true}
#    - CSH08_m: {atoms: [ID2], atomic: true}
#    - K_m: {atoms: [K], atomic: true}
#    - N_sol: {atoms: [Natm], atomic: true}
#    - NaCl: {Ninit: 100, atoms: [Na,Cl], atomic: true, Ninactive: 0}
#    - ID1_m: {atoms: [ID1], atomic: true}
#    - ID2_m: {atoms: [ID2], atomic: true}
#    - ID3_m: {atoms: [ID3], atomic: true}
#    - ID4_m: {atoms: [ID4], atomic: true}

insertmolecules:
    - OH_m: {N: 29494 }
    - Ca_m: {N: 14747, inactive: false }
    - H4SiO4_m: {N: 500, inactive: false }
    - H3SiO4_m: {N: 500, inactive: true } #481
    - H2SiO4_m: {N: 500, inactive: true }
    - CaOH2: {N: 1, inactive: true}
#    - CSH1_m: {N: 200, inactive: true}
#    - CSH08_m: {N: 200, inactive: true}
#    - Cl_m: {N: 100, inactive: false }
#    - N_sol: {N: 10 }
#    - ID1_m: {N: 100, inactive: false}
#    - ID2_m: {N: 10, inactive: false}
#    - ID3_m: {N: 100, inactive: true}
#    - ID4_m: {N: 10}
reactionlist:
    - S1: {log_k: 1.5, reactants: [], products: [H4SiO4_m], canonic: true, N_reservoir: 0}
#    - R2: {log_k: -1.0, reactants: [], products: [N_sol], canonic: true, N: 1000, N_reservoir: 0}
#    - R3: {log_k: 0.0, reactants: [Na_m, Na_m], products: [Ca_m], canonic: false}
#    - R4: {log_k: 0.0, products: [Na_m, Na_m], reactants: [Ca_m], canonic: false}
    - T1 :  {log_k: 4.16, reactants: [H4SiO4_m, OH_m], products: [H3SiO4_m], N_reservoir: 0, canonic: false}
#- T1 :  {log_k: 4.16, reactants: [OH_m], products: [H3SiO4_m], N_reservoir: 0, canonic: false}
    - T2 :  {log_k: 0.8, reactants: [H3SiO4_m, OH_m], products: [H2SiO4_m], N_reservoir: 0, canonic: false}
    - CS1: {log_k: 13.078, reactants: [], products: [Ca_m, OH_m, OH_m, H4SiO4_m], canonic: true, N_reservoir: 0}
    - CS08: {log_k: 57.595, reactants: [], products: [Ca_m, Ca_m, Ca_m, Ca_m, OH_m, OH_m, OH_m, OH_m, OH_m, OH_m, OH_m, OH_m,
             H4SiO4_m, H4SiO4_m, H4SiO4_m, H4SiO4_m, H4SiO4_m], canonic: true, N_reservoir: 0}
# Ideal test reactions devised for 2000 box length
#    - Test1R1: {log_k: -1.983, reactants: [], products: [ID1_m], N_reservoir: 0, canonic: true} #     - ID1_m: {N: 100, inactive: true}
#    - Test2R1: {log_k: -3.966, reactants: [], products: [ID1_m, ID1_m], N_reservoir: 0, canonic: true} #     - ID1_m: {N: 100, inactive: true}
#    - Test1: {log_k: 5.34, reactants: [ID1_m, ID4_m], products: [ID2_m], N_reservoir: 0, canonic: false}
#    - Test2: {log_k: 3.43, reactants: [ID2_m, ID4_m], products: [ID3_m], N_reservoir: 0, canonic: false}

energy:
    - nonbonded:
        default:
            - wca: {mixing: LB}
            - coulomb: {type: plain, epsr: 80, cutoff: 100}

moves:
    - transrot: {molecule: Ca_m, dp: 40, dprot: 0, repeat: 3, dir: [1,1,1] }
    - transrot: {molecule: OH_m, dp: 40, dprot: 0, repeat: 3 }
    - transrot: {molecule: H4SiO4_m, dp: 40, dprot: 0, repeat: 3 }
    - transrot: {molecule: H3SiO4_m, dp: 40, dprot: 0, repeat: 3 }
    - transrot: {molecule: H2SiO4_m, dp: 40, dprot: 0, repeat: 3 }
#    - transrot: {molecule: ID1_m, dp: 40, dprot: 0, repeat: 3 }
#    - transrot: {molecule: ID2_m, dp: 40, dprot: 0, repeat: 3 }
#    - transrot: {molecule: ID3_m, dp: 40, dprot: 0, repeat: 3 }
#    - transrot: {molecule: ID4_m, dp: 40, dprot: 0, repeat: 3 }
#    - volume: {dV: 0.1}
    - speciation: {repeat: 10}
analysis:
    - systemenergy: {file: energy.dat, nstep: 1000}
    - savestate: {file: confout.pqr}
    - savestate: {file: confout.state}
    - density: {nstep: 10}
#    - widom: {molecule: CaOH2, nstep: 50, ninsert: 10}
mcloop: {macro: 10, micro: 30}
geometry: {length: 2000} # sidelength(s): number OR array
temperature: 300
random: {seed: hardware }
' > input.yml
}

# $pH
for pH in 7.21; do
  #rm -f confout.state
  mkinput
  ../yason.py input.yml > input.json
  ../faunus -i input.json -o output.json -s confout.state
  ../yason.py output.json > output.yml
done
