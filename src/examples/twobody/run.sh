function mkinput() {
echo "macrosteps 10
microsteps $microsteps
cellradius $cellradius
bjerrum    7.1
LJeps      0.3
nion1      $nion1
nion2      $nion2
tion1      NA
tion2      CL
nprot1     1
protein1   lysozyme.aam
nprot2     1
protein2   lysozyme.aam
atomfile   ../../../misc/faunatoms.dat
gofr_pp_rf 0.95
pH         $pH
dp_salt    90
testsuite_testfile test.stable
testsuite_stable   yes
" > twobody.conf
}

exe="./twobody"
cellradius=100

for salt in 250 64
do
  for pH in 9 7 5
  do
    # Equilibrate
    rm -f confout.aam rdf*.dat
    prefix="salt${salt}.pH${pH}"
    nion1=$[salt+26]
    nion2=$[salt]
    microsteps=1000
    mkinput
    $exe

    # Production run
    microsteps=200000
    mkinput
    $exe > $prefix.out
    mv rdfprot.dat $prefix.rdfprot.dat
    mv rdfatomic.dat $prefix.rdfatomic.dat
    mv rdfsalt.sat $prefix.rdfsalt.dat
    mv distributions.dat $prefix.dists.dat
    mv confout.pqr $prefix.pqr
    mv coord.xtc $prefix.xtc
    mv confout.aam $prefix.aam
  done
done
