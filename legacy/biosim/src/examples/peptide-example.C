#include "../peptide.h"
#include<iostream>

using namespace std;

int main() {

  int i=aminoacid::NA;
  int j=aminoacid::CL;

  aminoacid aa;

  aa.loadpmf("/Users/mikael/biosim/pmf/NA-CL.dat");

  aa.pairpot[i][j].list();

  cout << "# " << aa.pairpot[i][j].xmin << " "
       << aa.pairpot[i][j].xmax << " "
       << aa.pairpot[i][j].res << endl;

  cout << aa.getId("GHOST");
};
