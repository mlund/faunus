#include "chain.h"

/*
 * CHAIN CONSTRUCTOR 
 *
 * input = String of residues that defines the chain, i.e. "3 CYS LYS NTR"
 *         The first word must a number, defining the chain length.
 * gp    = Graft point in particle vector (default is -1 == not grafted)
 *
 */
chain::chain(species &spc, string input, int gp) {
  int len,i=0;
  istringstream in(input);
  string type;
  graftpoint=gp;

  in >> len;
  v.resize(len);
  
  while (i<len) {
    in >> type;
    v[i].id=spc.id(type);
    v[i].radius=spc.d[v[i].id].radius;
    v[i].charge=spc.d[v[i].id].charge;
    i++;
  };
};

