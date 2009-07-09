/*! \file Coarse grain proteins 
 *  One big sphere with retained charge loci 
 *
 */
#include "faunus/faunus.h"
#include "faunus/energy.h"
#include "faunus/potentials/pot_minimage.h"
//#include "faunus/moves/markovmove.h"

using namespace Faunus;
using namespace std;

int main() {
  inputfile in("cg1.conf");
  cell con(in.getflt("radius"));         // We want a spherical cell
  //Varius objects and variables 
  macromolecule protein;
  ioaam aam;
  iopqr pqr;
  protein.add(con, aam.load(in.getstr("protein")));
  protein.center(con);
  protein.accept(con);
  hardsphere hs;
  particle probe;
  probe.radius=0;
  particle blob;
  blob.mw=0.1;
  int hit=0;
  int hitcgp=0;
  int comhit=0;
  bool cgpb, pb;
  vector<particle> p;
  macromolecule cgp;
  p.clear();
  //Compute the center of positive and negative charge
  int s;
  s=con.p.size(); 
  for (int i=0; i<s; i++) {
    if (con.p[i].charge != 0.0 ) { //&& std::abs(con.p[i].charge) >0.05) {
      p.push_back(con.p[i]);
    }
  }
  p.push_back(blob);
  con.p.push_back(blob);
  s=p.size()-1;
  cgp.beg=0, cgp.end=s;
  //Adjust the volume
  p[s].radius=in.getflt("blobradius");
    for (int macro=0; macro<in.getflt("macro"); macro++) {        
      for (int micro=0; micro<in.getflt("micro"); micro++) {
        cgpb=pb=false;
        con.randompos(probe);
        if(hs.overlap(con.p, protein, probe)==true) {
          hit++;
          pb=true;
        }
        if(hs.overlap(p, cgp, probe)==true) {
          hitcgp++;
          cgpb=true;
        }
        if((cgpb && pb)==true)
          comhit++;
      }
    }

  aam.save("CG2-R"+in.getstr("blobradius")+"-"+in.getstr("protein"), p);
  pqr.save("comb.pqr", con.p);

  cout << con.info()<<endl <<protein.info()<<in.info()<<endl
       << "#  Original protein radius       = "<< protein.radius(con.p) <<endl
       << "#  Coarse grained protein radius = "<< cgp.radius(p)<<endl
       << "#  Number of shots               = "<< in.getflt("macro")*in.getflt("micro") <<endl
       << "#  Number of hist in old         = "<< hit << endl
       << "#  Number of hist in coarse grain= "<< hitcgp << endl
       << "#  Number of common hits         = "<< comhit << endl
       << "#  Hit ratio        (old)        = "<< hit/double(in.getflt("macro")*in.getflt("micro")) <<endl
       << "#  Hit ratio        (new)        = "<< hitcgp/double(in.getflt("macro")*in.getflt("micro")) <<endl
       << "#  Hit ratio        (common/old) = "<< comhit/double(hit) <<endl
       << "#  Hit ratio        (common/new) = "<< comhit/double(hitcgp) <<endl
       << "#  _____________________________________" << endl
       << "#  Protein volume is " << con.getvolume()*hit/(in.getflt("macro")*in.getflt("micro"))<<" A^3"<<endl;
   
}

