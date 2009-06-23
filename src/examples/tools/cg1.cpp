/*! \file Coarse grain proteins 
 *  
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
  particle plus,minus;
  plus.x=plus.y=plus.z=0;
  minus.x=minus.y=minus.z=0;
  double   plussum, minussum;
  plussum=minussum=0;
  int hit=0;
  int hitcgp=0;
  int comhit=0;
  bool cgpb, pb;
  vector<particle> p;
  macromolecule cgp;
  p.clear();
  cgp.beg=0, cgp.end=1;
  //Compute the center of positive and negative charge
  int s;
  s=con.p.size(); 
  for (int i=0; i<s; i++) {
    if (con.p[i].charge > 0.0) {
      plus.x  +=con.p[i].x*con.p[i].charge;
      plus.y  +=con.p[i].y*con.p[i].charge;
      plus.z  +=con.p[i].z*con.p[i].charge;
      plussum +=con.p[i].charge;
    }
    if (con.p[i].charge < 0.0) {
      minus.x +=con.p[i].x*con.p[i].charge;
      minus.y +=con.p[i].y*con.p[i].charge;
      minus.z +=con.p[i].z*con.p[i].charge;
      minussum+=con.p[i].charge;
    }
  }
  plus.x=plus.x/plussum;
  plus.y=plus.y/plussum;
  plus.z=plus.z/plussum;
  plus.charge=plussum;
  p.push_back(plus);
  con.p.push_back(plus);
  minus.x=minus.x/minussum;
  minus.y=minus.y/minussum;
  minus.z=minus.z/minussum;
  minus.charge=minussum;
  p.push_back(minus);
  con.p.push_back(minus);
  //Adjust the volume
  p[0].radius=in.getflt("blobradius");
  p[1].radius=p[0].radius;

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

  aam.save("CG-"+in.getstr("protein"), p);
  pqr.save("comb.pqr", con.p);

  cout << con.info()<<endl <<protein.info()<<endl
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

