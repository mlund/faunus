#include "faunus/faunus.h"
#ifndef __cplusplus
#define __cplusplus
#endif
#ifndef CPLUSPLUS
#define CPLUSPLUS
#endif
#include "xdrfile/xdrfile_trr.h"
#include "xdrfile/xdrfile_xtc.h"

using namespace Faunus;
using namespace std;

typedef pot_minimage potT;

template<class T> class energy {
  public:
    T pair;
    energy(inputfile &in) : pair(in) {};
    double potential(vector<particle> &p, int i) {
      double phi=0;
      particle tmp=p[i];
      p[i].z=1e9;
      p[i].charge=0;
#pragma omp parallel for reduction (+:phi)
      for (int j=0; j<p.size(); j++)
        phi+=p[j].charge/sqrt(pair.sqdist(tmp,p[j]));
      p[i]=tmp;
      return pair.f*phi;
    }
};

class chargescale {
  public:
    average<double> expu;
    vector<int> sites;
    void analyze( container &c, energy<potT> pot, double dz) {
      double du=0;
      for (int i=0; i<sites.size(); i++)
        du+=pot.potential(c.p, sites[i]) * ( (c.p[sites[i]].charge+dz) - c.p[sites[i]].charge );
      for (int i=0; i<sites.size()-1; i++)
        for (int j=i+1; j<sites.size(); j++) {
          double r=sqrt( pot.pair.sqdist( c.p[i],c.p[j] ) );
          du-=c.p[i].charge*c.p[j].charge / r; // double count in pot!
          du+=dz*dz / r ;
        }
      expu+=exp(-du);
    }
    void info() {
      cout << "# Sites: ";
      for (int i=0; i<sites.size(); i++)
        cout << sites[i] << " ";
      cout << "\n# Excess chemical potential (kT) = " << -log(expu.avg()) << endl;
    }
};

int main() {
  inputfile in("g_widom.conf");
  double dz = in.getflt("dz",0.1);
  XDRFILE *xd;
  int rc,natoms_xtc,step_xtc;
  float time_xtc;
  matrix box_xtc;
  rvec *x_xtc;
  float prec_xtc = 1000.0;

  box con(100);
  atom.load("faunatoms.dat");
  energy<potT> pot(in);
  FAUrdf rdf(atom["NA"].id, atom["CL"].id,0.5,20);
  chargescale cs_one;
  cs_one.sites.push_back(21462);
  chargescale cs_two = cs_one;
  cs_two.sites.push_back(21463);

  // OPEN PQR FILE FOR NON-POS. INFO.
  iogro gro(in);
  iopqr pqr;
  ioaam aam;
  con.p = aam.load("conf.aam");
  //con.p = gro.load("conf.gro");
  //cout << con.p.size() << endl;
  //return 0;
  atom.reset_properties(con.p); // set all properties according to faunatoms.dat file!
  string s;

  // OPEN XTC TRAJECTORY FILE
  xd=xdrfile_open("coord.xtc", "r");
  if (xd!=NULL) {
    // Get length
    rc = read_xtc_natoms("coord.xtc", &natoms_xtc);
    if (rc==exdrOK && natoms_xtc==con.p.size()) {
      // Allocate memory
      x_xtc = (rvec*)calloc(natoms_xtc, sizeof(x_xtc[0]));
      cout << "# Number of atoms = " << natoms_xtc << endl;
      // Loop over frames
      while(1)
      {
        rc = read_xtc(xd, natoms_xtc, &step_xtc, &time_xtc,box_xtc, x_xtc, &prec_xtc);
        if (rc == 0) {
          con.setvolume( pow(box_xtc[0][0]*10,3) ); // store box size in container (cube assumed)
          pot.pair.setvolume(pow(con.len,3));       // ...and also in potential class.
          cout << "# Frame: " << time_xtc << "  Box = " << con.len << endl;
          for (int i=0; i<con.p.size(); i++) {
            con.p[i].x = x_xtc[i][0]*10.-con.len_half;  // store pos. in container.
            con.p[i].y = x_xtc[i][1]*10.-con.len_half;
            con.p[i].z = x_xtc[i][2]*10.-con.len_half;
          }
          // Analyse frame!!
          cs_one.analyze( con, pot, dz);
          cs_two.analyze( con, pot, dz);
          rdf.update(con);
        }
        else break;
      }
    }
    xdrfile_close(xd);
  }
  rdf.write("gofr.dat");
  cs_one.info();
  cs_two.info();
};
