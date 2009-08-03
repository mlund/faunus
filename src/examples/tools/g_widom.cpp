#include "faunus/faunus.h"
#include "faunus/average_vec.h"
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

class gmxenergy {
  private:
    Faunus::io inout;
  public:
    vector<double> u; // energy in kT
    gmxenergy(string file) {
      double kJ_per_kT=pyc.kT2kJ(1.0);
      vector<string> s;
      inout.readfile(file, s);
      u.resize(s.size());
      for (int i=0; i<s.size(); i++)
        u[i]=atof(s[i].c_str() ) / kJ_per_kT;
    }
};

template<class T> class energy {
  public:
    T pair;
    energy(inputfile &in) : pair(in) {};

    double potential(vector<particle> &p, int i) {
      int n=p.size();
      double phi=0;
      particle tmp=p[i];
      p[i].z=1e9;
      p[i].charge=0;
//#pragma omp parallel for reduction (+:phi)
      for (int j=0; j<n; j++)
        phi+=p[j].charge/sqrt(pair.sqdist(tmp,p[j]));
      p[i]=tmp;
      return pair.f*phi;
    }

    double system(vector<particle> &p) {
      int n=p.size();
      double uu=0;
#pragma omp parallel for reduction (+:uu) schedule (dynamic)
      for (int i=0; i<n-1; i++)
        for (int j=i+1; j<n; j++)
          uu += p[i].charge*p[j].charge/sqrt(pair.sqdist(p[i],p[j]));
      return pair.f*uu;
    }
};

class chargescale {
  private:
    unsigned int cnt, ucnt;
  public:
    average_vec<long double> expu, u, mu, a1,a2,a3;
    vector<int> sites;
    vector<double> uvec;
    chargescale() {
      cnt=ucnt=0;
      mu.maxcnt=1;
      u.maxcnt=500;
    }
    void analyze( container &c, energy<potT> pot, double dq) {
      cnt++;
      double du=0;
      for (int i=0; i<sites.size(); i++)
        du += pot.potential(c.p, sites[i]) * dq;
      for (int i=0; i<sites.size()-1; i++)
        for (int j=i+1; j<sites.size(); j++) {
          double r=sqrt( pot.pair.sqdist( c.p[sites[i]], c.p[sites[j]] ) );
          du += pot.pair.f * dq*dq/r ;
        }
      u+=du;
      expu+=exp(-du);
      if (ucnt<uvec.size()) {
        a1+=du*exp(-du);
        a2+=uvec[ucnt]*exp(-du);
        a3+=uvec[ucnt];
        //cout << "# GMX U = " << uvec[ucnt] << endl;
        ucnt++;
      }
      if (cnt>100) {
        mu+=-log(expu.avg());
        //expu.reset();
        cnt=0;
      }
    }
    double chempot() { return -log(expu.avg()); }
    double energy() { return ( a1.avg() + a2.avg() ) / expu.avg() - a3.avg(); }
    void info() {
      cout << "# Sites: ";
      for (int i=0; i<sites.size(); i++)
        cout << sites[i] << " ";
      cout << "\n# Excess chemical potential (kT) = " << -log(expu.avg()) << " " << mu.stdev()
           << "\n# <du> <du*exp(-du)> <Uexp(-du)> = " << u.avg() <<" "<< a1.avg() <<" "<< a2.avg()
           << "\n# <U> <exp(-du)>                 = " << a3.avg() <<" " << expu.avg()
           << "\n# Energy (kT)                    = " << energy() << endl;
    }
};

int main() {
  inputfile in("g_widom.conf");
  double dq = in.getflt("dq",0.1);
  XDRFILE *xd;
  int rc,natoms_xtc,step_xtc;
  float time_xtc;
  matrix box_xtc;
  rvec *x_xtc;
  float prec_xtc = 1000.0;
  string xtcname = in.getstr("xtcfile", "traj.xtc");

  box con(100);
  atom.load("faunatoms.dat");
  gmxenergy gmxu("../energy.dat");
  energy<potT> pot(in);

  chargescale cs_one;
  cs_one.uvec=gmxu.u; // copy energy from read energy.dat file (gmx generated)
  chargescale cs_two = cs_one;
  chargescale cs_both = cs_one;
  cs_one.sites.push_back(in.getint("site1", 0));
  cs_two.sites.push_back(in.getint("site2", 1));
  cs_both.sites.push_back(in.getint("site1", 0));
  cs_both.sites.push_back(in.getint("site2", 1));

  // FETCH ATOM NAMES FROM GRO FILE
  iogro gro(in);
  con.p = gro.load( in.getstr("grofile", "conf.gro") );
  atom.reset_properties(con.p); // set all properties according to faunatoms.dat file!

  // OPEN TRAJECTORY FILE
  xd=xdrfile_open(&xtcname[0], "r");
  if (xd!=NULL) {
    cout << "# xtc file opened: " << xtcname << endl;
    rc = read_xtc_natoms(&xtcname[0], &natoms_xtc);        // get length
    if (rc==exdrOK && natoms_xtc==con.p.size()) {
      //x_xtc = (rvec*)calloc(natoms_xtc, sizeof(x_xtc[0])); // allocate memory
      x_xtc = new rvec [natoms_xtc];
      cout << "# Number of atoms = " << natoms_xtc << endl;
      while(1) {                                           // loop over frames
        rc = read_xtc(xd, natoms_xtc, &step_xtc, &time_xtc,box_xtc, x_xtc, &prec_xtc);
        if (rc == 0) {
          con.setvolume( pow(box_xtc[0][0]*10,3) );        // store box size in container (cube assumed)
          pot.pair.setvolume(pow(con.len,3));              // ...and also in potential class.
          cout << "# Frame: " << time_xtc << "  Box = " << con.len << endl;
          for (int i=0; i<con.p.size(); i++) {
            con.p[i].x = x_xtc[i][0]*10.-con.len_half;     // store pos. in container.
            con.p[i].y = x_xtc[i][1]*10.-con.len_half;
            con.p[i].z = x_xtc[i][2]*10.-con.len_half;
          }
          // Analyse frame!!
          //cout << "# Utot = " << pot.system(con.p) << endl;
          cs_one.analyze(con, pot, dq);
          cs_two.analyze(con, pot, dq);
          cs_both.analyze(con, pot, dq);
        }
        else break;
      }
      delete x_xtc;
    }
    xdrfile_close(xd);
  }
  cs_one.info();
  cs_two.info();
  cs_both.info();
  cout << "#\n-----------------------------------------------" << endl
       << "# ddA   = " << cs_both.chempot() - cs_one.chempot() - cs_two.chempot() << endl
       << "# ddE   = " << cs_both.energy() - cs_one.energy() - cs_two.energy() << endl;
};
