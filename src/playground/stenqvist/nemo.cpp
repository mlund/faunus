#define DIPOLEPARTICLE
#include <faunus/faunus.h>
#include <faunus/multipole.h>
#include <functional>
#include <iostream>
using namespace Faunus;                     
using namespace Faunus::Move;
using namespace Faunus::Potential;

typedef Space<Geometry::Cuboid,DipoleParticle> Tspace; 
typedef CoulombWolf TpairDDW;
typedef LennardJonesLB TpairLJ;

typedef CombinedPairPotential<TpairLJ,TpairDDW> Tpair;

class DipoleDipoleTest : public DipoleDipole {
  private:
    string _brief() { return "DipoleDipole Test"; }
    double rc1, rc1i, rc2;
    int P;
  public:
    DipoleDipoleTest(InputMap &in, const string &dir="") : DipoleDipole(in) { 
      name += " Test"; 
      rc1  = in( "cutoff", pc::infty );
      rc1i = 1.0/rc1;
      rc2 = rc1*rc1;
      P  = in( "moments", 2 );
      _lB = 1.0;
    }

    template<class Tparticle>
      double operator()(const Tparticle &a, const Tparticle &b, const Point &r) const {
        double r2 = r.squaredNorm();
        if(r2 < rc2) {
          double q = sqrt(r2)/rc1;

          double Q = 1.0;
          double qk = q*q*q;
          for(int i = 1; i < P; i++) {
            Q = Q*(1 - qk);
            qk = qk*q;
          }

          return DipoleDipole::operator()(a,b,r)*Q;
        }
        return 0.0;
      }

    string info(char w) {
      using namespace textio;
      std::ostringstream o;
      o << DipoleDipole::info(w)
        << pad(SUB,w,"Cutoff") << rc1 << " "+angstrom << endl;
      return o.str();
    }
};

class IonIonTest : public Coulomb {
  private:
    string _brief() { return "Coulomb Test"; }
    double rc1, rc1i, rc2, _lB;
    int P;
  public:
    IonIonTest(InputMap &in, const string &dir="") : Coulomb(in) { 
      name += " Test"; 
      _lB = Coulomb(in,dir).bjerrumLength();
      rc1  = in( "cutoff", pc::infty );
      rc1i = 1.0/rc1;
      rc2 = rc1*rc1;
      P  = in( "moments", 2 );

      _lB = 1.0;
      //cout << "_lB: IonIonTest: " << _lB << endl;
    }

    template<class Tparticle>
      double operator()(const Tparticle &a, const Tparticle &b, double r2) const {
        double r1 = sqrt(r2);
        if(r2 < rc2) {
          double q = r1/rc1;

          double Q = 1.0;
          double qk = q;
          for(int i = 1; i < P; i++) {
            Q = Q*(1 - qk);
            qk = qk*q;
          }
          return _lB*(a.charge*b.charge/r1)*Q;
        }
        return 0.0;
      }

    template<class Tparticle>
      double operator()(const Tparticle &a, const Tparticle &b, const Point &r) const {
        double r2 = r.squaredNorm();
        return operator()(a,b,r2);
      }

    string info(char w) {
      using namespace textio;
      std::ostringstream o;
      o << Coulomb::info(w)
        << pad(SUB,w,"Cutoff") << rc1 << " "+angstrom << endl;
      return o.str();
    }
};






template<class Tpairpot, class Tid>
bool savePotential(Tpairpot pot, Tid ida, Tid idb, string file) {
  std::ofstream f(file.c_str());
  if (f) {
    DipoleParticle a,b;
    a=atom[ida];
    b=atom[idb];
    a.mu = Point(1,0,0);
    b.mu = Point(1,0,0);
    for (double r=0.5; r<=4.5; r+=0.05) {
      f << std::left << std::setw(10) << r << " "
        << pot(a,b,Point(r,0,0)) << endl;
    }
    return true;
  }
  return false;
}

bool checkRc(double RcC) {
  double rest = RcC;

  while(rest > 1.0) {
    rest = rest - 1.0;
  }

  if(rest > 0.4) {
    if(rest < 0.6) {
      return false;
    }
  }
  return true;
}

void generateLattice(Tspace &spc, double RcC, double side, double charge, int &cnt, int &closeCenter, int case0) {
  double RcC2 = RcC*RcC;

  if(case0 == 0) {
    for(double x = -RcC; x < RcC+side; x += side) {
      for(double y = -RcC; y < RcC+side; y += side) {
        for(double z = -RcC; z < RcC+side; z += side) {
          Point vec = Point(x,y,z);
          double r2 = vec.squaredNorm();
          if(r2 < RcC2) {
            if(r2 < 1e-6)
              closeCenter = cnt;
            spc.p[cnt] = Point(x,y,z);
            spc.p[cnt].mu = Point(charge,0,0);
            spc.p[cnt].muscalar = 1.0;
            spc.p[cnt++].charge = charge;
          }
          charge = -charge;
        }
      }
    }
    for(int n = cnt; n < spc.p.size(); n++) {
      spc.p[n] = Point(RcC,RcC,RcC);
      spc.p[n].mu = Point(1,0,0);
      spc.p[n].muscalar = 0.0;
      spc.p[n].charge = 0.0;
    }
    spc.trial = spc.p;
  } else if(case0 == 1) {
    for(double x = -RcC; x < RcC+side; x += side) {
      for(double y = -RcC; y < RcC+side; y += side) {
        double z = 0.0;
        Point vec = Point(x,y,z);
        double r2 = vec.squaredNorm();
        if(r2 < RcC2) {
          if(r2 < 1e-6)
            closeCenter = cnt;
          spc.p[cnt] = Point(x,y,z);
          spc.p[cnt].mu = Point(charge,0,0);
          spc.p[cnt].muscalar = 1.0;
          spc.p[cnt++].charge = charge;
        }
        charge = -charge;
      }
    }
    for(int n = cnt; n < spc.p.size(); n++) {
      spc.p[n] = Point(RcC,RcC,RcC);
      spc.p[n].mu = Point(1,0,0);
      spc.p[n].muscalar = 0.0;
      spc.p[n].charge = 0.0;
    }
    spc.trial = spc.p;
  } else {
    
    if(checkRc(RcC)) {
      charge = -charge;
    } else {

    }

    double extra = 0.0;
    for(double x = -RcC; x < RcC+side; x += 0.5*side) {
      charge = -charge;
      if(charge < 0) {
	  extra = side/2.0;
      } else {
	extra = 0.0;
      }

      for(double y = -RcC; y < RcC+side; y += side) {
        for(double z = -RcC; z < RcC+side; z += side) {
          Point vec = Point(x,y+extra,z+extra);
          double r2 = vec.squaredNorm();
          if(r2 < RcC2) {
            if(r2 < 1e-6)
              closeCenter = cnt;
            spc.p[cnt] = Point(x,y+extra,z+extra);
            spc.p[cnt].mu = Point(charge,0,0);
            spc.p[cnt].muscalar = 1.0;
            spc.p[cnt++].charge = charge;
          }
        }
      }
    }
    for(int n = cnt; n < spc.p.size(); n++) {
      spc.p[n] = Point(RcC,RcC,RcC);
      spc.p[n].mu = Point(1,0,0);
      spc.p[n].muscalar = 0.0;
      spc.p[n].charge = 0.0;
    }
    spc.trial = spc.p;

    
    
    
    /*
    if(checkRc(RcC)) {
      charge = -charge;
    } else {

    }

    double extra = 0.0;
    int cntX = 0;
    for(double x = -RcC; x < RcC+side; x += side) {
      if(cntX > 0) {
        cntX = -1;
        if(checkRc(RcC)) {
          extra = side/2.0;
        } else {
          extra = 0.0;
        }
      } else {
        if(checkRc(RcC)) {
          extra = 0.0;
        } else {
          extra = side/2.0;
        }
      }
      charge = -charge;

      for(double y = -RcC; y < RcC+side; y += side) {
        for(double z = -RcC; z < RcC+side; z += side) {
          Point vec = Point(x,y+extra,z+extra);
          double r2 = vec.squaredNorm();
          if(r2 < RcC2) {
            if(r2 < 1e-6)
              closeCenter = cnt;
            spc.p[cnt] = Point(x,y+extra,z+extra);
            spc.p[cnt].mu = Point(charge,0,0);
            spc.p[cnt].muscalar = 1.0;
            spc.p[cnt++].charge = charge;
          }
        }
      }
      cntX++;
    }
    for(int n = cnt; n < spc.p.size(); n++) {
      spc.p[n] = Point(RcC,RcC,RcC);
      spc.p[n].mu = Point(1,0,0);
      spc.p[n].muscalar = 0.0;
      spc.p[n].charge = 0.0;
    }
    spc.trial = spc.p;
    */


  }
}


int main() {
  InputMap in0("madelung0.json");  
  InputMap in("madelung.json");  

  int case0 = 1;    // 0 = SodiumChloride, 1 = Dipole, 2 = CesiumChloride
  double kappa = 0.8;
  int writePos = 0;

/*
  IonIonQ pairQ(in);
  IonIonFanourgakis pairFN(in);
  Coulomb pairT(in);
  MultipoleWolf<true,false,false,false> pairW(in0);
  MultipoleWolf<true,false,false,false> pairF(in0);
  MultipoleWolf<true,false,false,false> pairWD(in);
  MultipoleWolf<true,false,false,false> pairFD(in);
*/

     DipoleDipoleQ pairQ(in);
     DipoleDipoleFanourgakis pairFN(in);
     DipoleDipole pairT(in);
     MultipoleWolf<false,false,true,false,false> pairW(in0);
     MultipoleWolf<false,false,true,false,true> pairF(in0);
     MultipoleWolf<false,false,true,false,false> pairWD(in);
     MultipoleWolf<false,false,true,false,true> pairFD(in);
    


  Tspace spc(in);          
  Group sol;
  sol.addParticles(spc, in);

  int N = spc.p.size();
  int cnt = 0;
  int charge = 1;

  double side = 0.5;
  if(case0 == 2) {
      side = 2.0*side;
  }
  in.cd ("system/coulomb");
  double Rc = in("cutoff_A",pc::infty);

  in.cd();

  int steps = int(Rc/side) + 2;
  double RcC = double(steps)*side;
  int closeCenter = -1;
  
  
  RcC = 4.0;

  generateLattice(spc, RcC, side, charge, cnt, closeCenter, case0);

  if(writePos == 1) {
    for(int n = 0; n < cnt; n++)
      cout << spc.p[n].transpose() << " " << spc.p[n].charge << endl;
    return 0;
  }


  savePotential(Coulomb(in), atom["sol"].id, atom["sol"].id, "pot_dipdip_C.dat");
  savePotential(pairQ, atom["sol"].id, atom["sol"].id, "pot_dipdip_Q.dat");
  savePotential(pairFN, atom["sol"].id, atom["sol"].id, "pot_dipdip_FA.dat");
  savePotential(pairT, atom["sol"].id, atom["sol"].id, "pot_dipdip_T.dat");
  savePotential(pairW, atom["sol"].id, atom["sol"].id, "pot_dipdip_W.dat");
  savePotential(pairF, atom["sol"].id, atom["sol"].id, "pot_dipdip_F.dat");
  savePotential(pairWD, atom["sol"].id, atom["sol"].id, "pot_dipdip_WD.dat");
  savePotential(pairFD, atom["sol"].id, atom["sol"].id, "pot_dipdip_FD.dat");
  
  return 0;

  double energyQ = 0.0;
  double energyFN = 0.0;
  double energyWolf = 0.0;
  double energyFennel = 0.0;
  double energyWolfD = 0.0;
  double energyFennelD = 0.0;
  double energyTest = 0.0;
  Point r(0,0,0);



  for(int n = 0; n < N; n++) {
    r = spc.p[n];
    if(r.norm() < 1e-6)
      continue;
    energyQ += pairQ.operator()(spc.p[n],spc.p[closeCenter],r);
    energyFN += pairFN.operator()(spc.p[n],spc.p[closeCenter],r);
    energyWolf += pairW.operator()(spc.p[n],spc.p[closeCenter],r);
    energyFennel += pairF.operator()(spc.p[n],spc.p[closeCenter],r);
    energyWolfD += pairWD.operator()(spc.p[n],spc.p[closeCenter],r);
    energyFennelD += pairFD.operator()(spc.p[n],spc.p[closeCenter],r);
    energyTest += pairT.operator()(spc.p[n],spc.p[closeCenter],r);
  }
  //cout << "Five! " << closeCenter << ", RcC: " << RcC << ", side: " << side << ", bool: " << checkRc(RcC) << endl;

  double energyMadelung;

  if(case0 == 0) {
    energyMadelung = -(3.496/(2*side))*spc.p[closeCenter].charge*spc.p[closeCenter].charge;
  } else if(case0 == 1) {
    energyMadelung = 10.58354613*spc.p[closeCenter].charge*spc.p[closeCenter].charge;
  } else {
    energyMadelung = -(1.76267)*spc.p[closeCenter].charge*spc.p[closeCenter].charge;
  }

 

  
  Analysis::MultipoleAnalysis dian(spc,in);
  dian.setCutoff(in("cutoff",pc::infty));
  dian.sampleDP(spc);
  
 cout << in("cutoff",pc::infty) << " " << energyMadelung << " " << energyQ << " " << energyFN << " " << energyWolf << " " << energyWolfD << " " << energyFennel << " " << energyFennelD << " " << energyTest << endl;


  return 0;
}
