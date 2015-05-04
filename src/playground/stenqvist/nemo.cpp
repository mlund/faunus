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










/**
 * @brief Base class for original Wolf based interactions
 *
 * The idea is that this class has no dependencies and is
 * to be used as a helper class for other classes.
 */
class WolfOrgBase {

  struct wdata {
    double r1i_d, r2i, der_dT0c, T1, T1c_r1i, der_dT1c_r1i, T21, T22, T2c2_r2i, der_dT2c1, der_dT2c2_r2i;
  };
  private:
  double rc1, rc1i, rc1i_d, rc2i, kappa, kappa2, constant;
  double dT0c, T1c_rc1, dT1c_rc1, T2c1, T2c2_rc2, dT2c1, dT2c2_rc2;
  wdata data;

  public:
  void calcWolfData(const Point &r) {
    data.r2i = 1/r.squaredNorm();
    double r1i = sqrt(data.r2i);
    data.r1i_d = erfc_x(kappa/r1i)*r1i;
    double der = 0.0;
    double expK = constant*exp(-kappa2/data.r2i);
    data.der_dT0c = der*dT0c;
    data.T1 = (expK + data.r1i_d)*data.r2i;
    data.T1c_r1i = T1c_rc1*r1i;
    data.der_dT1c_r1i = der*dT1c_rc1*r1i;
    data.T21 = -(data.r1i_d + expK)*data.r2i;
    data.T22 = (3.*data.r1i_d*data.r2i + (3.*data.r2i + 2.*kappa2)*expK)*data.r2i;
    data.T2c2_r2i = T2c2_rc2*data.r2i;
    data.der_dT2c1 = der*dT2c1;
    data.der_dT2c2_r2i = der*dT2c2_rc2*data.r2i;
  }

  /**
   * @brief Constructor
   * @param alpha Dampening factor (inverse angstrom)
   * @param rcut Cutoff distance (angstrom)
   */
  WolfOrgBase(double alpha, double rcut) {
    kappa = alpha;
    kappa2 = kappa*kappa;
    constant = 2*kappa/sqrt(pc::pi);
    rc1 = rcut;
    double rc2 = rc1*rc1;
    rc2i = 1/rc2;
    rc1i = 1/rc1;
    double expKc = constant*exp(-kappa2/rc2i);
    rc1i_d = erfc_x(kappa*rc1)*rc1i;

    T1c_rc1 = (expKc + rc1i_d)*rc2i;
    T2c1 = -(expKc + rc1i_d)*rc2i;
    T2c2_rc2 = (3.*rc1i_d*rc2i*rc2i + (3.*rc2i + 2.*kappa2)*expKc*rc2i);

    dT0c = -(expKc + rc1i_d)*rc1i;
    dT1c_rc1 = (-2*T1c_rc1/rc1) - 2*kappa2*expKc*rc1i;
    dT2c1 = -(3*T2c1/rc1) + (2*kappa2*exp(-rc2*kappa2)*rc1i*constant);
    dT2c2_rc2 = (-3*T2c2_rc2/rc1) - (4*kappa2*kappa2*exp(-rc2*kappa2)*rc1i*constant);

    T1c_rc1 = T1c_rc1*rc1;
    dT1c_rc1 = dT1c_rc1*rc1;
    T2c2_rc2 = T2c2_rc2*rc2;
    dT2c2_rc2 = dT2c2_rc2*rc2;
  }

  /**
   * @brief Returns ion-ion interaction.
   * @param qA Charge of ion A
   * @param qB Charge of ion B
   * @param r Direction \f$ r_A - r_B \f$
   */
  template<bool useWdata=false, class Tvec>
    double q2q(double qA, double qB, const Tvec &r) const {
      /* Code to use if calcWolfData is not called */
      if (useWdata==false) {
        double r2i = 1/r.squaredNorm();
        if (r2i < rc2i)
          return 0;
        double r1i = sqrt(r2i);
        double der = 0.0;
        double r1i_d = erfc_x(kappa/r1i)*r1i;
        return (qA*qB*(r1i_d - rc1i_d - der*dT0c));
      }
      return (qA*qB*(data.r1i_d - rc1i_d - data.der_dT0c));
    }

  /**
   * @brief Returns ion-dipole interaction.
   * @param QBxMuA Product of ion B:s charge and dipole A:s scalar
   * @param muA Unit dipole moment vector of particel A
   * @param QAxMuB Product of ion A:s charge and dipole B:s scalar
   * @param muB Unit dipole moment vector of particel B
   * @param r Direction \f$ r_A - r_B \f$
   */
  template<bool useWdata=false, class Tvec>
    double q2mu(double QBxMuA, const Tvec &muA, double QAxMuB, const Tvec &muB, const Tvec &r) const {
      /* Code to use if calcWolfData is not called */
      if (useWdata==false) {
        double r2i = 1/r.squaredNorm();
        if (r2i < rc2i)
          return 0;
        double r1i = sqrt(r2i);
        double T1 = (constant*exp(-kappa2/r2i) + erfc_x(kappa/r1i)*r1i)*r2i;
        double der = 0.0;
        double W1 = QBxMuA*muA.dot(r)*(T1 - T1c_rc1*r1i - der*dT1c_rc1*r1i);
        double W2 = QAxMuB*muB.dot(-r)*(T1 - T1c_rc1*r1i - der*dT1c_rc1*r1i);
        return (W1 + W2);
      }
      double W1 = QBxMuA*muA.dot(r)*(data.T1 - data.T1c_r1i - data.der_dT1c_r1i);
      double W2 = QAxMuB*muB.dot(-r)*(data.T1 - data.T1c_r1i - data.der_dT1c_r1i);
      return (W1 + W2);
    }

  /**
   * @brief Dipole-dipole energy
   * @param muA Dipole moment (A) unit vector
   * @param muB Dipole moment (B) unit vector
   * @param muAxmuB Product of dipole moment scalars, |A|*|B|
   * @param r Direction \f$ r_A - r_B \f$
   */
  template<bool useWdata=false, class Tvec>
    double mu2mu(const Tvec &muA, const Tvec &muB, double muAxmuB, const Tvec &r) const {
      /* Code to use if calcWolfData is not called */
      if (useWdata==false) {
        double r2i = 1.0/r.squaredNorm();
        if (r2i < rc2i)
          return 0;
        double r1i = sqrt(r2i);
        double r1i_d = erfc_x(kappa/r1i)*r1i;
        double expK = constant*exp(-kappa2/r2i);
        double T2_1 = -(r1i_d + expK)*r2i;
        double T2_2 = (3.*r1i_d*r2i + (3.*r2i + 2.*kappa2)*expK)*r2i;
        double der = 0.0;
        double t3 = muA.dot(muB)*(T2_1 - T2c1 - der*dT2c1);
        double t5 = muA.dot(r)*muB.dot(r)*(T2_2 - T2c2_rc2*r2i - der*dT2c2_rc2*r2i);
        return -(t5 + t3)*muAxmuB;
      }
      double t3 = muA.dot(muB)*(data.T21 - T2c1 - data.der_dT2c1);
      double t5 = muA.dot(r)*muB.dot(r)*(data.T22 - data.T2c2_r2i - data.der_dT2c2_r2i);
      return -(t5 + t3)*muAxmuB;
    }

  /**
   * @brief Returns ion-quadrupole energy
   * @param qA Charge of particle A
   * @param quadB Quadrupole moment of particle B
   * @param qB Charge of particle B
   * @param quadA Quadrupole moment of particle A
   * @param r Direction @f$ r_A - r_B @f$
   */
  template<bool useWdata=false, class Tvec, class Tmat>
    double q2quad(double qA, const Tmat &quadB,double qB, const Tmat &quadA, const Tvec &r) const {
      /* Code to use if calcWolfData is not called */
      if (useWdata==false) {
        double r2i = 1/r.squaredNorm();
        if (r2i < rc2i)
          return 0;
        double r1i = sqrt(r2i);
        double r1i_d = erfc_x(kappa/r1i)*r1i;
        double expK = constant*exp(-kappa2/r2i);
        double T2_1 = -(r1i_d + expK)*r2i;
        double T2_2 = (3.*r1i_d*r2i + (3.*r2i + 2.*kappa2)*expK)*r2i;
        double der = 0.0;
        double WAB = r.transpose()*quadB*r;
        WAB = WAB*(T2_2 - T2c2_rc2*r2i - der*dT2c2_rc2*r2i) + quadB.trace()*(T2_1 - T2c1 - der*dT2c1);
        double WBA = r.transpose()*quadA*r;
        WBA = WBA*(T2_2 - T2c2_rc2*r2i - der*dT2c2_rc2*r2i) + quadA.trace()*(T2_1 - T2c1 - der*dT2c1);
        return (qA*WAB + qB*WBA);
      }
      double WAB = r.transpose()*quadB*r;
      WAB = WAB*(data.T22 - data.T2c2_r2i - data.der_dT2c2_r2i) + quadB.trace()*(data.T21 - T2c1 - data.der_dT2c1);
      double WBA = r.transpose()*quadA*r;
      WBA = WBA*(data.T22 - data.T2c2_r2i - data.der_dT2c2_r2i) + quadA.trace()*(data.T21 - T2c1 - data.der_dT2c1);
      return (qA*WAB + qB*WBA);
    }

  /** 
   * @brief Field at `r` due to charge `p` 
   * @param p Particles from which field arises
   * @param r Direction @f$ r_A - r_B @f$
   */
  template<bool useWdata=false, class Tparticle>
    Point fieldCharge(const Tparticle &p, const Point &r) const {
      /* Code to use if calcWolfData is not called */
      if (useWdata==false) {
        double r2i = 1/r.squaredNorm();
        if (r2i < rc2i)
          return Point(0,0,0);
        double r1i = sqrt(r2i);
        double T1 = (constant*exp(-kappa2/r2i) + erfc_x(kappa/r1i)*r1i)*r2i;
        double der = 0.0;
        return (T1 - T1c_rc1*r1i - der*dT1c_rc1*r1i)*r*p.charge;
      }
      return (data.T1 - data.T1c_r1i - data.der_dT1c_r1i)*r*p.charge;
    }

  /** 
   * @brief Field at `r` due to dipole `p` 
   * @param p Particles from which field arises
   * @param r Direction @f$ r_A - r_B @f$
   */
  template<bool useWdata=false, class Tparticle>
    Point fieldDipole(const Tparticle &p, const Point &r) const {
      /* Code to use if calcWolfData is not called */
      if (useWdata==false) {
        double r2i = 1.0/r.squaredNorm();
        if (r2i < rc2i)
          return Point(0,0,0);
        double r1i = sqrt(r2i);
        double r1i_d = erfc_x(kappa/r1i)*r1i;
        double expK = constant*exp(-kappa2/r2i);
        double T2_1 = -(r1i_d + expK)*r2i;
        double T2_2 = (3.*r1i_d*r2i + (3.*r2i + 2.*kappa2)*expK)*r2i;
        double der = 0.0;
        Point t3 = p.mu*(T2_1 - T2c1 - der*dT2c1);
        Point t5 = r*p.mu.dot(r)*(T2_2 - T2c2_rc2*r2i - der*dT2c2_rc2*r2i);
        return (t5 + t3)*p.muscalar;
      }
      Point t3 = p.mu*(data.T21 - T2c1 - data.der_dT2c1);
      Point t5 = r*p.mu.dot(r)*(data.T22 - data.T2c2_r2i - data.der_dT2c2_r2i);
      return (t5 + t3)*p.muscalar;
    }

  double getRc2i() const { return rc2i; }
  double getR2i() const { return data.r2i; }
  double getKappa() { return kappa; }
  double getCutoff() { return rc1; }
};










template<bool useIonIon=false, bool useIonDipole=false, bool useDipoleDipole=false, bool useIonQuadrupole=false>
class MultipoleWolf2 : public PairPotentialBase {
  private:
    WolfBase wolf;
    string _brief() {
      std::ostringstream o;
      o << "Multipole Wolf, lB=" << _lB << textio::_angstrom;
      return o.str();          
    }
  protected:
    double _lB;
  public:
    MultipoleWolf2(InputMap &in, double kappa=0.0, const string &dir="") : wolf(kappa,
        in("cutoff",pc::infty)) {
      name="Multipole Wolf2";
      _lB = pc::lB(in("epsr",1.0));
      _lB = 1.0;

      //cout << "_lB: MultipoleWolf2: " << _lB << endl;
    }
    template<class Tparticle>
      double operator()(const Tparticle &a, const Tparticle &b, const Point &r) {
        double U_total = 0;
        if((useIonIon? 1:0) + (useIonDipole? 1:0) + (useDipoleDipole? 1:0) + (useIonQuadrupole? 1:0) > 1) {
          wolf.calcWolfData(r);
          if (wolf.getR2i() < wolf.getRc2i())
            return 0;
          if(useIonIon == true) U_total += wolf.q2q<true>(a.charge,b.charge,r);
          if(useIonDipole == true) U_total += wolf.q2mu<true>(a.charge*b.muscalar,b.mu,b.charge*a.muscalar,a.mu,r);
          if(useDipoleDipole == true) U_total += wolf.mu2mu<true>(a.mu,b.mu, a.muscalar*b.muscalar, r);
          if(useIonQuadrupole == true) U_total += wolf.q2quad<true>(a.charge, b.theta,b.charge, a.theta,r);
          return _lB*U_total;
        }
        if(useIonIon == true) U_total += wolf.q2q(a.charge,b.charge,r);
        if(useIonDipole == true) U_total += wolf.q2mu(a.charge*b.muscalar,b.mu,b.charge*a.muscalar,a.mu,r);
        if(useDipoleDipole == true) U_total += wolf.mu2mu(a.mu,b.mu, a.muscalar*b.muscalar, r);
        if(useIonQuadrupole == true) U_total += wolf.q2quad(a.charge, b.theta,b.charge, a.theta,r);
        return _lB*U_total;
      }

    template<bool useIon=true, bool useDipole=true, class Tparticle>
      Point field(const Tparticle &p, const Point &r) {
        if(useIon && useDipole) {
          wolf.calcWolfData(r);
          if (wolf.getR2i() < wolf.getRc2i())
            return Point(0,0,0);
          Point E = wolf.fieldCharge<true>(p,r);
          E += wolf.fieldDipole<true>(p,r);
          return _lB*E;
        }
        if(useIon == true) return _lB*wolf.fieldCharge(p,r);
        if(useDipole == true) return _lB*wolf.fieldDipole(p,r);
        return Point(0,0,0);
      }

    /**
     * @brief Interaction of dipole `p` with field `E`, see 'Intermolecular and SUrface Forces' by J. Israelachvili, p. 97 eq. 5.15
     * @todo Needs to be tested!
     */
    template<class Tparticle>
      double fieldEnergy(const Tparticle &p, const Point &E) {
        return 0;
      }

    string info(char w) {
      using namespace textio;
      std::ostringstream o;
      o << pad(SUB,w,"Temperature") << pc::T() << " K" << endl
        << pad(SUB,w,"Bjerrum length") << _lB << " "+angstrom << endl
        << pad(SUB,w,"Cutoff") << wolf.getCutoff() << " "+angstrom << endl
        << pad(SUB,w,"Kappa") << wolf.getKappa() << " "+angstrom+"^-1" << endl;
      return o.str();
    }
};








template<bool useIonIon=false, bool useIonDipole=false, bool useDipoleDipole=false, bool useIonQuadrupole=false>
class MultipoleOrgWolf : public PairPotentialBase {
  private:
    WolfOrgBase wolf;
    string _brief() {
      std::ostringstream o;
      o << "Multipole Org Wolf, lB=" << _lB << textio::_angstrom;
      return o.str();          
    }
  protected:
    double _lB;
  public:
    //template<bool useIonIon=false, bool useIonDipole=false, bool useDipoleDipole=false, bool useIonQuadrupole=false>
    MultipoleOrgWolf(InputMap &in, double kappa=0.0, const string &dir="") : wolf(kappa,
        in("cutoff",pc::infty)) {
      name="Multipole Org Wolf";
      _lB = pc::lB(in("epsr",1.0));
      _lB = 1.0;

      //cout << "_lB: MultipoleOrgWolf: " << _lB << endl;
    }
    template<class Tparticle>
      double operator()(const Tparticle &a, const Tparticle &b, const Point &r) {
        double U_total = 0;
        if((useIonIon? 1:0) + (useIonDipole? 1:0) + (useDipoleDipole? 1:0) + (useIonQuadrupole? 1:0) > 1) {
          wolf.calcWolfData(r);
          if (wolf.getR2i() < wolf.getRc2i())
            return 0;
          if(useIonIon == true) U_total += wolf.q2q<true>(a.charge,b.charge,r);
          if(useIonDipole == true) U_total += wolf.q2mu<true>(a.charge*b.muscalar,b.mu,b.charge*a.muscalar,a.mu,r);
          if(useDipoleDipole == true) U_total += wolf.mu2mu<true>(a.mu,b.mu, a.muscalar*b.muscalar, r);
          if(useIonQuadrupole == true) U_total += wolf.q2quad<true>(a.charge, b.theta,b.charge, a.theta,r);
          return _lB*U_total;
        }
        if(useIonIon == true) U_total += wolf.q2q(a.charge,b.charge,r);
        if(useIonDipole == true) U_total += wolf.q2mu(a.charge*b.muscalar,b.mu,b.charge*a.muscalar,a.mu,r);
        if(useDipoleDipole == true) U_total += wolf.mu2mu(a.mu,b.mu, a.muscalar*b.muscalar, r);
        if(useIonQuadrupole == true) U_total += wolf.q2quad(a.charge, b.theta,b.charge, a.theta,r);
        return _lB*U_total;
      }

    template<bool useIon=true, bool useDipole=true, class Tparticle>
      Point field(const Tparticle &p, const Point &r) {
        if(useIon && useDipole) {
          wolf.calcWolfData(r);
          if (wolf.getR2i() < wolf.getRc2i())
            return Point(0,0,0);
          Point E = wolf.fieldCharge<true>(p,r);
          E += wolf.fieldDipole<true>(p,r);
          return _lB*E;
        }
        if(useIon == true) return _lB*wolf.fieldCharge(p,r);
        if(useDipole == true) return _lB*wolf.fieldDipole(p,r);
        return Point(0,0,0);
      }

    /**
     * @brief Interaction of dipole `p` with field `E`, see 'Intermolecular and SUrface Forces' by J. Israelachvili, p. 97 eq. 5.15
     * @todo Needs to be tested!
     */
    template<class Tparticle>
      double fieldEnergy(const Tparticle &p, const Point &E) {
        return 0;
      }

    string info(char w) {
      using namespace textio;
      std::ostringstream o;
      o << pad(SUB,w,"Temperature") << pc::T() << " K" << endl
        << pad(SUB,w,"Bjerrum length") << _lB << " "+angstrom << endl
        << pad(SUB,w,"Cutoff") << wolf.getCutoff() << " "+angstrom << endl
        << pad(SUB,w,"Kappa") << wolf.getKappa() << " "+angstrom+"^-1" << endl;
      return o.str();
    }
};

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
  InputMap in("madelung.json");  

  int case0 = 1;    // 0 = SodiumChloride, 1 = Dipole, 2 = CesiumChloride
  double kappa = 0.8;
  int writePos = 0;

/*
  IonIonQ pairQ(in);
  IonIonFanourgakis pairFN(in);
  Coulomb pairT(in);
  MultipoleOrgWolf<true,false,false,false> pairW(in,0.0);
  MultipoleWolf2<true,false,false,false> pairF(in,0.0);
  MultipoleOrgWolf<true,false,false,false> pairWD(in,kappa);
  MultipoleWolf2<true,false,false,false> pairFD(in,kappa);
*/

     DipoleDipoleQ pairQ(in);
     DipoleDipoleFanourgakis pairFN(in);
     DipoleDipole pairT(in);
     MultipoleOrgWolf<false,false,true,false> pairW(in,0.0);
     MultipoleWolf2<false,false,true,false> pairF(in,0.0);
     MultipoleOrgWolf<false,false,true,false> pairWD(in,kappa);
     MultipoleWolf2<false,false,true,false> pairFD(in,kappa);
    


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
