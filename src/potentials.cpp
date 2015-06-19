#ifndef FAU_POT_BASE_H
#define FAU_POT_BASE_H

#include <faunus/common.h>
#include <faunus/point.h>
#include <faunus/inputfile.h>
#include <faunus/physconst.h>
#include <faunus/geometry.h>
#include <faunus/potentials.h>
#include <faunus/textio.h>

namespace Faunus {

  namespace Potential {

    PairPotentialBase::PairPotentialBase() {
      if ( atom.empty() )
        std::cerr << "Warning: no atoms defined when initializing pair potential.\n";
      rcut2.resize(atom.size());
    }

    PairPotentialBase::PairPotentialBase( const std::string &sec ) : jsonsec(sec) {
      assert( ! jsonsec.empty() );
      if ( atom.empty() )
        std::cerr << "Warning: no atoms defined when initializing pair potential.\n";
      rcut2.resize(atom.size());
    }

    PairPotentialBase::~PairPotentialBase() {
    }

    /*
     * @param N Maximum number of atom types
     * @param rc Default cutoff distance (angstrom)
     *
     void PairPotentialBase::initCutoff(size_t N, float rcut) {
     rcut2.setConstant(N,N,rcut*rcut);
     }*/

    /*
     * @param i Particle type i
     * @param j Particle type j
     * @param rc Cutoff distance (angstrom)
     void PairPotentialBase::setCutoff(size_t i, size_t j, float rcut) {
     rcut2(i,j)=rcut2(j,i)=rcut*rcut;
     }*/

    std::string PairPotentialBase::_brief() {
      assert(!name.empty() && "Provide a descriptive name for the potential");
      return name;
    }

    std::string PairPotentialBase::info(char w) {
      return name+": N/A";
    }

    string PairPotentialBase::brief() {
      assert(!name.empty() && "Potential must have a name.");
      return _brief();
    }

    void PairPotentialBase::test(UnitTest&) {}

    string Harmonic::_brief() {
      using namespace Faunus::textio;
      std::ostringstream o;
      o << name + ": k=" << k << kT+"/"+angstrom+squared+" req=" << req << _angstrom; 
      return o.str();
    }

    Hertz::Hertz( Tmjson &j, const string &sec) : PairPotentialBase(sec) {
      name = "Hertz";
      E = j[sec]["_E"] = 0.0;
    }

    string Hertz::_brief() {
      using namespace Faunus::textio;
      std::ostringstream o;
      o <<name+" E:"<< E;
      return o.str();
    }

    string Hertz::info(char w) {
      using namespace Faunus::textio;
      std::ostringstream o;
      o <<name+" E: "<< E;
      return o.str();
    }

    YukawaGel::YukawaGel(Tmjson &j, const string &sec) : Coulomb(j, sec) {
      name = "YukawaGel";

      Z  = j[sec]["yukawagel_Z"]  | 0.0;
      nc = j[sec]["yukawagel_nc"] | 0.0;
      ns = j[sec]["yukawagel_ns"] | 0.0;
      v  = j[sec]["yukawagel_v"]  | 1.0;
      d  = j[sec]["yukawagel_d"]  | 1.0;

      k = sqrt(4.*pc::pi*(nc+2.*ns)*v*v*bjerrumLength());
      Z2e2overER = Z*Z*bjerrumLength();
      kd = k*d;
      k2d2 = k*k*d*d;
      ekd = exp(-kd);
      braket7 = (cosh(kd/2.) - (2.*sinh(kd/2.)/kd));
      cout << Z <<"   "<< bjerrumLength() << "  " << k << endl ;
    }

    string YukawaGel::_brief() {
      using namespace Faunus::textio;
      std::ostringstream o;
      o <<name;
      return o.str();
    }

    string YukawaGel::info(char w) {
      using namespace Faunus::textio;
      std::ostringstream o;
      o <<name+" Z: "<< Z << endl << " Counter ions:  "<< nc <<  endl
        <<" Salt: " << ns << endl <<" BjerrumL: " << bjerrumLength() << endl
        <<" kappa: " << k <<endl;

      return o.str();
    }

    CosAttract::CosAttract( Tmjson &j, const string &sec) : PairPotentialBase(sec) {
      name="CosAttract";
      eps = j[sec]["eps"] | 0.0;
      rc  = j[sec]["rc" ] | 0.0;
      wc  = j[sec]["wc" ] | 0.0;
      rc2=rc*rc;
      c=pc::pi/2/wc;
      rcwc2=pow((rc+wc),2);
    }

    string CosAttract::_brief() {
      using namespace Faunus::textio;
      std::ostringstream o;
      o << name+": "+epsilon+"=" << eps << kT+" rc=" << rc << _angstrom
        << " wc=" << wc << _angstrom;
      return o.str();
    }

    string CosAttract::info(char w) {
      using namespace textio;
      std::ostringstream o;
      o << pad(SUB,w,"Depth") << eps << kT << endl
        << pad(SUB,w,"Decay length") << wc << _angstrom << endl
        << pad(SUB,w,"Width") << rc << _angstrom << endl
        << pad(SUB,w,"More info") << "doi:10/chqzjk" << endl;
      return o.str();
    }

    HardSphere::HardSphere() {
      name="Hardsphere";
    }

    string HardSphere::info(char w) {
      using namespace Faunus::textio;
      return textio::indent(SUB)+name+"\n";
    }

    string LennardJones::_brief() {
      std::ostringstream o;
      o << name << ": " << textio::epsilon+"(LJ)=" << eps/4 << textio::kT;
      return o.str();
    }

    string LennardJones::info(char w) {
      using namespace Faunus::textio;
      std::ostringstream o;
      o << pad(SUB,w+1,epsilon+"(LJ)") << eps/4 << kT
        << " = " << eps/4.0_kJmol << " kJ/mol" << endl;
      return o.str();
    }

    LennardJonesTrunkShift::LennardJonesTrunkShift(Tmjson &j, const string &sec) : LennardJones( j,sec) {
      name+=" Truncated and shifted to sigma";
    }

    /**
     * @param in is scanned for the keywords `prefix_threshold` (angstrom)
     *        and `prefix_depth` (kT).
     * @param dir Root section for input
     */
    SquareWell::SquareWell( Tmjson &j, const string &sec) : PairPotentialBase(sec) {
      name="Square Well";
      threshold = j[sec]["threshold"] | 0.0;
      depth     = j[sec]["depth"] | 0.0;
    }

    string SquareWell::_brief() {
      std::ostringstream o;
      o << name << ": u=" << depth << textio::kT + " r=" << threshold;
      return o.str();
    }

    string SquareWell::info(char w) {
      using namespace Faunus::textio;
      std::ostringstream o;
      o << pad(SUB,w,"Threshold") << threshold << " "+angstrom+" (surface-surface)" << endl;
      o << pad(SUB,w,"Depth") << depth << kT << endl;
      return o.str();
    }

    /**
     * In addition to the keywords from `Potential::SquareWell` the json object is
     * searched for:
     * - `threshold_lower`
     */
    SquareWellShifted::SquareWellShifted(Tmjson &j, const string &sec): SquareWell( j, sec ) {
      name+=" Shifted";
      threshold_lower = j[sec]["threshold_lower"] | 0.0;
    }

    string SquareWellShifted::_brief() {
      std::ostringstream o;
      o << name << ": u=" << depth << textio::kT
        << " range=" << threshold_lower << "-" << threshold;
      return o.str();
    }

    string SquareWellShifted::info(char w) {
      using namespace Faunus::textio;
      std::ostringstream o;
      o << SquareWell::info(w)
        << pad(SUB,w,"Threshold_lower") << threshold_lower << _angstrom
        << " (surface-surface)" << endl;
      return o.str();
    }

    string R12Repulsion::_brief() {
      std::ostringstream o;
      o << name << ": " << textio::epsilon+"(r12)=" << eps/4 << textio::kT;
      return o.str();
    }

    string R12Repulsion::info(char w) {
      using namespace Faunus::textio;
      std::ostringstream o;
      o << pad(SUB,w+1,epsilon+"(r12_rep)") << eps/4 << kT << endl;
      return o.str();
    }

    /**
     * The following keywords are searched searched in section `coulomb`:
     * - `temperature` [Kelvin, default = 298.15]
     * - `epsilon_r` - relative dielectric constant. Default is 80.
     * - `depsdt` - temperature dependence of dielectric constant,
     *   \f$ \partial\epsilon_r/\partial T\approx-0.368\f$ for water.
     */
    Coulomb::Coulomb(Tmjson &j, const string &sec) : PairPotentialBase( sec ) {
      assert( ! j[sec].is_null() );
      name      = "Coulomb";
      epsilon_r = j[sec]["epsr"] | 80.0;
      depsdt    = ( j[sec]["depsdt"] | -0.368 ) * pc::T() / epsilon_r;
      lB        = pc::lB( epsilon_r );
    }

    string Coulomb::_brief() {
      std::ostringstream o;
      o << name << ": lB=" << lB << " eps_r=" << epsilon_r << " T=" << pc::T();
      return o.str();
    }

    double Coulomb::bjerrumLength() const {
      return lB;
    }

    string Coulomb::info(char w) {
      using namespace textio;
      std::ostringstream o;
      o << pad(SUB,w,"Temperature") << pc::T() << " K" << endl
        << pad(SUB,w,"Dielectric constant") << epsilon_r << endl
        << pad(SUB,w+6,"T"+partial+epsilon+"/"+epsilon+partial+"T") << depsdt << endl
        << pad(SUB,w,"Bjerrum length") << lB << " "+angstrom << endl;

      return o.str();
    }

    void Coulomb::test(UnitTest &t) {
      t("bjerrum", bjerrumLength(), 1e-6);
    }

    CoulombWolf::CoulombWolf( Tmjson &j, const string &sec ) : Coulomb( j, sec ) {
      double Rc = j[sec]["cutoff"] | 10.0;
      Rcinv = 1/Rc;
      Rc2 = Rc*Rc;
      name += "Wolf/Yonezawa";
    }

    string CoulombWolf::info(char w) {
      using namespace textio;
      std::ostringstream o;
      o << Coulomb::info(w)
        << pad(SUB,w,"More info") << "doi:10/j97\n"
        << pad(SUB,w,"Cut-off") << 1/Rcinv << _angstrom+"\n";
      return o.str();
    }

    ChargeNonpolar::ChargeNonpolar(Tmjson &j, const string &sec) : Coulomb(j, sec) {
      name="Charge-Nonpolar";
      c = 0.5 * bjerrumLength() * ( j[sec]["excess_polarization"] | -1.0 );
    }

    string ChargeNonpolar::info(char w) {
      std::ostringstream o;
      o << Coulomb::info(w)
        << textio::pad(textio::SUB,w,"Excess polarization") << 2*c*bjerrumLength() << endl;
      return o.str();
    }

    string DebyeHuckel::_brief() {
      std::ostringstream o;
      o << Coulomb::_brief() << " I=" << ionicStrength();
      return o.str();
    }

    double DebyeHuckel::ionicStrength() const {
      return k*k/c;
    }

    double DebyeHuckel::debyeLength() const {
      return 1/k;
    }

    /**
     * @details The Debye-Huckel potential is temperature dependent and contains entropy
     * contributions from both solvent and salt degrees of freedom.
     * This function return the entropy of interaction between a pair of
     * particles interacting with an effective Debye-Huckel potential. This is done by
     * taking the temperature derivate of w(R):
     *
     * @f[
     * S(r_{ij})/k_B = -\frac{ \partial w(r_{ij},T) } {k_B \partial T}
     *     = \beta w_{ij}\left [ \alpha - \frac{\kappa r_{ij}(\alpha+1)}{2}\right ]
     * @f]
     * where \f$ \alpha=T \partial \epsilon_r/\epsilon_r\partial T\f$
     * is determined experimentally for pure water. To get the entropy from salt ions
     * only, set \f$\alpha=0\f$ via the following keywords:
     *
     * @param  betaw    Inter particle free energy, \f$\beta w\f$, in units of kT.
     * @param  r        Inter particle distance
     * @return Interaction entropy \f$ S(r_{ij})/k_B = \beta TS(r_{ij})\f$
     * @todo   Optimize
     */
    double DebyeHuckel::entropy(double betaw, double r) const {
      return betaw * (depsdt - 0.5*k*r*(depsdt+1));
    }

    /**
     * @return \f$\beta \mu_{\mbox{excess}} = -\frac{l_Bz^2\kappa}{2(1+\kappa a)}\f$
     * @param z Charge number
     * @param a Particle diameter (angstrom)
     */
    double DebyeHuckel::excessChemPot(double z, double a) const {
      return -lB*z*z*k / ( 2 * (1+k*a) );
    }

    /**
     * @return \f$\exp {(\beta\mu_{\mbox{excess}})}\f$
     * @param z Charge number
     * @param a Particle diameter (angstrom)
     */
    double DebyeHuckel::activityCoeff(double z, double a) const {
      return exp( excessChemPot(z,a) );
    }

    string DebyeHuckel::info(char w) {
      using namespace Faunus::textio;
      std::ostringstream o;
      o << Coulomb::info(w)
        << pad(SUB,w,"Ionic strength") << ionicStrength() << " mol/l" << endl
        << pad(SUB,w+1,"Debye length, 1/"+textio::kappa) << debyeLength() << " "
        << 1/k << " "+angstrom << endl;
      if ( k2_count_avg.cnt > 0 ) {
        double k2_s = k*k - k2_count;
        o << pad( SUB, w+1, "Debye length, 1/"+textio::kappa )
          << 1/sqrt( k2_count_avg.avg()) << " "+angstrom+" (counter ions)\n"
          << pad( SUB, w+1,"Debye length, 1/"+textio::kappa )
          << 1/sqrt( k2_s ) << " "+angstrom+" (salt)" << endl;
      }
      return o.str();
    }

  } //Potential namespace

} //Faunus namespace
#endif
