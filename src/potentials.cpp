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

    PairPotentialBase::PairPotentialBase() {}

    PairPotentialBase::~PairPotentialBase() { }

    /**
     * @param N Maximum number of atom types
     * @param rc Default cutoff distance (angstrom)
     */
    void PairPotentialBase::initCutoff(size_t N, float rc) {
      cutoff2.setConstant(N,N,rc*rc);
    }

    /**
     * @param i Particle type i
     * @param j Particle type j
     * @param rc Distance cutoff (angstrom) 
     */
    void PairPotentialBase::setCutoff(size_t i, size_t j, float rc) {
      cutoff2(i,j)=cutoff2(j,i)=rc*rc;
    }

    std::string PairPotentialBase::info(char w) {
      return name+": N/A";
    }
    
    /**
     * @param a First particle
     * @param b Second particle
     * @param r2 Squared distance between them (angstrom squared)
     */
    double PairPotentialBase::operator() (const particle &a, const particle &b, double r2) const {
      assert(!"Pair energy not defined!");
      return pc::infty;
    }

    double PairPotentialBase::operator() (const particle &a, const particle &b, const Point &r2) const {
      return operator()(a,b,r2.squaredNorm());
    }

    /**
     * @param a First particle
     * @param b Second particle
     * @param r2 Squared distance between them (angstrom squared)
     * @param p Vector from: p=b-a
     */
    Point PairPotentialBase::force(const particle &a, const particle &b, double r2, const Point &p) {
      assert(!"Force not overrided!");
      return Point(0,0,0);
    }

    /**
     * @param a Particle emanating the field
     * @param r Position in which to calculate the field
     */
    Point PairPotentialBase::field(const particle &a, const Point &r) const {
      return Point(0,0,0);
    }

    /**
     * This will reset the temperature to the specified value. By default this function
     * does nothing, although in Debug mode it will throw an exception if derived classes
     * do not implement it (and is called).
     */
    void PairPotentialBase::setTemperature(double) {
      assert(!"Not implemented.");
    }

    string PairPotentialBase::brief() {
      assert(!name.empty() && "Potential must have a name.");
      return _brief();
    }

    bool PairPotentialBase::save(string filename, particle::Tid ida, particle::Tid idb) {
      std::ofstream f(filename.c_str());
      if (f) {
        double min=0.9 * (atom[ida].radius+atom[idb].radius);
        particle a,b;
        a = atom[ida];
        b = atom[idb];
        f << "# Pair potential: " << brief() << endl
          << "# Atoms: " << atom[ida].name << "<->" << atom[idb].name << endl;
        for (double r=min; r<=150; r+=0.5)
          f << std::left << std::setw(10) << r << " " << operator()(a,b,r*r) << endl; 
        return true;
      }
      return false;
    }

    void PairPotentialBase::test(UnitTest&) {}

    Harmonic::Harmonic(double forceconst, double eqdist) : k(forceconst), req(eqdist) {
      name="Harmonic";
    }

    Harmonic::Harmonic(InputMap &in, string pfx) {
      name="Harmonic";
      k  = in.get<double>( pfx+"forceconst", 0);
      req = in.get<double>( pfx+"eqdist", 0);
    }

    string Harmonic::_brief() {
      using namespace Faunus::textio;
      std::ostringstream o;
      o << name << ": k=" << k << kT+"/"+angstrom+squared+" req=" << req << _angstrom; 
      return o.str();
    }

    /**
     * @param k_kT   Stiffness or bond strength [kT]
     * @param rmax_A Maximum length after which the energy goes to infinity [angstrom]
     */
    FENE::FENE(double k_kT, double rmax_A) : k(k_kT) {
      name="FENE";
      r02=rmax_A*rmax_A;
      r02inv=1/r02;
    }

    FENE::FENE(InputMap &in, string pfx) {
      name="FENE";
      k  = in.get<double>( pfx+"stiffness", 0);
      r02 = pow( in.get<double>( pfx+"maxsep", 0), 2);
      r02inv = 1/r02;
    }

    string FENE::_brief() {
      using namespace Faunus::textio;
      std::ostringstream o;
      o << name+": k=" << k << kT+"/"+angstrom+squared+" r0=" << sqrt(r02) << _angstrom; 
      return o.str();
    }

    CosAttract::CosAttract(InputMap &in, string pfx) {
      name="CosAttract";
      eps = in.get<double>( pfx+"eps", 0);
      rc  = in.get<double>( pfx+"rc", 0);
      wc  = in.get<double>( pfx+"wc", 0);
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

    string HardSphere::_brief() {
      return name;
    }

    string HardSphere::info(char w) {
      using namespace Faunus::textio;
      return textio::indent(SUB)+name+"\n";
    }

    LennardJones::LennardJones() : eps(0) {
      name="Lennard-Jones";
    }

    /**
     * @param in InputMap is scanned for the `lj_eps` and should be in units of kT
     * @param pfx Prefix for InputMap - default is `ls_`
     */
    LennardJones::LennardJones(InputMap &in, string pfx) {
      name="Lennard-Jones";
      eps = 4*in.get<double>( pfx+"eps", 0);
      string unit = in.get<string>(pfx+"unit", "kT");
      if (unit=="kJ/mol")
        eps=eps/pc::kT2kJ(1.);
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
        << " = " << pc::kT2kJ(eps/4) << " kJ/mol" << endl;
      return o.str();
    }

    WeeksChandlerAndersen::WeeksChandlerAndersen(InputMap &in) :
      Tbase(in), onefourth(1/4.), twototwosixth(std::pow(2,2/6.))  {
        name="WeeksChandlerAnderson";
      }

    LorentzBerthelot::LorentzBerthelot() : name("Lorentz-Berthelot Mixing Rule") {}

    double LorentzBerthelot::mixSigma(double sigma1, double sigma2) const {
      return 0.5*(sigma1+sigma2);
    }

    double LorentzBerthelot::mixEpsilon(double eps1, double eps2) const {
      return sqrt(eps1*eps2);
    }

    LennardJonesR12::LennardJonesR12(InputMap &in, string pfx) : LennardJones(in,pfx) {
      name+="R12";
    }

    LennardJonesTrunkShift::LennardJonesTrunkShift(InputMap &in, string pfx) : LennardJones(in,pfx) {
      name+=" Truncated and shifted to sigma";
    }

    /**
     * @param in is scanned for the keywords `prefix_threshold` (angstrom)
     *        and `prefix_depth` (kT).
     * @param prefix InputMap keyword prefix. Default is `squarewell`
     */
    SquareWell::SquareWell(InputMap &in, string prefix) {
      name="Square Well";
      threshold = in.get<double>(prefix+"_threshold", 0, name+" upper threshold (AA)");
      depth     = in.get<double>(prefix+"_depth", 0, name+" depth (kT)");
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
     * In addition to the keywords from `Potential::SquareWell` the InputMap is
     * searched for:
     * - `prefix_threshold_lower`
     */
    SquareWellShifted::SquareWellShifted(InputMap &in, string prefix): SquareWell(in,prefix) {
      name=name + " Shifted";
      threshold_lower = in.get<double>(prefix+"_threshold_lower", 0, name+" lower threshold (AA)");
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

    SquareWellHydrophobic::SquareWellHydrophobic(InputMap &in, string prefix) : SquareWell(in,prefix) {
      name="Hydrophobic " + name;
    }

    /**
     * @param in InputMap is scanned for the keyword `softrep_sigma` which should be in angstrom
     */
    SoftRepulsion::SoftRepulsion(InputMap &in) {
      name="Repulsive r6";
      sigma6 = pow( in.get<double>( "softrep_sigma", 5 ), 6);
    }

    string SoftRepulsion::_brief() {
      std::ostringstream o;
      o << name << ": " << textio::sigma  << pow(sigma6,1/6.) << textio::_angstrom;
      return o.str();
    }

    string SoftRepulsion::info(char w) {
      using namespace Faunus::textio;
      std::ostringstream o;
      o << textio::pad(SUB,w+1,textio::sigma) << pow(sigma6,1/6.) << _angstrom << endl;
      return o.str();
    }

    R12Repulsion::R12Repulsion() {
      name="r12-Repulsion";
    }

    /**
     * @param in InputMap is scanned for the keyword `lj_eps` and should be in units of kT
     * @param pfx InputMap prefix
     */
    R12Repulsion::R12Repulsion(InputMap &in, string pfx) {
      name="r12-Repulsion";
      eps = 4*in.get<double>( pfx+"eps", 0.05, name+" epsilon (kT)" );
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
     * The following input keywords are searched searched:
     * - `temperature` [Kelvin, default = 298.15]
     * - `epsilon_r` - relative dielectric constant. Default is 80.
     * - `depsdt` - temperature dependence of dielectric constant,
     *   \f$ \partial\epsilon_r/\partial T\approx-0.368\f$ for water.
     */
    Coulomb::Coulomb(InputMap &in) {
      name="Coulomb";
      pc::setT ( in.get<double>("temperature", 298.15, "Absolute temperature (K)") );
      epsilon_r = in.get<double>("epsilon_r",80., "Dielectric constant");
      depsdt = in.get<double>("depsdt", -0.368, "See documentation") * pc::T() / epsilon_r;
      lB=pc::lB( epsilon_r );
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

    CoulombWolf::CoulombWolf(InputMap &in) : Coulomb(in) {
      double Rc=in.get<double>("coulomb_cut", 10.);
      Rcinv=1/Rc;
      Rc2=Rc*Rc;
      name+="Wolf/Yonezawa";
    }

    string CoulombWolf::info(char w) {
      using namespace textio;
      std::ostringstream o;
      o << Coulomb::info(w)
        << pad(SUB,w,"More info") << "doi:10.1063/1.4729748\n"
        << pad(SUB,w,"Cut-off") << 1/Rcinv << _angstrom+"\n";
      return o.str();
    }

    ChargeNonpolar::ChargeNonpolar(InputMap &in) : Coulomb(in) {
      name="Charge-Nonpolar";
      c=bjerrumLength()/2*in.get<double>("excess_polarization", -1);
    }

    string ChargeNonpolar::info(char w) {
      std::ostringstream o;
      o << Coulomb::info(w)
        << textio::pad(textio::SUB,w,"Excess polarization") << 2*c*bjerrumLength() << endl;
      return o.str();
    }

    /**
     * In addition to the keywords from Potential::Coulomb, InputMap is searched for:
     *
     * - `dh_ionicstrength` [mol/l] 
     * - `dh_debyelength` [angstrom] (only if I=0, default)
     */
    DebyeHuckel::DebyeHuckel(InputMap &in) : Coulomb(in) {
      double I;
      const double zero=1e-10;
      name="Debye-Huckel";
      c=8 * lB * pc::pi * pc::Nav / 1e27;
      I=in.get<double>("dh_ionicstrength",0, "Ionic strength (mol/l)");  // [mol/l]
      k=sqrt( I*c );
      if (k<zero)
        k=1/in.get<double>("dh_debyelength", 1/zero, "Debye length (AA)"); // [A]
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
     * only, set \f$\alpha=0\f$ via the InputMap.
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
      o << Coulomb::info(w);
      o << pad(SUB,w,"Ionic strength") << ionicStrength() << " mol/l" << endl;
      o << pad(SUB,w+1,"Debye length, 1/"+textio::kappa) << debyeLength() << " "+angstrom << endl;
      return o.str();
    }

    DebyeHuckelShift::DebyeHuckelShift(InputMap &in) : DebyeHuckel(in) {
      double rc=in.get<double>("pairpot_cutoff",pc::infty);
      sqcutoff=rc*rc;
      shift = exp(-rc*k)/rc;
      std::ostringstream o;
      o << " (shifted, rcut=" << rc << textio::_angstrom << ")";
      name+=o.str();
    }

    /**
     * \f$ \beta u(r) = l_B \frac{ z_1 z_2 }{r}\f$
     */
    double MultipoleEnergy::ionion(double z1, double z2, double r) {
      return lB*z1*z2/r;
    }

    /**
     * \f$ \beta u(r) = -l_B \frac{ z a_z }{r^2}\f$
     */
    double MultipoleEnergy::iondip(double z, const Point &a, double r) {
      return -lB*z*a.z()/(r*r);
    }

    /**
     * \f$ \beta u(r) = l_B \frac{a_x b_x + a_y b_y - 2a_z b_z  }{r^3}\f$
     */
    double MultipoleEnergy::dipdip(const Point &a, const Point &b, double r) {
      return lB*( a.x()*b.x() + a.y()*b.y() - 2*a.z()*b.z() ) / (r*r*r);
    }

    /* 
       double PatchSCsphere::eattractive_sc_sphere(const Point &a, const Point &b, const Point r_cm) {
       double atrenergy, a, b, f0, halfl;
       struct vector vec1;

    //TODO if we dont have calculate distance segment to point - distvec a dist a contt point

    //calculate closest distance attractive energy
    if (dist < interact->param->pdis)
    atrenergy = -interact->param->epsilon;
    else {
    atrenergy = cos(PIH*(interact->dist-interact->param->pdis)/interact->param->pswitch);
    atrenergy *= -atrenergy*interact->param->epsilon ;
    }
    //scaling function: angular dependence of patch1
    if (b.halfl < 1e-6) {
    which = 0;
    vec1=vec_perpproject(distvec, a.dir);
    vec1.normalize();
    a = vec1.dot(a.patchdir);
    halfl=a.halfl;
    } else {
    which = 1;
    vec1=vec_perpproject(distvec, b.dir);
    vec1.normalize();
    a = vec1.dot(b.patchdir);
    halfl=b.halfl;
    }
    //scaling function for the length of spherocylinder within cutoff

    b = sqrt(rcut*rcut-dist*dist);
    if ( contt + b > halfl ) 
    f0 = halfl;
    else 
    f0 = contt + b;
    if ( contt - b < -halfl ) 
    f0 -= -halfl;
    else 
    f0 -= interact->contt - b;
    atrenergy *= fanglscale(a,interact->param, which)*(f0+1.0);

    return atrenergy;
    }//TODO atrenergy from ndist, cutoff rcut at the beginning, epsilon atd
    */

  } //Potential namespace

} //Faunus namespace
#endif
