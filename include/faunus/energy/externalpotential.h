#ifndef FAUNUS_EXTERNALPOTENTIAL_H
#define FAUNUS_EXTERNALPOTENTIAL_H

#include "faunus/energy/base.h"
#include "faunus/energy/springinteraction.h"
#include "faunus/energy/penaltyfunction.h"
#include <faunus/physconst.h>

namespace Faunus {
  /*!
   * \brief Mean field electric potential from outside rectangular simulation box.
   * \author Mikael Lund
   * \date Asljunga, December 2010.
   * \note Currently works with the cuboid simulation container
   *
   * This class will calculate the average potential outside a simulation box due to ion
   * densities inside the box, extended to infinity.
   * The update() function calculates the average charge densities in slits in the xy-plane
   * and with a spacing dz. This is used to evaluate the electric potential along
   * the z-axis by using the method by Torbjorn and CO (Mol. Phys. 1996, 87:407). To avoid
   * energy drifts, update() returns the energy change brought about by updating the charge profile.
   * This class should be used in conjunction with an interaction class (energybase derivative) to
   * take into account explicit interactions within the container.
   */

  class expot_akesson {
    private:
      unsigned int phiupdate;                 //!< Only update phi every phiupdate'th time
      unsigned int cnt;                       //!< Number of charge density updates
      bool enabled;                           //!< Set to true to enable potential (default)
      double dz;                              //!< z spacing between slits (A)
      double lB;                              //!< Bjerrum length (A)
      xytable< double, average<double> > rho; //!< Charge density at z (unit A^-2)
      xytable< double, double > phi;          //!< External potential at z (unit: beta*e)

      double getPotential(const point &a) { return phi(a.z); }

      //!< This is Eq. 15 of the mol. phys. 1996 paper by Greberg et al.
      //!< (sign typo in paper: phi^infty(z) should be "-2*pi*z" on page 413, middle)
      double phi_ext(double z, double a) {
        return  -2.0*pyc.pi*z-8.0*a*log((sqrt(2.0*a*a+z*z)+a)/sqrt(a*a+z*z))+2.0*z*(pyc.pi/
            2.0+std::asin((a*a*a*a-z*z*z*z-2.0*a*a*z*z)/pow(a*a+z*z,2.0)));
      }

    public:
      expot_akesson(inputfile &in) {
        cnt=0;
        phiupdate=10;
        dz = in.getflt("akesson_dz", 0.1);
        double zmin = -in.getflt("cuboid_zlen",0)/2 - dz;
        lB = in.getflt("bjerrum",0);
        rho.init(dz, zmin, -zmin);
        phi.init(dz, zmin, -zmin);
        enabled = in.getboo("akesson_enabled", true);
        if (enabled==false)
          lB=0;
      }

      string info() {
        std::ostringstream o;
        if (enabled==true || lB>1e-6)
          o << "#   Akesson external potential:\n"
            << "#     Bjerrum length         = " << lB << " A\n"
            << "#     xy-slit resolution     = " << dz << " A (" << phi.y.size() << " slits)\n"
            << "#     Number of pot. updates = " << cnt/phiupdate << endl
            << "#     More information:        doi:10.1080/00268979600100271\n";
        return o.str();
      }

      void save(cuboid &c) {
        std::ofstream f("akesson.dat");
        for (double z=-c.len_half.z; z<=c.len_half.z; z+=dz)
          if (rho(z).cnt>0)
            f << z << " " << rho(z).avg() << " " << phi(z) << endl;
      }

      void dump(string filename) {
        phi.dumptodisk(filename);
      }

      double energy_particle( const particle &a ) { return a.charge*getPotential(a); }

      double energy_vector( const vector<particle> &p ) {
        double u=0;
        for (int i=0; i<p.size(); i++)
          if (p[i].charge!=0)
            u += p[i].charge * getPotential(p[i]);
        return u; // in kT
      }

      double energy_group( const vector<particle> &p, const group &g ) {
        double u=0;
        for (int i=g.beg; i<=g.end; i++)
          if (p[i].charge!=0)
            u += p[i].charge * getPotential(p[i]);
        return u; // in kT
      }

      /*!
       * \brief Updated xy-slit charge densities as well as the potential along z.
       * \param c Cuboid container where charged particles are sought out and averaged
       * \return Energy change of updating the external potential (use to avoid energy drifts)
       */
      double update(cuboid &c) {
        cnt++;
        double area=c.len.x*c.len.y;
        for (double z=-c.len_half.z; z<=c.len_half.z; z+=dz) { // update rho(z)
          double Q=0;
          for (int i=0; i<c.p.size(); i++)
            if (c.p[i].z>=z)
              if (c.p[i].z<z+dz)
                Q+=c.p[i].charge;
          rho(z) += Q / area; 
        }

        if (cnt % phiupdate == 0) {  // update phi(z) - but not so often
          double a=c.len_half.x, uold=energy_vector(c.p);
          for (double z=-c.len_half.z; z<=c.len_half.z; z+=dz) {
            double s=0;
            for (double zn=-c.len_half.z; zn<=c.len_half.z; zn+=dz)
              if (rho(zn).cnt>0)
                s += rho(zn).avg() * phi_ext( std::abs(z-zn), a );  // Eq. 14 in Greberg paper
              else s+=0;
            phi(z) = lB*s;
          }
          return energy_vector(c.p) - uold;
        }
        return 0;
      }
  };

  /*!
   * \brief Tabulated external potential
   * \author Chris Evers
   * \date Lund, February 2011
   * \note Currently works with the cuboid simulation container
   *
   * This class is used to introduce a potential in a simulation box which
   * depends on the z-position of each particle and is stored in a table.
   * 
   * This class should be used as a basis for a specific external potential and 
   * in conjunction with an interaction class (energybase derivative) to
   * take into account explicit interactions within the container.
   */
  class expot_table {
  private:

  public:
    bool enabled;                           //!< Set to true to enable potential (default)
    xytable< double, double > phi;          //!< External potential at z [kT]
    double dz;                              //!< Table resolution
    double len_halfz;                       //!< Half cuboid length in z-direction
    unsigned int cnt;                       //!< Number of potential updates
    string name;                            //!< Name of external potential

    expot_table(inputfile &in) {
      name = "No external potential";
      enabled = in.getboo("expot_enabled", true);
      cnt=0;
      dz = in.getflt("expot_dz", 0.1);
      len_halfz = in.getflt("cuboid_zlen",0)/2;
      if (len_halfz<1e-6)
        len_halfz = in.getflt("cuboid_len",0)/2;
      double zmin = -len_halfz - dz;
      phi.init(dz, zmin, -zmin);
    }

    string info() {
      std::ostringstream o;
      if (enabled==true)
        o << "#   " << name << ":" << endl
          << "#     Potential resolution   = " << dz   << " A (" << phi.y.size() << " slits)" << endl
          << "#     Number of pot. updates = " << cnt  << endl;
        return o.str();
    }

    void dump(string filename) {
      phi.dumptodisk(filename);
    }

    double getPotential(const point &a) { 
      return phi.interpolate(a.z); 
    }

    virtual double energy_particle( const particle &a ) { 
      return getPotential(a); // in kT
    }

    double energy_vector( const vector<particle> &p ) {
      double u=0;
      for (int i=0; i<p.size(); i++)
        u += energy_particle(p[i]);
      return u; // in kT
    }

    double energy_group( const vector<particle> &p, const group &g ) {
      double u=0;
      for (int i=g.beg; i<=g.end; i++)
        u += energy_particle(p[i]);
      return u; // in kT
    }

    double update(cuboid &c) {
      cnt++;
      if (enabled==true) {
        double uold=energy_vector(c.p);
        for (double z=-len_halfz; z<=len_halfz; z+=dz)
          phi(z)=calcPotential(len_halfz-z);
        return energy_vector(c.p) - uold;
      }
      return 0;
    }

    virtual double calcPotential(double z) { 
      return 0;
    }
  };

  /*!
   * \brief Gouy-Chapman interaction with a planar charged surface
   * \author Chris Evers
   * \date Lund, February 2011
   * \note Currently works with the cuboid simulation container
   *
   * This class will calculate the potential in a simulation box due to a charged planar surface
   * and screend by an electrical double layer with a corresponding Debye screening length.
   * 
   * This class is based on expot_table and should be used in conjunction with an interaction class 
   * (energybase derivative) to take into account explicit interactions within the container.
   */

  class expot_gouychapman : public expot_table {
  private:
    double lB;                              //!< Bjerrum length (A)
    double k;                               //!< Inverse debye length (A-1)
    double I;                               //!< Ionic strenght (M)
    double c0;                              //!< Ion concentration (A-3)
    double rho;                             //!< Surface charge density (e A-2)
    double phi0;                            //!< Unitless surface potential \frac{\phi0 e}{kT}
    double gamma0;                          //!< Gouy-chapman coefficient ()

  public:
    expot_gouychapman(inputfile &in) : expot_table(in) {
      name = "Gouy-Chapman external potential";

      lB = in.getflt("bjerrum",0);
      if (enabled==false)
        lB=0;

      k=1/in.getflt("debyelen",1.1e6);              // Inverse debye length
      if ( 1/k<1e6)
        I=k*k*1e27/(8*pyc.pi*lB*pyc.Nav);
      else {
        I=in.getflt("ionicstr",0);
        k=sqrt( 8*pyc.pi*lB*pyc.Nav/1e27*I );       // k=\sqrt{8 \pi \lambda_B[A] I[A-3]}
      }
      c0=I*pyc.Nav/1e27;                            // assuming 1:1 salt, so c0=I

      phi0=in.getflt("expot_phi0",1.1e6);           // Surface potential [V]
      if ( fabs(phi0)<1e6 ) {
        phi0=phi0*pyc.e/(pyc.kB*pyc.T);             // Unitless surface potential \frac{\phi_0 e}{kT}
        rho=sqrt(2*c0/(pyc.pi*lB))*sinh(.5*phi0);   // \rho = \sqrt\frac{2 c_0}{\pi l_B}  \sinh(\frac{\phi0 e}{2kT}) [Evans & Wennerström, 1999, Colloidal Domain p 138-140]
      }
      else {
        rho=1/in.getflt("expot_qarea",0);
        if (rho>1e20)
          rho=in.getflt("expot_rho",0);
        phi0=2.*asinh(rho * sqrt(.5*lB*pyc.pi/c0 ));// \frac{\phi0 e}{kT}=2arcsinh(\rho \sqrt\frac{{\pi l_B}{2 c_0}}) [Evans..]
      }

      gamma0=tanh(phi0/4);                          // assuming z=1 \Gamma_0=\tanh{\frac{z\phi_0 e}{4kT}} [Evans..]
    }

    string info() {
      std::ostringstream o;
      if (enabled==true || lB>1e-6)
        o << expot_table::info()
          << "#     Bjerrum length         = " << lB   << " A "  << endl
          << "#     Debye length           = " << 1./k << " A "  << endl
          << "#     Ionic strenght         = " << I*1e3<< " mM " << endl
          << "#     Bulk ion concentration = " << c0   << " A-3 " << endl
          << "#     Surface potential      = " << phi0*pyc.kB*pyc.T/pyc.e << " V  " << endl
          << "#     Unitless surface pot   = " << phi0 << "   " << endl
          << "#     Area per surface charge= " << 1/rho<< " A2 " << endl
          << "#     Surface charge density = " << rho*pyc.e*1e20  << " C m-2 " << endl
          << "#     GC-coefficient Gamma_0 = " << gamma0  << "  " << endl;
        return o.str();
    }

    double energy_particle( const particle &a ) { 
      return a.charge*getPotential(a); // in kT
    }

    double calcPotential(double z) { 
      double phiz;
      double exponent=exp(-k*z);        //\exp{-\kappa z}
      phiz=2 * log((1+gamma0*exponent)/(1-gamma0*exponent));  // \frac{Phi z e}{kT}=2\ln{\frac{1+\Gamma_0 \exp{-\kappa z}}{1-\Gamma_0 \exp{-\kappa z}}}
      return phiz;
    }
  };
  
  /*!
   * \brief Hydrophobic interaction with a planar surface
   * \author Chris Evers
   * \date Lund, February 2011
   * \note Currently works with the cuboid simulation container
   *
   * This class will calculate the potential in a simulation box due to a hydrophoic planar surface
   * 
   * This class is based on expot_table and should be used in conjunction with an interaction class 
   * (energybase derivative) to take into account explicit interactions within the container.
   */
  class expot_hydrophobic : public expot_table {
  private:
    double u;                               //!< Hydrophobic energy (kT)
    double r;                               //!< Hydrophic interaction distance (A)
  public:
    expot_hydrophobic(inputfile &in) : expot_table(in) {
      name = "Hydrophobic external potential";
      u = in.getflt("wallphob_u",0);
      if (enabled==false)
        u=0;
      r = in.getflt("wallphob_r",0);
    }

    string info() {
      std::ostringstream o;
      if (enabled==true)
        o << expot_table::info()
          << "#     Interaction strength   = " << u    << " kT " << endl
          << "#     Interaction length     = " << r    << " A "  << endl;
        return o.str();
    }

    double energy_particle( const particle &a ) { 
      if (a.hydrophobic==true)
        return getPotential(a); // in kT
      else
        return 0;
    }

    double calcPotential(double z) { 
      if (z <= r)
        return u;
      else 
        return 0;
    }
  };

  /*!
   * \brief Spring interaction class with arbitrary external potential
   * \author Mikael Lund
   * \date Asljunga, December 2010.
   *
   * In addition to the normal springinteraction class, the particles interact with
   * an arbitrary external potential given as a template parameter (Texpot). Texpot is
   * expected to have functions that describe how the system, particles and groups interact
   * with the external potential.
   */
  template<class Tpairpot, class Texpot>
    class springinteraction_expot : public springinteraction<Tpairpot> {
      public:
        using springinteraction<Tpairpot>::energy;
        Texpot expot;
        springinteraction_expot(inputfile &in) : springinteraction<Tpairpot>(in), expot(in) {
          interaction<Tpairpot>::name+=" + external potential";
        }

        string info() {
          std::ostringstream o;
          o << springinteraction<Tpairpot>::info() << expot.info();
          return o.str();
        }

        double energy(vector<particle> &p, const particle &a) {
          return springinteraction<Tpairpot>::energy(p,a) + expot.energy_particle(a);
        }

        double u_monomer(vector<particle> &p, const polymer &g, unsigned int i ) {
          return springinteraction<Tpairpot>::u_monomer(p,g,i) + expot.energy_particle( p[i] );
        }

        double uself_polymer(vector<particle> &p, const polymer &g) {
          return springinteraction<Tpairpot>::uself_polymer(p,g);// + expot.energy(p,g);
        }

        double internal(vector<particle> &p, const group &g, int step=1) { 
          return springinteraction<Tpairpot>::internal(p,g,step) + expot.energy_group(p,g);
        }

        double energy(vector<particle> &p) {
          return springinteraction<Tpairpot>::energy(p) + expot.energy_vector(p);
        }

        double energy(vector<particle> &p, int i) {
          return springinteraction<Tpairpot>::energy(p,i) + expot.energy_particle(p[i]);
        }

        double energy(vector<particle> &p, const group &g) {
          return springinteraction<Tpairpot>::energy(p,g) + expot.energy_group(p,g);
        }
    };

  /*!
   * \brief Spring interaction class with arbitrary external potential and penaltyfunction
   * \author Chris Evers
   * \date Lund, March 2011.
   *
   * In addition to the normal springinteraction class, the particles interact with
   * an arbitrary external potential and a penaltyfunction. The external potential is
   * given as a template parameter (Texpot). Texpot is expected to have functions that 
   * describe how the system, particles and groups interact with the external potential.
   */
  template<class Tpairpot, class Texpot>
  class springinteraction_expot_penalty : public springinteraction_expot<Tpairpot,Texpot> {
    private:
      container* conPtr;
      group* groupPtr;
      double zlen_half;
    public:
      using springinteraction_expot<Tpairpot,Texpot>::energy;
      penaltyfunction pen;

      springinteraction_expot_penalty(inputfile &in, container &c, polymer &p)
        : springinteraction_expot<Tpairpot,Texpot>(in),
        pen( 0, in.getflt("cuboid_zlen",100), 0.25, in.getflt("penalty_energy",.1), in.getflt("penalty_scalingfactor",1)) {
            zlen_half=.5*in.getflt("cuboid_zlen",100);
            conPtr=&c;
            groupPtr=&p;
            interaction<Tpairpot>::name+=" + penalty function";
      }

      string info() {
        std::ostringstream o;
        o << springinteraction_expot<Tpairpot,Texpot>::info() << pen.info();
        return o.str();
      }

      // Energy in momomer and crankshaft functions
      double u_monomer(vector<particle> &p, const polymer &g, unsigned int i ) {
        double u=0;
        if (&g==groupPtr) {
          point cm = g.masscenter(*conPtr, p);
          u=pen.energy(zlen_half-cm.z);
        }
        return springinteraction_expot<Tpairpot,Texpot>::u_monomer(p,g,i) + u;
      }

      // Energy in translation function
      double energy(vector<particle> &p, const group &g) {  
        double u=0;
        if (&g==groupPtr) {
          point cm = g.masscenter(*conPtr, p);
          u=pen.energy(zlen_half-cm.z);
        }
        return springinteraction_expot<Tpairpot,Texpot>::energy(p,g) + u;
      }
  };
}//namespace 
#endif
