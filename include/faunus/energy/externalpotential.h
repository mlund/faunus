#ifndef FAUNUS_EXTERNALPOTENTIAL_H
#define FAUNUS_EXTERNALPOTENTIAL_H

#include "faunus/energy/base.h"
#include "faunus/energy/springinteraction.h"

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
   * take inte account explicit interactions within the container.
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
      //!< (sign typo in manuscript: phi^infty(z) should be "-2*pi*z" on page 413, middle)
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
            << "#     More information:        Mol. Phys. 1996, 87:407\n";
        return o.str();
      }

      void save(cuboid &c) {
        std::ofstream f("akesson.dat");
        for (double z=-c.len_half.z; z<=c.len_half.z; z+=dz)
          if (rho(z).cnt>0)
            f << z << " " << rho(z).avg() << " " << phi(z) << endl;
      }

      double energy( const particle &a ) { return a.charge*getPotential(a); }

      double energy( const vector<particle> &p ) {
        double u=0;
        for (int i=0; i<p.size(); i++)
          if (p[i].charge!=0)
            u += p[i].charge * getPotential(p[i]);
        return u; // in kT
      }

      double energy( const vector<particle> &p, const group &g ) {
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
          double a=c.len_half.x, uold=energy(c.p);
          for (double z=-c.len_half.z; z<=c.len_half.z; z+=dz) {
            double s=0;
            for (double zn=-c.len_half.z; zn<=c.len_half.z; zn+=dz)
              if (rho(zn).cnt>0)
                s += rho(zn).avg() * phi_ext( std::abs(z-zn), a );  // Eq. 14 in Greberg paper
              else s+=0;
            phi(z) = lB*s;
          }
          return energy(c.p) - uold;
        }
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
          return springinteraction<Tpairpot>::energy(p,a) + expot.energy(a);
        }

        double u_monomer(vector<particle> &p, const polymer &g, unsigned int i ) {
          return springinteraction<Tpairpot>::u_monomer(p,g,i) + expot.energy( p[i] );
        }

        double uself_polymer(vector<particle> &p, const polymer &g) {
          return springinteraction<Tpairpot>::uself_polymer(p,g) + expot.energy(p,g);
        }

        double internal(vector<particle> &p, const group &g, int step=1) { 
          return springinteraction<Tpairpot>::internal(p,g,step) + expot.energy(p,g);
        }

        double energy(vector<particle> &p) {
          return springinteraction<Tpairpot>::energy(p) + expot.energy(p);
        }

        double energy(vector<particle> &p, int i) {
          return springinteraction<Tpairpot>::energy(p,i) + expot.energy(p[i]);
        }

        double energy(vector<particle> &p, const group &g) {
          return springinteraction<Tpairpot>::energy(p,g) + expot.energy(p,g);
        }
    };

}//namespace 
#endif
