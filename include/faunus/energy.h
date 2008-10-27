#ifndef FAU_ENERGY_H
#define FAU_ENERGY_H

#include "faunus/common.h"
#include "faunus/group.h"
#include "faunus/inputfile.h"
#include "faunus/hardsphere.h"
#include "faunus/potentials/pot_coulomb.h"
#include "faunus/potentials/pot_minimage.h"
#include "faunus/potentials/pot_hsminimage.h"
#include "faunus/potentials/pot_hydrophobic.h"
#include "faunus/potentials/pot_datapmf.h"
#include "faunus/potentials/pot_debyehuckel.h"
#include "faunus/potentials/pot_debyehuckelP3.h"
#include "faunus/potentials/pot_hccoulomb.h"
#include "faunus/potentials/pot_barecoulomb.h"
#include "faunus/potentials/pot_netz.h"

namespace Faunus {
  /*!
   *  \brief Base class for interactions between particles and groups
   *  \author Mikael Lund
   *
   *  Calculates interaction energies between particles and groups. The
   *  pair potential is specified as a template type which allows inlining.
   *  Unless otherwise specified, all energies will be returned in units of \b kT.
   */
  class energybase {
    protected:
      string name;
    public:
      double tokT;
      energybase(double f) { tokT=f; }
      virtual double energy(const particle &, const particle &)=0;                       //!< particle<->particle (slow!)
      virtual double energy(const vector<particle> &, const particle &)=0;               //!< all<->external particle
      virtual double energy(const vector<particle> &, int)=0;                            //!< all<->particle i.
      virtual double energy(const vector<particle> &, const group &)=0;                  //!< all<->group.
      virtual double energy(const vector<particle> &)=0;                                 //!< all<->all (System energy).
      virtual double energy(const vector<particle> &, const group &, const group &)=0;   //!< group<->group.
      virtual double energy(const vector<particle> &, const group &, int)=0;             //!< group<->particle i.
      virtual double energy(const vector<particle> &, const group &, const particle &)=0;//!< group<->external particle.
      virtual double energy(const vector<particle> &, const vector<macromolecule> &)=0;  //!< vector<group> <-> vector<group>
      virtual double potential(const vector<particle> &, unsigned short)=0;              //!< Electric potential at j'th particle
      virtual double internal(const vector<particle> &, const group &)=0;                //!< internal energy in group
      virtual double pot(const vector<particle> &, const point &)=0;                     //!< Electrostatic potential in a point
      virtual double dipdip(const point &, const point &, double)=0;                     //!< Dipole-dipole energy.
      virtual double iondip(const point &, double, double)=0;                            //!< Ion-dipole energy.

      string info() {
        std::ostringstream o;
        o << endl
          << "# ENERGY EVALUATION:" << endl
          << "#   Scheme:             " << name << endl;
        return o.str();
      }
  };

  /*!
   * \brief Implementation of all energy functions
   * \author Mikael Lund
   *
   * This is an expression template that constructs the
   * energy functions of a particular pair potential. All
   * returned energies are in units of kT.
   *
   * \code
   * inputfile in("input.conf");
   * interaction<pot_coulomb> pot(in);
   * pot.energy(...);
   * \endcode
   */
  template<class T> class interaction : public energybase {
    public:
      T pair; //!< An instance of the pair-potential.
      interaction(inputfile const &in) : pair(in), energybase(pair.f) {
        name="Standard";
        tokT=pair.f;
      };

      double energy(const particle &a, const particle &b) {
        return pair.pairpot(a,b) * pair.f;
      }

      double energy(const vector<particle> &p, int j) {
        int i,ps=p.size();
        double u=0;
        //#pragma omp parallel for reduction (+:u)
        for (i=0; i<j; i++)
          u+=pair.pairpot( p[i],p[j] );
        //#pragma omp parallel for reduction (+:u)
        for (i=j+1; i<ps; i++)
          u+=pair.pairpot( p[i],p[j] );
        return pair.f*u;
      }

      double energy(const vector<particle> &p, const group &g) {
        int n=g.end+1, psize=p.size();
        double u=0;
#pragma omp parallel for reduction (+:u)
        for (int i=g.beg; i<n; i++) {
          for (int j=0; j<g.beg; j++)
            u += pair.pairpot(p[i],p[j]);
          for (int j=n; j<psize; j++)
            u += pair.pairpot(p[i],p[j]);
        }
        return pair.f*u;
      }

      double energy(const vector<particle> &p, const group &g, int j) {
        double u=0;
        int len=g.end+1;
        if (g.find(j)==true) {   //avoid self-interaction...
          for (int i=g.beg; i<j; i++)
            u+=pair.pairpot(p[i],p[j]);
          for (int i=j+1; i<len; i++)
            u+=pair.pairpot(p[i],p[j]);
        } else                        //simple - j not in g
          for (int i=g.beg; i<len; i++)
            u+=pair.pairpot(p[i],p[j]);
        return pair.f*u;  
      }

      double energy(const vector<particle> &p, const group &g, const particle &a) {
        if (g.beg==-1) return 0;
        double u=0;
        int i,n=g.end+1;
        //#pragma omp parallel for reduction (+:u)
        for (i=g.beg; i<n; i++)
          u+=pair.pairpot(a, p[i]); 
        return pair.f*u;
      }

      double energy(const vector<particle> &p, const vector<macromolecule> &g) {
        double u=0;
        int k,j=g.size(),t=p.size();
        for (int l=0; l<j; l++) {
          k=g[l].end+1;
          for (int i=g[l].beg; i<k; i++) {
            for (int s=(g[l].end+1); s<t; s++) {
              u+=pair.pairpot(p[i],p[s]);
            }
          }
        }
        return pair.f*u;
      }

      double energy(const vector<particle> &p) {
        double u=0;
        int n = p.size();
#pragma omp parallel for reduction (+:u)
        for (int i=0; i<n-1; ++i)
          for (int j=i+1; j<n; ++j)
            u+=pair.pairpot(p[i], p[j]);
        return pair.f*u;
      }

      double energy(const vector<particle> &p, const group &g1, const group &g2) {
        double u=0;
        int ilen=g1.end+1, jlen=g2.end+1;
#pragma omp parallel for reduction (+:u)
        for (int i=g1.beg; i<ilen; i++)
          for (int j=g2.beg; j<jlen; j++)
            u += pair.pairpot(p[i],p[j]);
        return pair.f*u;
      }

      /*!
       * ...between the two dipoles a and b, separated by the
       * distance r.
       * \f$ \beta u(r) = l_B \frac{a_x b_x + a_y b_y - 2a_z b_z  }{r^3}\f$
       */
      double dipdip(const point &a, const point &b, double r) {
        return pair.f*( a.x*b.x + a.y*b.y - 2*a.z*b.z )/(r*r*r);
      }
      double iondip(const point &a, double q, double r) { return -pair.f*q*a.z/(r*r); }

      // Total electrostatic potential in a point
      double pot(const vector<particle> &p, const point &a) {
        double u=0;
        int i,n=p.size();  
        for (i=0; i<n; i++)
          u+=p[i].charge/p[i].dist(a);
        return pair.f*u;
      }

      // Internal (NON)-electrostatic energy in group
      double internal(const vector<particle> &p, const group &g) {
        if (g.beg==-1) return 0;
        double u=0;
        int i,j,glen=g.end+1;
        for (i=g.beg; i<glen-1; i++)
          for (j=i+1; j<glen; j++)
            u+=pair.pairpot(p[i],p[j]);
        return pair.f*u;
      }

      double energy(const vector<particle> &p, const particle &a) {
        double u=0;
        int i,n=p.size();
        for (i=0; i<n; i++)
          u+=pair.pairpot(p[i], a);
        return pair.f*u;
      }

      /*! \note If the charge of the j'th particle is 0, ZERO will be returned!
       *  \return \f$ \phi_j = \sum_{i\neq j}^{N} \frac{l_B z_i}{r_{ij}} \f$
       *  \param p Particle vector
       *  \param j The electric potential will be calculated in the point of this particle
       *  \warning Uses simple particle->dist() function!
       *  \todo Respect cell boundaries
       */
      double potential(const vector<particle> &p, unsigned short j) {
        if (p[j].charge==0) return 0;
        double u=0;
        int i,n=p.size();
        for (i=0; i<j; ++i) u+=p[i].charge/p[i].dist(p[j]);
        for (i=j+1; i<n; ++i) u+=p[i].charge/p[i].dist(p[j]);
        return u;
      }

      string info() {
        std::ostringstream o;
        o << energybase::info()
          << "#   Pair potential:" << endl << pair.info();
        return o.str();
      }
  }; //end of interaction class

  /*!
   * \brief Hydrophobic interaction between ions and molecular surfaces
   * \author Mikael Lund
   * \date Canberra 2008
   * \todo Not optimized - inelegant "end_of_protein_one" hack. Use vector instead...
   *
   * This class will use the specified pair potential as usual but in addition add
   * a hydrophobic interaction between ions (or any specified species) and the
   * nearest hydrophobic particle. This requires an expanded pair-potential that
   * contains a function hypairpot(). If you need the ions to interact with the
   * hydrophobic groups on TWO proteins, you must set the end_of_protein_one variable.
   * In this way the minimum distance search is repeated on the remaining particles.
   */
  template<class T> class int_hydrophobic : public interaction<T> {
    private:
      vector<unsigned short> hy,pa;
      double hyenergy(const vector<particle> &);
      double hyenergy(const vector<particle> &, int);
    public:
      int_hydrophobic(const inputfile &in) : interaction<T>(in) { end_of_protein_one=int(1e7); }
      unsigned int end_of_protein_one;              //!< Last particle in protein one (set if appropriate)
      void search(const vector<particle> &);        //!< Locate hydrophobic groups and ions
      double energy(const vector<particle> &p ) { return interaction<T>::energy(p) + hyenergy(p);}
      double energy(const vector<particle> &p, int i) { return interaction<T>::energy(p,i) + hyenergy(p);}
      double energy(const  vector<particle> &p, const group &g ) { return interaction<T>::energy(p,g) + hyenergy(p);}
  };

  template<class T> void int_hydrophobic<T>::search(const vector<particle> &p) {
    pa.resize(0);
    hy.resize(0);
    for (int i=0; i<p.size(); i++)
      if (p[i].hydrophobic==true)
        hy.push_back(i);
      else if (p[i].id==particle::NA || p[i].id==particle::CL || p[i].id==particle::I)
        pa.push_back(i);
  }

  template<class T> double int_hydrophobic<T>::hyenergy(const vector<particle> &p) {
    double u=0;
#pragma omp parallel for reduction (+:u)
    for (int i=0; i<pa.size(); i++)  // loop over ions
      u+=hyenergy(p, pa[i]);                    // energy with hydrophobic groups
    return u; // in kT
  }

  template<class T> double int_hydrophobic<T>::hyenergy(const vector<particle> &p, int i) {
    if (p[i].hydrophobic==true) return 0;
    int j,hymin;
    double d,dmin=1e7,u=0;
    for (j=0; j<hy.size(); j++) {     // loop over hydrophobic groups
      if (hy[j]>end_of_protein_one) { // test if we move into second protein
        u=interaction<T>::pair.hypairpot( p[i], p[hymin], sqrt(dmin) );
        dmin=1e7;                     // reset min. distance
      }
      d=p[i].sqdist( p[hy[j]]);       // find min. distance
      if (d<dmin) {
        dmin=d;      // save min dist.
        hymin=hy[j]; // ...and particle number
      }
    }
    return interaction<T>::pair.f *
      (u + interaction<T>::pair.hypairpot( p[i], p[hymin], sqrt(dmin) ) );
  }



  /*!
   * \brief Hardsphere check, then normal potential function
   * \author Mikael Lund
   * \date Prague 2008
   * \todo Add overlap check in the system energy function.
   * \warning Untested!
   *
   * This interaction class first check for hardsphere overlap
   * and - if none found - proceeds with normal energy summation
   * according the specified pair potential.
   */
  template<class T>
    class interaction_hs : public interaction<T>, private hardsphere {
      private:
        double infty;
      public:
        interaction_hs(const inputfile &in) : interaction<T>(in) { infty=1000.; }
        double energy(const particle &a, const particle &b) {
          return (a.overlap(b)==true) ? infty  : interaction<T>::energy(a,b);
        }
        double energy(const vector<particle> &p) {
          return interaction<T>::energy(p);
        }
        double energy(const vector<particle> &p, int i) {
          return (overlap(p,i)==true) ? infty : interaction<T>::energy(p,i);
        }
        double energy(const vector<particle> &p, const particle &a) {
          return (overlap(p,a)==true) ? infty : interaction<T>::energy(p,a);
        }
        double energy(const vector<particle> &p, const group &g) {
          return (overlap(p,g)==true) ? infty : interaction<T>::energy(p,g);
        }
    };
}
#endif
