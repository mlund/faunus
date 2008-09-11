#ifndef FAU_potential_h
#define FAU_potential_h

#include <cmath>
#include "species.h"
#include "xydata.h"
#include "group.h"
#include "lennardjones.h"
#include "inputfile.h"

namespace Faunus {
  /*!
   *  \brief Setup for potentials.
   *  \author Mikael Lund
   *  \date Prague, 2007
   *
   *  This class is used to pass parameters to classes
   *  that handles particle pair-potentials.
   */
  class pot_setup {
    public:
      pot_setup();
      double kappa,        //!< Inverse Debye screening length
             lB,           //!< Bjerrum length
             eps,          //!< L-J parameter
             r0,           //!< Bond eq. distance
             epsi,         //!< Internal dielectric constant
             epso,         //!< External dielectric constant
             hydroscale,   //!< LJ scaling factor for hydrophobic interactions
             box,          //!< Cubic box length
             a,            //!< Cavity radius
             A,B,C,D;      //!< Empirical parameters for "pot_netz"
      pot_setup(inputfile &);
  };

  /*! \brief Coulomb potential
   *  \author Mikael Lund
   *  \date Prague, 2007
   */
  class pot_hscoulomb : public pot_hs {
    public:
      /*! \param pot.lB Bjerrum length
       *  \param pot.eps L-J epsilon parameter (in kT) */
      pot_hscoulomb ( pot_setup &pot) : pot_hs() {
        f=pot.lB;
        name+="/Coulomb";
      }
      /*! \brief Return Coulomb energy between a pair of particles
       *  \return Energy in units of kT/f (f=lB).
       *  \f$ \beta u/f = \frac{z_1 z_2}{r} + u_{HS}/f \f$
       */
      inline double pairpot(particle &p1, particle &p2) {
        register double r2=p1.sqdist(p2);
        return hs(p1,p2,r2) + p1.charge*p2.charge/sqrt(r2);
      }
      string info();
  };


  /*! \brief Coulomb potential
   *  \author Mikael Lund
   *  \date Prague, 2007
   */
  class pot_coulomb : public pot_lj {
    public:
      /*! \param pot.lB Bjerrum length
       *  \param pot.eps L-J epsilon parameter (in kT) */
      pot_coulomb ( pot_setup &pot) : pot_lj(pot.eps/pot.lB) {
        f=pot.lB;
        name+="/Coulomb";
      }
      /*! \brief Return Coulomb energy between a pair of particles
       *  \return Energy in units of kT/f (f=lB).
       *  \f$ \beta u/f = \frac{z_1 z_2}{r} + u_{LJ}/f \f$
       */
      inline double pairpot(particle &p1, particle &p2) {
        register double r2=p1.sqdist(p2);
        return lj(p1,p2,r2) + p1.charge*p2.charge/sqrt(r2);
      }
      string info();
  };

  /*!
   * \brief Coulomb pot. with minimum image.
   * \author Mikael Lund
   * \date 2007
   */
  class pot_minimage : public pot_lj {
    private:
      double invbox,box;
    public:
      pot_minimage(pot_setup &pot) : pot_lj(pot.eps/pot.lB) {
        f=pot.lB;
        box=pot.box;
        invbox=1./box;
        name+="/Coulomb w. minimum image";
      }
      string info();
      void setvolume(double vol) {
        box=pow(vol, 1./3);;
        invbox=1./box;
      }
      inline double pairpot(particle &p1, particle &p2) {
        register double r2=p1.sqdist(p2,box,invbox);
        return lj(p1,p2,r2) + p1.charge*p2.charge/sqrt(r2);
      }
  };

  class pot_test : public pot_lj {
    public:
      pot_test( pot_setup &pot ) : pot_lj(pot.eps/pot.lB) { f=pot.lB; }
      string info();
      inline double pairpot(particle &p1, particle &p2) {
        register double r2=p1.sqdist(p2);
        register double a=p1.radius+p2.radius;
        a=a*a;
        if (r2<4*a) {
          a=a/r2;
          a=a*a*a;
          a=a*a/f;
        } else a=0;
        register double qq=p1.charge*p2.charge;
        return (qq!=0) ? qq/sqrt(r2)+a : a;
      }
  };

  /*! \brief Debye-Huckel potential
   *  \author Mikael Lund
   */
  class pot_debyehuckel : public pot_lj {
    private:
      double k;
    public:
      //! \param pot.lB Bjerrum length
      //! \param pot.eps L-J parameter
      //! \param pot.kappa Inverse Debye screening length
      pot_debyehuckel( pot_setup &pot ) : pot_lj(pot.eps/pot.lB) {
        f=pot.lB; 
        k=pot.kappa; 
        name+="/Debye-Huckel";
      };
      string info();
      //! \f$ \beta u/f = \frac{z_1z_2}{r}\exp(-\kappa r) + u_{lj}/f \f$
      //! \return Energy in kT/f (f=lB)
      inline double pairpot( particle &p1, particle &p2 ) {
        double r2=p1.sqdist(p2),
               r=sqrt(r2);
        return lj(p1,p2,r2) + p1.charge*p2.charge/r*exp(-k*r);
      }
  };

  /*! \brief Debye-Huckel potential for periodic boundry 
   *         conditions in 3D, it is extended to preform 
   *         under conditions of constant pressure.
   See class isobaric->markovmove.h
   *  \author Mikael Lund/Bjoern Persson
   *  \date Lund/Prag 2008
   */
  class pot_debyehuckelP3 : public pot_lj {
    private:
      double k;
    public:
      double box, invbox;
      //! \param pot.lB Bjerrum length
      //! \param pot.eps L-J parameter
      //! \param pot.kappa Inverse Debye screening length
      pot_debyehuckelP3( pot_setup &pot ) : pot_lj(pot.eps/pot.lB) {
        f=pot.lB; 
        k=pot.kappa;
        box=pot.box;
        invbox=1./box;
        name+="/Debye-Huckel w. minimum image";
      };
      string info();
      //! \f$ \beta u/f = \frac{z_1z_2}{r}\exp(-\kappa r) + u_{lj}/f \f$
      //! \return Energy in kT/f (f=lB)
      inline double pairpot( particle &p1, particle &p2 ) {
        register double r2=p1.sqdist(p2,box,invbox), r=sqrt(r2);
        return lj(p1,p2,r2) + p1.charge*p2.charge/r*exp(-k*r);
      }
      void setvolume(double vol) {
        box=pow(vol, 1./3);;
        invbox=1./box;
      }
  }; 

  /*!
   * \brief Load pair potentials from disk
   * \author Mikael Lund
   * \date Prague, 2007
   *
   * Class to load potential of mean force (PMF) data from file(s). If
   * the distance between the particles is outside the data range a 
   * simple Coulomb + LJ potential will be applied.
   */
  class pot_datapmf : public pot_lj {
    private:
      xydata<double> pmfd[particle::LAST][particle::LAST];
    public:
      //! \param pot.lB Bjerrum length
      //! \param pot.eps L-J parameter
      pot_datapmf(pot_setup &pot) : pot_lj(pot.eps/pot.lB) {
        f=pot.lB;
        name+="/Empirical data potential";
      }
      bool loadpmf(species &, string);          // load pmf's from disk
      void loadpmf(species &, string,string);   // -//-
      void showpmf(species &);                  //!< Lists loaded pmf's.
      string info();
      double pairpot (particle &p1, particle &p2) {
        unsigned short i=p1.id,j=p2.id;
        if (i>j) swap(i,j);
        double r2=p1.sqdist(p2);
        if (pmfd[i][j].xmax==0) {               // if no data exists:
          double u,x=p1.charge*p2.charge;       // use Coulomb + lj pot.
          u=lj(p1,p2,r2);
          return (x!=0) ? u+x/sqrt(r2) : u; 
        }; 
        r2=sqrt(r2);
        return (r2>pmfd[i][j].xmax) ? 
          p1.charge*p2.charge/r2 : // use Coulomb pot. outside data 
          pmfd[i][j].x2y(r2);      // ...else use table. 
      }
  };

  /*!
   *  \brief Interaction between particles and groups
   *  \author Mikael Lund
   *  \param T pairpotential
   *
   *  Calculates interaction energies between particles and groups. The
   *  pair potential is specified as a template type which allows inlining.
   *  Unless otherwise specified, all energies will be returned in units of \b kT.
   */
  template<class T>
    class interaction {
      public:
        T pair;     //!< Pair potential class
        interaction(pot_setup &pot) : pair(pot) {};
        string info();
        double energy(vector<particle> &, particle &);            //!< all<->external particle
        virtual double energy(vector<particle> &, int);           //!< all<->particle i.
        virtual double energy(vector<particle> &, group &);               //!< all<->group.
        virtual double energy(vector<particle> &);                          //!< all<->all (System energy).
        double energy(vector<particle> &, group &, group &);     //!< group<->group.
        double energy(vector<particle> &, group &, int);            //!< group<->particle i.
        double energy(vector<particle> &, group &, particle &);     //!< group<->external particle.
        double energy(vector<particle> &, vector<group> &, int,...);
        double energy(vector<particle> &, int, vector<short int> &);//!< particle<->list of particles.
        double energy(vector<particle> &, vector<macromolecule> &); //!< vector<group> <-> vector<group>
        double potential(vector<particle> &, unsigned short);       //!< Electric potential at j'th particle
        double internal(vector<particle> &, group &);               //!< internal energy in group
        double pot(vector<particle> &, point &);                    //!< Electrostatic potential in a point
        double quadratic(point &, point &);
        double graft(vector<particle> &, group &);
        double chain(vector<particle> &, group &, int);
        double dipdip(point &, point &, double);                    //!< Dipole-dipole energy.
        double iondip(point &, double, double);                     //!< Ion-dipole energy.
    };

  template<class T>
    double interaction<T>::energy(vector<particle> &p, int j) {
      int i,ps=p.size();
      double u=0;
#pragma omp parallel for reduction (+:u) num_threads(2)
      for (i=0; i<j; i++)
        u+=pair.pairpot( p[i],p[j] );

#pragma omp parallel for reduction (+:u) num_threads(2)
      for (i=j+1; i<ps; i++)
        u+=pair.pairpot( p[i],p[j] );
      return pair.f*u;
    }
  template<class T> double interaction<T>::energy(vector<particle> &p, group &g) {
    int i,j,n=g.end+1, psize=p.size();
    double u=0;
#pragma omp parallel for reduction (+:u) num_threads(2)
    for (i=g.beg; i<n; i++) {
      for (j=0; j<g.beg; j++)
        u += pair.pairpot(p[i],p[j]);
      for (j=n; j<psize; j++)
        u += pair.pairpot(p[i],p[j]);
    };
    return pair.f*u;
  }
  template<class T> double interaction<T>::energy(vector<particle> &p, group &g, int j) {
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
  template<class T> double interaction<T>::energy(vector<particle> &p, group &g, particle &a) {
    if (g.beg==-1) return 0;
    double u=0;
    int i,n=g.end+1;
    //#pragma omp parallel for reduction (+:u)
    for (i=g.beg; i<n; i++)
      u+=pair.pairpot(a, p[i]); 
    return pair.f*u;
  }
  template<class T> double interaction<T>::energy(vector<particle> &p, vector<macromolecule> &g) {
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
  template<class T> double interaction<T>::energy(vector<particle> &p) {
    double u=0;
    int n = p.size();
#pragma omp parallel for reduction (+:u) num_threads(2)
    for (int i=0; i<n-1; ++i)
      for (int j=i+1; j<n; ++j)
        u+=pair.pairpot(p[i], p[j]);
    return pair.f*u;
  }
  template<class T> double interaction<T>::energy(vector<particle> &p, group &g1, group &g2) {
    double u=0;
    int i,j,ilen=g1.end+1, jlen=g2.end+1;
#pragma omp parallel for reduction (+:u) num_threads(2)
    for (i=g1.beg; i<ilen; i++)
      for (j=g2.beg; j<jlen; j++)
        u += pair.pairpot(p[i],p[j]);
    return pair.f*u;
  }

  /*!
   * ...between the two dipoles a and b, separated by the
   * distance r.
   * \f$ \beta u(r) = l_B \frac{a_x b_x + a_y b_y - 2a_z b_z  }{r^3}\f$
   */
  template<class T> double interaction<T>::dipdip(point &a, point &b, double r) {
    return pair.f*( a.x*b.x + a.y*b.y - 2*a.z*b.z )/(r*r*r);
  }
  template<class T>
    double interaction<T>::iondip(point &a, double q, double r) { return -pair.f*q*a.z/(r*r); }

  // Total electrostatic potential in a point
  template<class T> double interaction<T>::pot(vector<particle> &p, point &a) {
    double u=0;
    int i,n=p.size();  
    for (i=0; i<n; i++)
      u+=p[i].charge/p[i].dist(a);
    return pair.f*u;
  }
  // Internal (NON)-electrostatic energy in group
  template<class T> double interaction<T>::internal(vector<particle> &p, group &g) {
    if (g.beg==-1) return 0;
    double u=0;
    int i,j,glen=g.end+1;
    for (i=g.beg; i<glen-1; i++)
      for (j=i+1; j<glen; j++)
        u+=pair.pairpot(p[i],p[j]);
    return pair.f*u;
  }
  template<class T> double interaction<T>::energy(vector<particle> &p, particle &a) {
    double u=0;
    int i,n=p.size();
    for (i=0; i<n; i++)
      u+=pair.pairpot(p[i], a);
    return pair.f*u;
  }
  /*! \note If the charge of the j'th particle is 0, ZERO will be returned!
   *  \return \f$ \phi_j = \sum_{i\neq j}^{N} \frac{l_B z_i}{r_{ij}} \f$
   *  \param j The electric potential will be calculated in the point of this particle
   *  \warning Uses simple particle->dist() function!
   *  \todo Respect cell boundaries
   */
  template<class T> double interaction<T>::potential(vector<particle> &p, unsigned short j) {
    if (p[j].charge==0) return 0;
    double u=0;
    int i,n=p.size();
    for (i=0; i<j; ++i) u+=p[i].charge/p[i].dist(p[j]);
    for (i=j+1; i<n; ++i) u+=p[i].charge/p[i].dist(p[j]);
    return u;
  }
  template<class T> string interaction<T>::info() {
    ostringstream o;
    o << endl
      << "# POTENTIAL ENERGY FUNCTION:" << endl
      << pair.info();
    return o.str();
  }

  // -------------------------------------- INTERACTION W. FORCES ---------
  /*!
   * \brief Interaction class w. additional force functions
   * \author Mikael Lund
   * \date Amaroo, 2008
   * \todo Unfinished
   */
  template<class T>
    class interaction_force : public interaction<T> {
      public:
        interaction_force(pot_setup &pot) : interaction<T>(pot) {}
        point force( vector<particle> &, group & );     //!< Calculate the force acting on a group
    }; 

  template<class T> point interaction_force<T>::force( vector<particle> &p, group &g ) {
    point f = interaction<T>::pair.force(p[1],p[2]);
    return f;
  }

  // --------------------------------- HYDROPHOBIC INTERACTION ------------
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
      double hyenergy(vector<particle> &);
      double hyenergy(vector<particle> &, int);
    public:
      unsigned int end_of_protein_one;                  //!< Last particle in protein one (set if appropriate)
      int_hydrophobic(pot_setup &pot) : interaction<T>(pot) { end_of_protein_one=int(1e7); }
      void search(vector<particle> &);                  //!< Locate hydrophobic groups and ions
      double energy(vector<particle> &);                //!< all<->all
      double energy(vector<particle> &, int);           //!< all<->particle i.
      double energy(vector<particle> &, group &);       //!< all<->group.
  };
  template<class T> void int_hydrophobic<T>::search(vector<particle> &p) {
    pa.resize(0);
    hy.resize(0);
    for (int i=0; i<p.size(); i++)
      if (p[i].hydrophobic==true)
        hy.push_back(i);
      else if (p[i].id==particle::NA || p[i].id==particle::CL || p[i].id==particle::I)
        pa.push_back(i);
  }
  template<class T> double int_hydrophobic<T>::hyenergy(vector<particle> &p) {
    double u=0;
#pragma omp parallel for reduction (+:u)
    for (int i=0; i<pa.size(); i++)  // loop over ions
      u+=hyenergy(p, pa[i]);                    // energy with hydrophobic groups
    return u; // in kT
  }
  template<class T> double int_hydrophobic<T>::hyenergy(vector<particle> &p, int i) {
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
  template<class T> double int_hydrophobic<T>::energy( vector<particle> &p ) {
    return interaction<T>::energy(p) + hyenergy(p);
  }
  template<class T> double int_hydrophobic<T>::energy( vector<particle> &p, int i) {
    return interaction<T>::energy(p,i) + hyenergy(p);
  }
  template<class T> double int_hydrophobic<T>::energy( vector<particle> &p, group &g ) {
    return interaction<T>::energy(p,g) + hyenergy(p);
  }
}//namespace
#endif
