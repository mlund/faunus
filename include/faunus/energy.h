#ifndef FAU_ENERGY_H
#define FAU_ENERGY_H

#include "faunus/common.h"
#include "faunus/group.h"
#include "faunus/inputfile.h"
#include "faunus/hardsphere.h"
#include "faunus/potentials/pot_coulomb.h"
#include "faunus/potentials/pot_minimage.h"
#include "faunus/potentials/pot_hsminimage.h"
#include "faunus/potentials/pot_hscoulomb.h"
#include "faunus/potentials/pot_hydrophobic.h"
#include "faunus/potentials/pot_datapmf.h"
#include "faunus/potentials/pot_table.h"
#include "faunus/potentials/pot_debyehuckel.h"
#include "faunus/potentials/pot_debyehuckelP3.h"
#include "faunus/potentials/pot_barecoulomb.h"
#include "faunus/potentials/pot_netz.h"
#include "faunus/potentials/pot_test.h"
#include "faunus/potentials/pot_hypersphere.h"

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
      virtual double energy(vector<particle> &, const particle &)=0;               //!< all<->external particle
      virtual double energy(vector<particle> &, int)=0;                            //!< all<->particle i.
      virtual double energy(vector<particle> &, const group &)=0;                  //!< all<->group.
      virtual double energy(vector<particle> &)=0;                                 //!< all<->all (System energy).
      virtual double energy(vector<particle> &, const group &, const group &)=0;   //!< group<->group.
      virtual double energy(vector<particle> &, const group &, int)=0;             //!< group<->particle i.
      virtual double energy(vector<particle> &, const group &, const particle &)=0;//!< group<->external particle.
      virtual double energy(vector<particle> &, molecules &, vector<int> &)=0;     //!< subset[molecules]<->all
      virtual double potential(vector<particle> &, unsigned short)=0;              //!< Electric potential at j'th particle
      virtual double potential(vector<particle> &, point)=0;                       //!< Electric potential in point
      virtual double potential(const vector<particle>&,unsigned short, macromolecule&, double& )=0; //! For the AT titration
      virtual double internal(vector<particle> &, const group &, int=1)=0;         //!< internal energy in group
      virtual double pot(vector<particle> &, const point &)=0;                     //!< Electrostatic potential in a point
      virtual double dipdip(const point &, const point &, double)=0;                     //!< Dipole-dipole energy.
      virtual double iondip(const point &, double, double)=0;                            //!< Ion-dipole energy.
      virtual double force(container &, particle, particle, point, double, double=.5)=0;        //!< Force vector
      virtual void forceall(container &, vector<point> &)=0;
      //virtual void set_kappa( double& )=0;
      virtual double u_monomer(vector<particle> &, const polymer &, unsigned int)=0; //!< all<->monomer in polymer
      virtual double uself_polymer(vector<particle> &, const polymer&)=0;          //!< internal polymer energy

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
      interaction(inputfile &in) : pair(in), energybase(pair.f) {
        name="Full N^2";
        tokT=pair.f;
      };

      double energy(const particle &a, const particle &b) {
        return pair.pairpot(a,b) * pair.f;
      }

      double energy(vector<particle> &p, int j) {
        int n=p.size();
        double u=0;
        particle bak=p[j]; //avoid self interaction
        p[j].charge=0;     //in a slightly dirty manner!
        p[j].y=1e9;
#pragma omp parallel for reduction (+:u)
        for (int i=0; i<n; ++i)
          u+=pair.pairpot( p[i],p[j] );
        p[j]=bak;
        return pair.f*u;
      }

      double energy(vector<particle> &p, const group &g) {
        int n=g.end+1, psize=p.size();
        double u=0;
#pragma omp parallel for reduction (+:u)
        for (int i=g.beg; i<n; ++i) {
          for (int j=0; j<g.beg; j++)
            u += pair.pairpot(p[i],p[j]);
          for (int j=n; j<psize; j++)
            u += pair.pairpot(p[i],p[j]);
        }
        return pair.f*u;
      }

      double energy(vector<particle> &p, const group &g, int j) {
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

      double energy(vector<particle> &p, const group &g, const particle &a) {
        if (g.beg==-1) return 0;
        double u=0;
        int i,n=g.end+1;
        //#pragma omp parallel for reduction (+:u)
        for (i=g.beg; i<n; i++)
          u+=pair.pairpot(a, p[i]); 
        return pair.f*u;
      }

      double energy(vector<particle> &p, const vector<macromolecule> &g) {
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

      double energy(vector<particle> &p) {
        double u=0;
        int n = p.size();
#pragma omp parallel for reduction (+:u) schedule (dynamic)
        for (int i=0; i<n-1; ++i)
          for (int j=i+1; j<n; ++j)
            u+=pair.pairpot(p[i], p[j]);
        return pair.f*u;
      }

      double energy(vector<particle> &p, const group &g1, const group &g2) {
        double u=0;
        int ilen=g1.end+1, jlen=g2.end+1;
//#pragma omp parallel for reduction (+:u) schedule (dynamic)
        for (int i=g1.beg; i<ilen; i++)
          for (int j=g2.beg; j<jlen; j++)
            u += pair.pairpot(p[i],p[j]);
        return pair.f*u;
      }
      double energy(vector<particle> &p, molecules &m, vector<int> &i) {
        double du=0;
        return du;
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
      double pot(vector<particle> &p, const point &a) {
        double u=0;
        int i,n=p.size();  
        for (i=0; i<n; i++)
          u+=p[i].charge/sqrt(pair.sqdist(a,p[i]));
        return pair.f*u;
      }

      virtual double internal(vector<particle> &p, const group &g, int step=1) {
        if (g.beg==-1) return 0;
        double u=0;
        int n=g.end+1;
        for (int i=g.beg; i<n-step; i++)
          for (int j=g.beg+step*((i-g.beg)/step+1); j<n; j++)
            u+=pair.pairpot(p[i],p[j]);
        return pair.f*u;
      }

      double energy(vector<particle> &p, const particle &a) {
        double u=0;
        int i,n=p.size();
        for (i=0; i<n; i++)
          u+=pair.pairpot(p[i], a);
        return pair.f*u;
      }

      /*! 
       *  \return \f$ \phi_j = \sum_{i\neq j}^{N} \frac{l_B z_i}{r_{ij}} \f$
       *  \param p Particle vector
       *  \param j The electric potential will be calculated in the point of this particle
       */
      double potential(vector<particle> &p, unsigned short j) {
        double phi=0;
        int n=p.size();
        particle tmp=p[j]; // avoid self-interaction
        p[j].charge=0;     // in a slightly dirty manner
        p[j].y=1e9;        // (for symmetric parallel load)
#pragma omp parallel for reduction (+:phi)
        for (int i=0; i<n; ++i)
          phi+=p[i].charge/sqrt(pair.sqdist(p[i],tmp));
        p[j]=tmp;
        return pair.f*phi;
      }

      double potential(vector<particle> &p, point a) {
        double phi=0;
        for (int i=0; i<p.size(); ++i)
          phi+=p[i].charge/sqrt( pair.sqdist(a,p[i]) );
        return pair.f*phi;
      }


    double potential(const vector<particle> &p, unsigned short j, macromolecule& g, double& k) {
        if (abs(p[j].charge)<1e-6)
          return 0;
        double u=0;
        double r;
        int i,n=p.size();
        
        /* Inside the protein */
        for (i=g.beg; i<j; ++i)
          u+=p[i].charge/sqrt(pair.sqdist(p[i],p[j]));
        for (i=j+1; i<=g.end; ++i)
          u+=p[i].charge/sqrt(pair.sqdist(p[i],p[j]));
   
        
        /* Outside the protein */
          if (g.beg>0) {
            for (int i=0 ; i<g.beg ; i++) {
              r=sqrt(pair.sqdist(p[i],p[j]));
              u+=p[i].charge/r*exp(-k*r);
            };
          };
          for (int i=g.end+1; i<n; i++) {
            r=sqrt(pair.sqdist(p[i],p[j]));
            u+=p[i].charge/r*exp(-k*r);
          };
        
        return pair.f*u;
      }


      void forceall(container &c, vector<point> &f) {
        point r;
        double ff;
        int n=c.p.size();
        f.resize(n);
        for (int i=0; i<n; i++)
          f[i].clear();
        for (int i=0; i<n-1; i++)
          for (int j=i+1; j<n; j++)
          {
            r=c.vdist(c.p[i], c.p[j]);
            ff=c.p[i].charge*c.p[j].charge/(r.x*r.x+r.y*r.y+r.z*r.z);
            f[i].x+=ff*r.x;
            f[i].y+=ff*r.y;
            f[i].z+=ff*r.z;
            f[j].x-=ff*r.x;
            f[j].y-=ff*r.y;
            f[j].z-=ff*r.z;
          }
      }
      
      /*
      void set_kappa( double& kappa) {
      	try { pair.k = kappa; }
      	catch (...) {};
      };
      */

      double force(container &c, particle a, particle b, point rij, double r, double dr) {
        double forward,center,f;
        point unit;
        unit=rij, forward=1./r, unit=unit*forward*dr;  // Dirty lending of variable
        a.x=a.y=a.z=0;
        b.x=rij.x+unit.x, b.y=rij.y+unit.y, b.z=rij.z+unit.z;
        forward=pair.pairpot(a,b);
        b.x=rij.x-unit.x, b.y=rij.y-unit.y, b.z=rij.z-unit.z;
        center=pair.pairpot(a,b);
        f= -(forward-center)/(2*dr);
        return f;
      }
      ;
      double u_monomer(vector<particle> &p, const polymer &g, unsigned int i) {
        return 0;
      }
      double uself_polymer(vector<particle> &p, const polymer &g) {
        return 0;
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
      double hyenergy(vector<particle> &);
      double hyenergy(vector<particle> &, int);
    public:
      int_hydrophobic(inputfile &in) : interaction<T>(in) { end_of_protein_one=int(1e7); }
      unsigned int end_of_protein_one;              //!< Last particle in protein one (set if appropriate)
      void search(vector<particle> &);        //!< Locate hydrophobic groups and ions
      double energy(vector<particle> &p ) { return interaction<T>::energy(p) + hyenergy(p);}
      double energy(vector<particle> &p, int i) { return interaction<T>::energy(p,i) + hyenergy(p);}
      double energy(const  vector<particle> &p, const group &g ) { return interaction<T>::energy(p,g) + hyenergy(p);}
  };

  template<class T> void int_hydrophobic<T>::search(vector<particle> &p) {
    pa.resize(0);
    hy.resize(0);
    for (int i=0; i<p.size(); i++)
      if (p[i].hydrophobic==true)
        hy.push_back(i);
      else if (p[i].id==0 || p[i].id==1 || p[i].id==3)
        pa.push_back(i);
    //else if (p[i].id==particle::NA || p[i].id==particle::CL || p[i].id==particle::I)
    //  pa.push_back(i);
  }

  template<class T> double int_hydrophobic<T>::hyenergy(vector<particle> &p) {
    double u=0;
    int n=pa.size();
#pragma omp parallel for reduction (+:u)
    for (int i=0; i<n; i++)  // loop over ions
      u+=hyenergy(p, pa[i]);                    // energy with hydrophobic groups
    return u; // in kT
  }

  template<class T> double int_hydrophobic<T>::hyenergy(vector<particle> &p, int i) {
    if (p[i].hydrophobic==true) return 0;
    int j,hymin=0;
    double d,dmin=1e7,u=0;
    for (j=0; j<hy.size(); j++) {     // loop over hydrophobic groups
      if (hy[j]>end_of_protein_one) { // test if we move into second protein
        u=interaction<T>::pair.hypairpot( p[i], p[hymin], sqrt(dmin) );
        dmin=1e7;                     // reset min. distane
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
        interaction_hs(inputfile &in) : interaction<T>(in) { infty=1000.; }
        double energy(const particle &a, const particle &b) {
          return (a.overlap(b)==true) ? infty  : interaction<T>::energy(a,b);
        }
        double energy(vector<particle> &p) {
          return interaction<T>::energy(p);
        }
        double energy(vector<particle> &p, int i) {
          return (overlap(p,i)==true) ? infty : interaction<T>::energy(p,i);
        }
        double energy(vector<particle> &p, const particle &a) {
          return (overlap(p,a)==true) ? infty : interaction<T>::energy(p,a);
        }
        double energy(vector<particle> &p, const group &g) {
          return (overlap(p,g)==true) ? infty : interaction<T>::energy(p,g);
        }
    };

  template<class T>
    class interaction_vector : public interaction<T> {
      private:
        int len;
        double infty;
        double r2[4000], qq[4000];
      public:
        interaction_vector(inputfile &in) : interaction<T>(in) {
          interaction<T>::name+="Vectorized, Full N^2";
        }
        double energy(vector<particle> &p) {
          len=1;
          int n=p.size();
          for (int i=0; i<n-1; i++)
            for (int j=i+1; j<n; j++) {
              r2[len]   = interaction<T>::pair.sqdist(p[i],p[j]);
              qq[len++] = p[i].charge*p[j].charge;
            }
          return interaction<T>::pair.VectorEnergy(r2,qq,&len);
        }
        double energy(vector<particle> &p, const group &g) {
          len=1;
          int n=g.end+1, psize=p.size();
          for (int i=g.beg; i<n; ++i) {
            for (int j=0; j<g.beg; j++) {
              r2[len]   = interaction<T>::pair.sqdist(p[i],p[j]);
              qq[len++] = p[i].charge*p[j].charge;
            }
            for (int j=n; j<psize; j++) {
              r2[len]   = interaction<T>::pair.sqdist(p[i],p[j]);
              qq[len++] = p[i].charge*p[j].charge;
            }
          }
          return interaction<T>::pair.VectorEnergy(r2,qq,&len);
        }
        double energy(vector<particle> &p, const group &g1, const group &g2) {
          len=1;
          int ilen=g1.end+1, jlen=g2.end+1;
          //#pragma omp parallel for reduction (+:u)
          for (int i=g1.beg; i<ilen; i++)
            for (int j=g2.beg; j<jlen; j++) {
              r2[len]   = interaction<T>::pair.sqdist(p[i],p[j]);
              qq[len++]= p[i].charge*p[j].charge;
            }
          return interaction<T>::pair.VectorEnergy(r2,qq,&len);
        }
    };

  /*!
   * \brief Treats far-away groups as monopoles for faster energy evaluation
   * \author Mikael Lund
   * \date Lund 2009
   * \warning Be careful when you have non-molecular groups such as salt.
   *
   * This interaction class redefines group-group and group-particle interactions
   * so that groups (i.e. molecules) far away will be seen as monopoles centered
   * at their charge-centers. I.e. the interaction between two charged proteins beyond
   * some threshold distance will interact as two point charges. Likewise a particle
   * far way from a group will see it only as a single particle
   */
  template<class T>
    class interaction_monopole : public interaction<T> {
      private:
        particle monopole(vector<particle> &p, const group &g) {
          double sum=0;
          particle cm;
          point t, o = p[g.beg]; // set origo to first particle
          for (unsigned short i=g.beg; i<=g.end; i++) {
            t = p[i]-o;        // translate to origo
            interaction<T>::pair.boundary(t);       // periodic boundary (if any)
            cm += t * p[i].mw;
            sum += p[i].mw; 
            cm.charge+=p[i].charge;
          }
          cm=cm*(1./sum) + o;
          interaction<T>::pair.boundary(cm);
          return cm;
        }
      public:
        container *cPtr;
        double cut_g2g; //!< Cut-off distance for group-group interactions
        double cut_g2p; //!< Cut-off distance for group-particle interactions
        interaction_monopole(inputfile &in, container &con) : interaction<T>(in) {
          cPtr=&con;
          interaction<T>::name+=" w. monopole cut-offs";
          cut_g2g = in.getflt( "threshold_g2g", 1e6 );
          cut_g2p = in.getflt( "threshold_g2p", 1e6 );
        }
        double energy(vector<particle> &p ) { return interaction<T>::energy(p); } 
        double energy(vector<particle> &p, const group &g1, const group &g2) {
          particle mp1=monopole(p,g1), mp2=monopole(p,g2);
          return ( sqrt(interaction<T>::pair.sqdist(mp1, mp2)) > cut_g2g ) ?
            interaction<T>::energy( mp1, mp2 ) : interaction<T>::energy(p, g1, g2);
        }
        double energy(vector<particle> &p, const group &g) {
          double u=0;
          particle mp=monopole(p,g);
          for (int i=0; i<g.beg; i++)
            u+= (sqrt(interaction<T>::pair.sqdist(mp,p[i]))>cut_g2p) ? interaction<T>::energy(mp,p[i]) : interaction<T>::energy(p,g,i);
          for (int i=g.beg+1; i<p.size(); i++)
            u+= (sqrt(interaction<T>::pair.sqdist(mp,p[i]))>cut_g2p) ? interaction<T>::energy(mp,p[i]) : interaction<T>::energy(p,g,i);
          return u;
        }
        string info() {
          std::ostringstream o;
          o << interaction<T>::info()
            << "#   Group-Group threshold    = " << cut_g2g << endl
            << "#   Group-Particle threshold = " << cut_g2p << endl;
          return o.str();
        }
    };

  /*!
   * \brief Interaction class that includes image charges outside a spherical cell
   * \todo Not optimized
   * \warning Untested! (behaves properly but not independently confirmed correct)
   *          This class is only INTENDED to be used with the potetnential test_pot!
   * \author Mikael Lund, Bjorn Persson
   * \date Lund 2008-2009
   *
   * Calculates inter-particle interactions in a dielectric sphere (epsi) and
   * adds image charge interactions with the dielectric surroundings (epso). Useful
   * for simulating explicit water in a spherical container. Note that after each
   * MC move -- both accepted and rejected -- the img vector *must* be updated
   * with the updateimg() function. The cavity origin is assumed to be 0,0,0. Note that
   * image charge and distance are divergent at (0,0,0). Also keep in mind that the energy
   * is NOT the product of the potential at a charge times the charge it self but rather 
   * half of the product since it is an reaction potential (as what is conserned for
   * the reaction part). There are there for additional 
   * functions sphericalimage::elenergy to obtain the electrostatic energy for particles 
   * and groups. The functions sphericalimage::potential returns the total potential.
   *
   * For reference: Harold L Friedman, Molecular Physics, 1975 vol. 29 ppg. 1533-1543.
   * Extended method: Wei Cai, XXXX.
   */
  template<class T> class sphericalimage : public interaction<T> {
    private:
      average<double> ratio;
      double scale, radius, radius2, epso, epsi, ui, ur;
      vector<point> img;   //!< Contain image charge particles
      vector<double> ich;  //!< Complementary vector for image charges
    public:
      sphericalimage(inputfile &in) : interaction<T>(in) {
        interaction<T>::name+=" w. spherical image charges";
        epso=in.getflt("epso",80);   //permitivity of outside medium
        epsi=in.getflt("epsi",1);  //permitivity of inner medium
        radius=in.getflt("cellradius",0); 
        radius2=radius*radius;
        scale=-(epso-epsi)/(epso+epsi)*radius/2*in.getflt("bjerrum", 560.2);  // Friedman 1975, Mol. Phys 29 pp. 1533
      }                                                                 // warning!!! temperature units are given through scale!

      // RE-CALC IMAGE POSITIONS
      inline void updateimg(const particle &a, int i) {
        img[i] = a* (radius2 / a.dot(a)); 
        ich[i] = a.charge/a.len();
      }
      void updateimg(vector<particle> &p) {
        img.resize( p.size() );
        ich.resize( p.size() );
        for (int i=0; i<img.size(); ++i)
          updateimg(p[i],i);
      }
      void updateimg(vector<particle> &p, const group &g) {
        for (int i=g.beg; i<=g.end; ++i)
          updateimg(p[i],i);
      }
      void updateimg(vector<particle> &p, molecules &m, vector<int> &n) {
        for (int i=0; i<n.size(); i++) {
          updateimg(p,m[i]);
          updateimg(p,m[n[i]]);
        }
      }
      // IMAGE POTENTIAL  (half the potential)
      double impot(double &ich, const point &pr, point &pi) {
        return ich/pr.dist(pi);
      }

      // IMAGE ENERGY
      double image(vector<particle> &p) {
        int n=p.size();
        double u=0;
#pragma omp parallel for reduction (+:u) schedule (dynamic)
        for (int i=0; i<n; ++i)
          u+=image(p,i);
        return u;
      }
      double image(vector<particle> &p, int i) {
        double u=0;
        int t=img.size();
#pragma omp parallel for reduction (+:u) schedule (dynamic)
        for (int j=0; j<t; ++j)
          u += impot(ich[j], p[i], img[j] );
        return p[i].charge*scale*u*2. - p[i].charge*scale*impot(ich[i],p[i],img[i]) ;        //make sure too and not too double count
      }
      double image(vector<particle> &p, const group &g) {
        double u=0;
        for (int i=g.beg; i<=g.end; ++i)
          u+=image(p, i, g.beg, g.end);
        return u;
      }
      double image(vector<particle> &p, int i, int j, int k) {
        int n=p.size();
        double uin=0,uex=0;
#pragma omp parallel for reduction (+:uex) schedule (dynamic)
        for (int s=0; s<j; s++)
          uex += impot(ich[s], p[i], img[s] );  //make sure to double count
#pragma omp parallel for reduction (+:uex) schedule (dynamic)
        for (int t=k+1; t<n; t++)
          uex += impot(ich[t], p[i], img[t] );  //make sure to double count
        for (int u=j; u<=k; u++)
          uin += impot(ich[u], p[i], img[u] );  // internal interactions will be double counted implicitly
                                                // not paralleized on purpose
        // the self term will not be double counted
        return p[i].charge*scale*(uex*2+uin); 
      }
      double image(vector<particle> &p, const group &g1, const group &g2) {
        double u=0;
        for (int i=g1.beg; i<=g1.end; i++)               //Dielectric and g1
          for (int j=g1.beg; j<=g1.end; j++)
            u += p[i].charge*impot(ich[j], p[i], img[j]);
        for (int i=g2.beg; i<=g2.end; i++)               //Delectric and g2
          for (int j=g2.beg; j<=g2.end; j++)
            u += p[i].charge*impot(ich[j], p[i], img[j]);
        for (int i=g1.beg; i<=g1.end; i++)               //g1 and g2 through dielectric
          for (int j=g2.beg; j<=g2.end; j++)
            u += 2*p[i].charge*impot(ich[j], p[i], img[j]);
        return u*scale;
      } 
      double imageint(vector<particle> &p, group) {
        double u=0;
        return u;
      }
      // TOTAL ENERGY
      double potential(vector<particle> &p, int i) {
        int n=p.size();
        double ur=0,ui=0;
        updateimg(p[i],i);
        ur=interaction<T>::potential(p,i);
#pragma omp parallel for reduction (+:ui) schedule (dynamic)
        for (int s=0; s<n; s++)
          ui += impot(ich[s], p[i], img[s] );  
        ui*=2;
        return ur+ui*scale;  // This is the POTENTIAL, this should not be used to calculate the 
      }                      // interaction since that would 'double' count the self term.
      double potential(vector<particle> &p, point i) {
        ur=ui=0;
        ur=interaction<T>::potential(p,i);
//#pragma omp parallel for reduction (+:ui) schedule (dynamic)
        for (int s=0; s<p.size(); s++)
          ui += impot(ich[s], i, img[s] );  
        ui*=2;               // Due to the definition of scale
        return ur+ui*scale;  // This is the POTENTIAL in any given point inside the cavity
      }                     
      double elenergy(vector<particle> &p, int i) {
        ur=ui=0;
        updateimg(p[i],i);
        ur=p[i].charge*interaction<T>::potential(p,i);
        ui=image(p,i);
        return ur+ui;        // Returns the electrostatic interaction energy of i with p
      }
      double elenergy(vector<particle> &p, const group &g) {
        ur=ui=0;
        updateimg(p,g);
        for (int i=g.beg; i<=g.end; i++) {
          ur+=p[i].charge*interaction<T>::potential(p,i);
          ui+=image(p,i);
        }
        return ur+ui;        // Returns the electrostatic energy of g with p
      }
      double energy(vector<particle> &p ) {
        updateimg(p);
        ur=interaction<T>::energy(p);
        ui=image(p);
        ratio+=std::abs(ui/(ur+ui));
        return ur+ui;
      }
      double energy(vector<particle> &p, int i) {
        updateimg(p[i],i);
        ur=interaction<T>::energy(p,i);
        ui=image(p,i);
        ratio+=std::abs(ui/(ur+ui));
        return ur+ui;
      }
      double energy(vector<particle> &p, const group &g ) {
        updateimg(p,g); // update all images in group
        ur=interaction<T>::energy(p,g);
        ui=image(p,g);
        ratio+=std::abs(ui/(ur+ui));
        return ur+ui;
      }
      double energy(vector<particle> &p, const group&g1, const group &g2) {
        updateimg(p,g1), updateimg(p,g2);
        ur=interaction<T>::energy(p, g1, g2);
        ui=image(p,g1,g2);
        return ur+ui;
      }
      double energy(vector<particle> &p, molecules &m, vector<int> &n) {
        ur=ui=0;
        group g;
        g.beg=m[0].beg;
        g.end=m[n.size()-1].end;
        updateimg(p,m,n);
        ur=interaction<T>::energy(p,g);
        for (int i=0; i<n.size()-1; i++)
          for (int j=i+1; j<n.size(); j++)
            ur+=interaction<T>::energy(p,m[i],m[j]);
        //for (int i=0; i<n.size(); i++)
        //  ur+=interaction<T>::energy(p, m[i]);
        ui=image(p,g);
        ratio+=std::abs(ui/(ur+ui));
        return ur+ui;
      }
      // INFO
      string info() {
        std::ostringstream o;
        o << interaction<T>::info()
          << "#   Dielectric const. (in, out) = " << epsi << " " << epso << endl
          << "#   Cavity radius               = " << radius << endl
          << "#   Avg. image energy ratio     = " << ratio.avg() <<" , stdev "<< ratio.stdev()<< endl
          << "#   Number of images            = " << img.size() << endl
          << "#   Scaling const.              = " << scale << endl
          << "#      -(epso-1)/(epso+1)*radius/2*lb "<<endl;
        return o.str();
      }
      string printimg() {
        std::ostringstream o;
        for (int i =0; i<img.size(); i++) 
          o << img[i] <<"  "<<ich[i]<<endl;
        o << endl;
        return o.str();
      }
  };

}//namespace
#endif
