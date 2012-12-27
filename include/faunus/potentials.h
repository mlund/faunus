#ifndef FAU_POTENTIAL_H
#define FAU_POTENTIAL_H

#ifndef SWIG
#include <faunus/common.h>
#include <faunus/point.h>
#include <faunus/textio.h>
#include <faunus/physconst.h>
#include <faunus/auxiliary.h>
#include <faunus/inputfile.h>
#include <faunus/species.h>
#include <faunus/geometry.h>
#endif

namespace Faunus {

  /*!
   * \brief Namespace for pair potentials
   *
   * This namespace contains classes and templates that calculates the
   * pair potential between particles. The majority of these classes/templates
   * are derived from the base class Potential::PairPotentialBase and thus
   * have a common interface. Several pair potentials can be combined into
   * one by the template Potential::CombinedPairPotential and a number of
   * common combinations are already defined as typedefs.
   */
  namespace Potential {

    class DebyeHuckel;

    /*!
     * \brief Base class for pair potential classes
     *
     * This is a base class for all pair potentials which must implement the function
     * operator so that the potential can work as a class function.
     * To make a new pair potential you must implement
     * 1) a function that takes two particles as arguments as well as the
     *    squared distance between them (i.e. the function operator), and
     * 2) a brief information string.
     * The unit of the returned energy is arbitrary but you *must* ensure that
     * when multiplied by tokT() that it is converted to kT units (thermal energy).
     * By default _tokT=1 and it is a good policy to return energies in kT.
     *
     * \todo remove tokT() function
     */
    class PairPotentialBase {
      private:
        virtual string _brief()=0;
        virtual void _setScale(double);
      protected:
        double _tokT;
      public:  
        PairPotentialBase();
        virtual ~PairPotentialBase();
        string name;             //!< Short (preferably one-word) description of the core potential
        string brief();          //!< Brief, one-lined information string
        void setScale(double=1);
        /*! \brief Convert returned energy to kT.*/
        inline double tokT() {
          //        inline double tokT() __attribute__ ((deprecated)) {
          return _tokT;
        }
        virtual void setTemperature(double); //!< Set temperature [K]

        /*!
         * \brief Particle-particle energy in units of \c kT
         * \param a First particle
         * \param b Second particle
         * \param r2 Squared distance between them (angstrom squared)
         */
        virtual double operator() (const particle&, const particle&, double) const;
      
        /*!
         * \brief Particle-particle force in units of \c kT/Å
         * \param a First particle
         * \param b Second particle
         * \param r2 Squared distance between them (angstrom squared)
         */
        virtual Point force(const particle&, const particle&, double) {
          assert(!"Force not overrided!");
          return Point(0.0, 0.0, 0.0);
        }
      
        bool save(string, particle::Tid, particle::Tid); //!< Save table of pair potential to disk
        virtual void test(UnitTest&);                    //!< Perform unit test
        };

        /*!
         * \brief Harmonic pair potential
         * \details The harmonic potential has the form
         * \f$ \beta u_{ij} = k(r_{ij}-r_{eq})^2 \f$ where k is the force constant
         * (kT/angstrom^2) and req is the equilibrium distance (angstrom).
         * \note We do not multiply with 1/2 which must be included in the supplied force constant, k
         */
        class Harmonic : public PairPotentialBase {
          private:
            string _brief();
            void _setScale(double);
          public:
            double k;   //!< Force constant (kT/A^2) - Did you rember to divide by two? See note.
            double req; //!< Equilibrium distance (angstrom)
            Harmonic(double=0, double=0);
            Harmonic(InputMap&, string="harmonic_");
            inline double operator() (const particle &a, const particle &b, double r2) const FOVERRIDE {
              double d=sqrt(r2)-req;
              return k*d*d;
            }
        };

        /*!
         * \brief Cosine attraction ala Cooke and Deserno
         * \details This is an attractive potential used for coarse grained lipids
         * and has the form:
         * \f[
         *     \beta u(r) = -\epsilon \cos^2 [ \pi(r-r_c)/2w_c ]
         * \f]
         * for \f$r_c\leq r \leq r_c+w_c\f$. For \f$r<r_c\f$, \f$\beta u=-\epsilon\f$, while
         * zero for \f$r>r_c+w_c\f$.
         * The InputMap parameters are:
         * \li \c cosattract_eps  Depth, \f$\epsilon\f$ [kT]
         * \li \c cosattract_rc   Width, r_c [angstrom]
         * \li \c cosattract_wc   Decay range, w_c [angstrom]
         *
         * \warning Untested!
         */
        class CosAttract : public PairPotentialBase {
          private:
            double eps, wc, rc, rc2, c, rcwc2;
            string _brief();
          public:
            CosAttract(InputMap&, string="cosattract_");
            inline double operator() (const particle &a, const particle &b, double r2) const FOVERRIDE {
              if (r2<rc2)
                return -eps;
              if (r2>rcwc2)
                return 0;
              double x=cos( c*( sqrt(r2)-rc ) );
              return -eps*x*x;
            }
            string info(char); // More verbose information
        };

        /*!
         * \brief Finite Extensible nonlinear elastic (FENE) potential
         * \details This is an anharmonic bonding potential with the form:
         * \f[
         *     \beta u(r) = -\frac{k r_0^2}{2}\ln \left [ 1-(r/r_0)^2 \right ]
         * \f]
         * for $r<r_0$, otherwise infinity. The input parameters read by InputMap
         * are as follows:
         * \li \c fene_stiffness Bond stiffness, k [kT]
         * \li \c fene_maxsep Maximum separation, r_0 [angstrom]
         *
         * \warning Untested!
         */
        class FENE : public PairPotentialBase {
          private:
            double k,r02,r02inv;
            string _brief();
          public:
            FENE(double,double); // Constructor
            FENE(InputMap&, string="fene_"); // Constructor
            inline double operator() (const particle &a, const particle &b, double r2) const FOVERRIDE {
              return (r2>r02) ? pc::infty : -0.5*k*r02*std::log(1-r2*r02inv);
            }
        };

        /*!
         * \brief Hard sphere pair potential
         */
        class HardSphere : public PairPotentialBase {
          private:
            string _brief();
          public:
            HardSphere();
            HardSphere(InputMap&);
            inline double operator() (const particle &a, const particle &b, double r2) const FOVERRIDE {
              double mindist=a.radius+b.radius;
              if (r2<mindist*mindist)
                return pc::infty;
              return 0;
            }
            string info(char w);
        };

        /*!
         * \brief Hard pair potential for spherocylinders
         */
        class HardSpheroCylinder : public PairPotentialBase {
          private:
            string _brief();
            Geometry::Geometrybase *geoPtr;
          public:
            struct prop {
              double halfl;
            };

            std::map<particle::Tid, prop> m;

            HardSpheroCylinder(InputMap&);
            inline double operator() (const CigarParticle &p1, const CigarParticle &p2, double r2) {
              Point r_cm = geoPtr->vdist(p1,p2);
              Point distvec = Geometry::mindist_segment2segment(p1.dir, m[p1.id].halfl, p2.dir, m[p2.id].halfl, r_cm );
              double mindist=p1.radius+p2.radius;
              if ( distvec.dot(distvec) < mindist*mindist)
                return pc::infty;
              return 0;		
            }
            string info(char w);
        };

        /*!
         * \brief Lennard-Jones (12-6) pair potential
         * \details The Lennard-Jones potential has the form:
         * \f$
         * \beta u = 4\epsilon_{lj} \left (  (\sigma_{ij}/r_{ij})^{12} - (\sigma_{ij}/r_{ij})^6    \right )
         * \f$
         * where \f$\sigma_{ij} = (\sigma_i+\sigma_j)/2\f$ and \f$\epsilon_{lj}=\epsilon\f$ is the same for all pairs in this class.
         */
        class LennardJones : public PairPotentialBase {
          private:
            string _brief();
            void _setScale(double);
          protected:
            inline double r6(double sigma, double r2) const {
              double x(sigma*sigma/r2);  // 2
              return x*x*x;             // 6
            }
            double eps;
          public:
            LennardJones();
            LennardJones(InputMap&, string="lj_");
            inline double operator() (const particle &a, const particle &b, double r2) const FOVERRIDE {
              double x(r6(a.radius+b.radius,r2));
              return eps*(x*x - x);
            }
            string info(char);
        };

        /*! \brief Lorentz-Berthelot Mixing Rule for Lennard-Jones parameters
        */
        class LorentzBerthelot {
          public:
            string name;
            LorentzBerthelot();
            double mixSigma(double,double) const;
            double mixEpsilon(double,double) const;
        };

        /*!
         * \brief Lennard-Jones potential with arbitrary mixing rules between particle types
         * \details This is a template for Lennard-Jones pair interactions where the template parameter
         * must be a class for the epsilon and sigma mixed rules. The atomic values for 
         * sigma and epsilon are taken from the AtomTypes class via the global instance
         * \c atom. In your InputMap configuration file you would typically set the atom list file using
         * the keyword \c atomlist. Note that sigma for each atom is set to two times the radius found in
         * AtomTypes.
         *
         * For example:
         * \code
         * InputMap mcp("myconfig");
         * LennardJonesMixed<LorentzBerthelot> lj(mcp);
         * \endcode
         * \todo Prettify output
         */
        template<class Tmixingrule>
          class LennardJonesMixed : public PairPotentialBase {
            private:
              Tmixingrule mixer; // mixing rule class for sigma and epsilon
              string _brief() { return name + " w. " + mixer.name; }
            protected:
              vector< vector<double> > s2, eps; // lookup tables for the mixed parameter
            public:
              LennardJonesMixed(InputMap &in) {
                name="Lennard-Jones";
                int n=(int)atom.list.size(); // number of atom types
                eps.resize(n);
                s2.resize(n);
                for (int i=0; i<n; i++) {    // loop over atom types
                  eps[i].resize(n);
                  s2[i].resize(n);
                  for (int j=0; j<n; j++) {
                    s2[i][j] = pow( mixer.mixSigma( atom.list[i].sigma, atom.list[j].sigma), 2);
                    eps[i][j] = 4*mixer.mixEpsilon( atom.list[i].eps, atom.list[j].eps );
                    eps[i][j] = pc::kJ2kT(eps[i][j]); // convert to kT
                  }
                }
              }

              inline double operator() (const particle &a, const particle &b, double r2) const FOVERRIDE {
                double x( s2[a.id][b.id] / r2 ); //s2/r2
                x=x*x*x; // s6/r6
                return eps[a.id][b.id] * (x*x - x);
              }

              /*!
               * \brief This will set a custom epsilon for a pair of particles
               * \param i Particle id of first particle
               * \param j Particle id of second particle
               * \param eps_kT epsilon in units of kT
               */
              void customEpsilon(particle::Tid i, particle::Tid j, double eps_kT) {
                eps[i][j]=4*eps_kT;
                eps[j][i]=eps[i][j];
              }

              void customSigma(particle::Tid i, particle::Tid j, double sigma) {
                s2[i][j]=sigma*sigma;
                s2[j][i]=s2[i][j];
              }

              string info(char w=0) {
                using namespace Faunus::textio;
                std::ostringstream o;
                o << indent(SUB) << name+" pair parameters:\n";
                int n=(int)atom.list.size();
                for (int i=0; i<n; i++)
                  for (int j=0; j<n; j++)
                    if (i>=j)
                      o << indent(SUBSUB) << setw(12) << atom[i].name+"<->"+atom[j].name
                        << indent(SUB) << sigma+" = " << sqrt( s2[i][j] ) << _angstrom
                        << indent(SUB) << epsilon+" = " << eps[i][j]/4 << kT+" = " << pc::kT2kJ(eps[i][j]/4) << " kJ/mol"
                        << endl;
                return o.str();
              }
          };

        /*!
         * \brief Weeks-Chandler-Andersen pair potential
         * \details This is a Lennard-Jones type potential, cut and shifted to zero
         * at \f$r_c=2^{1/6}\sigma\f$. More info can be found in DOI: 10.1063/1.1674820
         * and the functional form is:
         * \f[
         * \beta u = 4 \epsilon \left ( (b/r)^{12} - (b/r)^6 + \frac{1}{4} \right )
         * \f]
         * where sigma, epsilon per default are set from using Lorentz-Berthelot mixing rules.
         *
         * \warning Untested!
         */
        class WeeksChandlerAndersen : public LennardJonesMixed<LorentzBerthelot> {
          protected:
            const double onefourth, twototwosixth;
          public:
            typedef LennardJonesMixed<LorentzBerthelot> Tbase;
            WeeksChandlerAndersen(InputMap&);
            inline double operator() (const particle &a, const particle &b, double r2) const FOVERRIDE {
              double x=s2[a.id][b.id]; // s^2
              if (r2>x*twototwosixth)
                return 0;
              x=x/r2;  // (s/r)^2
              x=x*x*x;// (s/r)^6
              return eps[a.id][b.id]*(x*x - x + onefourth);
            }
            inline double operator() (const particle &a, const particle &b, const Point &r) const {
              return operator()(a,b,r.squaredNorm());
            }
        };
        


        /*
           class cos2 : public PairPotentialBase {
           private:
           string _brief();
           void _setScale(double);
           protected:

           public:
        //TODO cutoff, epsilon, pdist, pswitch!!!
        cos2();
        cos2(InputMap&, string="cos2_");
        inline double operator() (const particle &a, const particle &b, double r2) const FOVERRIDE {
        double atrenergy;
        double ndist=sqrt(r2);
        if (ndist < atom.list[a.id].pdis.....) 
        atrenergy = -interact->param->epsilon;
        else {
        atrenergy = cos(PIH*(ndist-interact->param->pdis)/interact->param->pswitch);
        atrenergy *= -interact->param->epsilon*atrenergy;

        }
        return atrenergy;
        }
        string info(char);
        };  
        */
        //TODO spherocylinder at external wall!!!


        inline double fanglscale(double a, const CigarParticle &p){
          // a = r_ij * n_i
          double f;
          if (a <= p.pcanglsw)
            f=0.0;
          else {
            if (a >= p.pcangl)
              f=1.0;
            else {
              f = 0.5 - ((p.pcanglsw + p.pcangl)*0.5 - a )/(p.pcangl - p.pcanglsw);
            }
          }
          return f;
        }

        template<typename Tcigarsphere >
          class PatchyCigarSphere : public PairPotentialBase {
            private:
              string _brief() {
                return pairpot.brief();
              }
            public:
              Tcigarsphere pairpot;
              vector< vector<double> >* rcutPtr;

              PatchyCigarSphere(InputMap &in) : pairpot(in) {
                rcutPtr=nullptr;
              }

              double operator() (const CigarParticle &a, const CigarParticle &b, const Point &r_cm) {
                assert(rcutPtr!=nullptr);
                //0- isotropic, 1-PSC all-way patch,2 -CPSC cylindrical patch
                //b is sphere, a is spherocylinder
                double s, t, f0, f1, contt;
                double rcut=(*rcutPtr)[a.id][b.id];

                Point vec1,distvec;

                assert( a.halfl < 1e-6 && "First should be cigar then sphere, not opposite!");
                double c = a.dir.dot(r_cm);
                if (c > a.halfl)
                  contt = a.halfl;
                else {
                  if (c > -a.halfl) contt = c;
                  else contt = -a.halfl;
                }
                distvec = -r_cm + (a.dir*contt);

                if (atom.list[a.id].patchtype ==0 ) {
                  if (atom.list[b.id].patchtype == 0) {
                    return pairpot(a,b,distvec.dot(distvec));
                  }
                }

                // scaling function: angular dependence of patch1
                vec1=Geometry::vec_perpproject(distvec, a.dir);
                vec1.normalize();
                s = vec1.dot(a.patchdir);
                f1 = fanglscale(s,a);

                // scaling function for the length of spherocylinder within cutoff
                double ndistsq = distvec.dot(distvec);
                t = sqrt(rcut*rcut- ndistsq);//TODO cutoff
                if ( contt + t > a.halfl )
                  f0 = a.halfl;
                else
                  f0 = contt + t;
                if ( contt - t < -a.halfl )
                  f0 -= -a.halfl;
                else
                  f0 -= contt - t;

                return pairpot(a,b,ndistsq)*f1*(f0+1.0);

              }

              string info(char w) { return pairpot.info(w); }
          };


        template<typename Tcigarcigar>
          class PatchyCigarCigar : public PairPotentialBase {
            private:
              string _brief() {
                return pairpot.brief();
              }
            public:
              vector< vector<double> >* rcutPtr;
              Tcigarcigar pairpot;

              PatchyCigarCigar(InputMap &in) : pairpot(in) {
                rcutPtr=nullptr;
              }

              double operator() (const CigarParticle &a, const CigarParticle &b, const Point &r_cm) {
                assert(rcutPtr!=nullptr);
                //0- isotropic, 1-PSC all-way patch,2 -CPSC cylindrical patch
                if (atom.list[a.id].patchtype >0 ) {
                  if (atom.list[b.id].patchtype > 0) {
                    //patchy sc with patchy sc
                    int i, intrs;
                    double rcut=(*rcutPtr)[a.id][b.id];
                    double ndistsq;
                    double v1, v2, f0, f1, f2, T1, T2, S1, S2,s;
                    double intersections[5];
                    Point vec1, vec2, vec_intrs, vec_mindist;

                    for(i=0;i<5;i++)
                      intersections[i]=0;
                    //1- do intersections of spherocylinder2 with patch of spherocylinder1 at.
                    // cut distance C
                    if (atom.list[a.id].patchtype == 1) {
                      intrs=Geometry::psc_intersect(a,b,r_cm, intersections, rcut);
                    } else {
                      if (atom.list[a.id].patchtype == 2) {
                        intrs=Geometry::cpsc_intersect(a,b,r_cm, intersections, rcut);
                      } else {
                        //we dont have anything like this
                        assert(!"Patchtype not implemented!");
                      }
                    }
                    if (intrs ==0){
                      return 0.0; //sc is all outside patch, attractive energy is 0
                    }
                    T1=intersections[0]; //points on sc2
                    T2=intersections[1];
                    //2- now do the same oposite way psc1 in patch of psc2
                    for(i=0;i<5;i++)
                      intersections[i]=0;
                    if (atom.list[a.id].patchtype == 1) {
                      intrs=Geometry::psc_intersect(b,a,-r_cm, intersections, rcut);
                    } else {
                      if (atom.list[a.id].patchtype == 2) {
                        intrs=Geometry::cpsc_intersect(b,a,-r_cm, intersections, rcut);
                      } else {
                        assert(!"Patchtype not implemented!");
                        //we dont have anything like this
                      }
                    }
                    if (intrs ==0) {
                      return 0.0; //sc is all outside patch, attractive energy is 0
                    }
                    S1=intersections[0]; //points on sc1
                    S2=intersections[1];

                    //3a- with two intersection pices calculate vector between their CM
                    //-this is for angular orientation
                    v1=fabs(S1-S2);
                    v2=fabs(T1-T2);
                    vec1=a.dir*(S1+S2)*0.5;
                    vec2=b.dir*(T1+T2)*0.5;
                    vec_intrs=vec2-vec1-r_cm;
                    //vec_intrs should be from sc1 t sc2

                    //3b - calculate closest distance attractive energy from it
                    vec_mindist = Geometry::mindist_segment2segment(a.dir,v1,b.dir,v2,vec_intrs);
                    ndistsq=vec_mindist.dot(vec_mindist);

                    //4- scaling function1: dependence on the length of intersetions
                    //F0=(V1+V2)*0.5/interact->param->sigma;
                    f0=(v1+v2)*0.5+1.0;
                    //5- scaling function2: angular dependence of patch1
                    vec1=Geometry::vec_perpproject(vec_intrs, a.dir);
                    vec1.normalize();
                    s = vec1.dot(a.patchdir);
                    f1 = fanglscale(s,a);

                    //6- scaling function3: angular dependence of patch2
                    vec1=Geometry::vec_perpproject(-vec_intrs, b.dir);
                    vec1.normalize();
                    s = vec1.dot(b.patchdir);
                    f2 = fanglscale(s,b);

                    //7- put it all together and output scale
                    return f0*f1*f2*pairpot(a,b,ndistsq);
                  }
                  else {
                    assert(!"Patchy cigar w. isotropic cigar not implemented!");
                    //patchy sc with isotropic sc
                    //we dont have at the moment
                  }
                } else {
                  if (atom.list[b.id].patchtype > 0) {
                    assert(!"Patchy cigar w. isotropic cigar not implemented!");
                    //isotropic sc with patchy sc
                    //we dont have at the moment
                  }
                  else {
                    //isotropic sc with isotropic sc
                    Point rclose=Geometry::mindist_segment2segment(a.dir, a.halfl, b.dir, b.halfl, r_cm);
                    return pairpot(a,b,rclose.squaredNorm());
                  }

                }
                //something we have not implemented
                return 0.0;
              }
              string info(char w) { return pairpot.info(w); }
          };

        template<typename Tcigarcigar, typename Tspheresphere, typename Tcigarsphere>
          class CigarSphereSplit : public PairPotentialBase {
            private:
              string _brief() {
                return pairpot_cc.brief() + " "
                  + pairpot_ss.brief() + " "
                  + pairpot_cs.brief();
              }

            public:
              PatchyCigarCigar<Tcigarcigar> pairpot_cc;
              Tspheresphere pairpot_ss;
              PatchyCigarSphere<Tcigarsphere> pairpot_cs;

              CigarSphereSplit(InputMap &in) : pairpot_cc(in), pairpot_ss(in), pairpot_cs(in){
                name="CigarSphereSplit";
              }

              double operator() (const particle &a, const particle &b, double r2) const {
                assert(!"Undefined!");
                return 0;
              }

              double operator() (const CigarParticle &a, const CigarParticle &b, const Point &r_cm)
              {
                if (a.halfl<1e-6) {
                  // a sphere - b sphere
                  if (b.halfl<1e-6) {
                    return pairpot_ss(a,b,r_cm.squaredNorm());
                  }
                  // a sphere - b cigar
                  else {
                    //PatchyCigarSphere(b,a)
                    return pairpot_cs(b,a,r_cm);
                  }
                } else {
                  // a cigar - b sphere
                  if (b.halfl<1e-6) {
                    //PatchyCigarSphere(a,b)
                    return pairpot_cs(a,b,r_cm);
                  }
                  // a cigar - b cigar
                  else {
                    //PatchyCigarCigar
                    return pairpot_cc(a,b,r_cm);
                  }
                }
                return 0;
              }

              string info(char w) {
                return pairpot_cc.info(w) + pairpot_ss.info(w) + pairpot_cs.info(w);
              }
          };

        /*!
         * \brief Square well pair potential
         */
        class SquareWell : public PairPotentialBase {
          private:
            string _brief();
            void _setScale(double);
          public:
            double threshold;                           //!< Threshold between particle *surface* [A]
            double depth;                               //!< Energy depth [kT] (positive number)
            SquareWell(InputMap&, string="squarewell"); //!< Constructor
            inline double operator() (const particle &a, const particle &b, double r2) const FOVERRIDE {
              double d=a.radius+b.radius+threshold;
              if ( r2 < d*d )
                return -depth;
              return 0;
            }
            string info(char);
        };

        /*!
         * \brief Square well pair potential shifted
         * \author Anil Kurut
         */
        class SquareWellShifted : public SquareWell {
          private:
            string _brief();
          public:
            double threshold_lower;
            SquareWellShifted(InputMap&, string="squarewell"); //!< Constructor
            inline double operator() (const particle &a, const particle &b, double r2) const FOVERRIDE {
              double d=a.radius+b.radius+threshold_lower;
              if ( r2 > d*d )
                return SquareWell::operator()(a,b,r2);
              return 0;
            }
            string info(char);
        };

        /*!
         * \brief Hydrophobic pair potential based on SASA and surface tension
         * \todo Documentation is incorrect.
         * \details The potential is not zero if the distance between hydrophobic particles
         * is smaller than size of solvent molecule (2*Rs)  
         * Potential has the form:
         *
         * \f$ u = Surface tension * (\Delta SASA_i + \Delta SASA_j) \f$
         *
         * Surface area which is not accesible for solvent \f$ \Delta SASA_i = (SASA_i(r_{ij})-SASA_i(\inf)) \f$ is calculated based on surface of a sphere cap
         *
         * \f$ SA_{cap,i}=2\pi(R_i+R_s)h_i \f$ where h is dependent on distance between the particles as 
         *
         * \f$ h_i=(R_i+R_s)*(\frac{(R_i+R_s)^2-(R_j+R_s)^2)+r_{ij}^2}{2r_{ij}(R_i+R_s)}\f$
         *
         */
        class SquareWellHydrophobic : public SquareWell {
          public:
            SquareWellHydrophobic(InputMap&, string="squarewell"); //!< Constructor
            inline double operator() (const particle &a, const particle &b, double r2) const FOVERRIDE {
              if (a.hydrophobic)
                if (b.hydrophobic)
                  return SquareWell::operator()(a,b,r2);
              return 0;
            }
        };

        /*!
         * \brief Soft repulsion of the form \f$ \beta u = \sigma^6 / (r_{ij}-r_i-r_j)^6 \f$
         * \todo This applies sqrt() and thus may be slow. Also remove floating point comparison.
         */
        class SoftRepulsion : public PairPotentialBase {
          private:
            string _brief();
            void _setScale(double);
            double sigma6;
          public:
            SoftRepulsion(InputMap&);
            inline double operator() (const particle &a, const particle &b, double r2) const FOVERRIDE {
              double d=a.radius+b.radius;
              if (r2<=d*d) //fp comparison!
                return pc::infty;
              d = sqrt(r2)-d;
              d=d*d; //d^2
              return sigma6 / (d*d*d); // s^6/d^6
            }
            string info(char);
        };

        /*!
         * \brief r12-Repulsion of the form
         * \details \f$ \beta u = 4\epsilon_{lj} \left (  (\sigma_{ij}/r_{ij})^{12}  \right ) \f$
         * where \f$\sigma_{ij} = (\sigma_i+\sigma_j)/2\f$ and \f$\epsilon_{lj}\f$ is fixed for this class.
         * \todo Same as LennardJonesR12. Remove?
         */
        class R12Repulsion : public PairPotentialBase {
          private:
            string _brief();
            void _setScale(double);
          protected:
            double eps;
          public:
            R12Repulsion(); // __attribute__ ((deprecated));
            R12Repulsion(InputMap&, string="r12rep_"); // __attribute__ ((deprecated));
            inline double operator() (const particle &a, const particle &b, double r2) const FOVERRIDE {
              double x=(a.radius+b.radius);
              x=x*x/r2; // r2
              x=x*x*x; // r6
              return eps*x*x;
            }
            string info(char);
        };
    
    
    

        /*!
         * \brief As LennardJones but only the repulsive R12 part.
         */
        class LennardJonesR12 : public LennardJones {
          public:
            LennardJonesR12(InputMap&, string="r12rep_");
            inline double operator() (const particle &a, const particle &b, double r2) const FOVERRIDE {
              double x=r6(a.radius+b.radius,r2);
              return eps*x*x;
            }
        };
    
        /*!
         * \brief As LennardJones but truncated and shifted to sigma.
         */
        class LennardJonesTrunkShift : public LennardJones {
        public:
          LennardJonesTrunkShift(InputMap&, string="ljts_");
          inline double operator() (const particle &a, const particle &b, double r2) const FOVERRIDE {
            double sigma = a.radius+b.radius;
            if (r2 > sigma*sigma)
              return 0.0;
            
            double x=r6(sigma,r2)*0.5;
            return eps*(x*x - x + 0.25);
          }
        };


        /*!
         * \brief Coulomb pair potential between charges in a dielectric medium.
         * \details The Coulomb potential has the form
         * \f[
         * \beta u_{ij} = \frac{e^2}{4\pi\epsilon_0\epsilon_rk_BT} \frac{z_i z_j}{r_{ij}}
         *              = \lambda_B \frac{z_i z_j}{r_{ij}}
         * \f]
         * where \f$\lambda_B\f$ is the Bjerrum length and \c z are the valencies.
         */
        class Coulomb : public PairPotentialBase {
          friend class Potential::DebyeHuckel;
          friend class Energy::GouyChapman;
          private:
          string _brief();
          void _setScale(double);
          double epsilon_r;
          protected:
          double depsdt;      //!< \f$ T\partial \epsilon_r / \epsilon_r \partial T = -1.37 \f$
          double lB;          //!< Bjerrum length (angstrom)

          public:
          Coulomb(InputMap&); //!< Construction from InputMap
          double bjerrumLength() const;  //!< Returns Bjerrum length [AA]

          inline double operator() (const particle &a, const particle &b, double r2) const FOVERRIDE {
#ifdef FAU_APPROXMATH
            return lB*a.charge*b.charge * invsqrtQuake(r2);
#else
            return lB*a.charge*b.charge / sqrt(r2);
#endif
          }
          
          string info(char);
          void test(UnitTest&); //!< Perform unit test
        };

        /*!
         * \brief Coulomb pair potential shifted according to Wolf/Yonezawaa -- doi:10.1063/1.4729748
         * \details The Coulomb potential has the form:
         * \f[
         * \beta u_{ij} = \frac{e^2}{4\pi\epsilon_0\epsilon_rk_BT}
         * z_i z_j \left (
         * \frac{1}{r} - \frac{1}{R_c} + \frac{r-R_c}{R_c^2}
         * \right )
         * \f]
         */
        class CoulombWolf : public Coulomb {
          private:
            double Rc2, Rcinv;
          public:
            CoulombWolf(InputMap&); //!< Construction from InputMap
            inline double operator() (const particle &a, const particle &b, double r2) const FOVERRIDE {
              if (r2>Rc2)
                return 0;
#ifdef FAU_APPROXMATH
              r2=invsqrtQuake(r2);  // 1/r
              return lB * a.charge * b.charge * (r2 - Rcinv + (Rcinv/r2-1)*Rcinv );
#else
              r2=sqrt(r2); // r
              return lB * a.charge * b.charge * (1/r2 - Rcinv + (r2*Rcinv-1)*Rcinv );
#endif
            }
            string info(char);
        };

        /*!
         * \brief Charge-nonpolar pair interaction
         * \details This accounts for polarization of
         * \f[
         * \beta u_{ij} = -\frac{\lambda_B z_i^2 \delta a_j^3}{2r_{ij}^4}
         * \f]
         * where a is the radius of the nonpolar particle. Note that this version requires that one of the particles
         * is charged, while the other is neutral. Delta is a unitless scaling parameter of the excess
         * polarizability. For non-polar particles in a polar medium, this is a negative number. For more information, see
         * Israelachvili, Chapter 5. The InputMap is scanned for
         * \li The parameters from Potential::Coulomb
         * \li \c excess_polarization for the delta value
         */
        class ChargeNonpolar : public Coulomb {
          private:
            double c;
          public:
            ChargeNonpolar(InputMap&); //!< Construction from InputMap
            inline double operator() (const particle &a, const particle &b, double r2) const FOVERRIDE {
              if ( abs(a.charge)>1e-6 )
                return -c * a.charge * a.charge / (r2*r2) * (b.radius*b.radius*b.radius);
              else if ( abs(b.charge)>1e-6 )
                return -c * b.charge * b.charge / (r2*r2) * (a.radius*a.radius*a.radius);
              return 0;
            }
            string info(char);
        };

        /*!
         * \brief Debye-Huckel/Yukawa pair potential
         * \details This potential is similar to the plain Coulomb potential but with an extra exponential term to described salt screening:
         * \f[ \beta w_{ij} = \frac{e^2}{4\pi\epsilon_0\epsilon_rk_BT} \frac{z_i z_j}{r_{ij}} \exp(-\kappa r_{ij}) \f]
         * where \f$\kappa=1/D\f$ is the inverse Debye screening length.
         */
        class DebyeHuckel : public Coulomb {
          private:
            string _brief();
          protected:
            double c,k;
          public:
            DebyeHuckel(InputMap&);                       //!< Construction from InputMap
            inline double operator() (const particle &a, const particle &b, double r2) const FOVERRIDE {
#ifdef FAU_APPROXMATH
              double rinv = invsqrtQuake(r2);
              return lB * a.charge * b.charge * rinv * exp_cawley(-k/rinv);
#else
              double r=sqrt(r2);
              return lB * a.charge * b.charge / r * exp(-k*r);
#endif
            }
            double entropy(double, double) const;         //!< Returns the interaction entropy 
            double ionicStrength() const;                 //!< Returns the ionic strength (mol/l)
            double debyeLength() const;                   //!< Returns the Debye screening length (angstrom)
            double excessChemPot(double, double=0) const; //!< Single ion excess chemical potential (kT)
            double activityCoeff(double, double=0) const; //!< Single ion activity coefficient (molar scale) 
            string info(char);
        };

        /*!
         * \brief as DebyeHuckel but shifted to reaach zero at a user specified cut-off distance
         * \details The cut-off distance is read from the InputMap with the following keyword:
         * \li \c pairpot_cutoff Spherical cut-off in angstroms
         */
        class DebyeHuckelShift : public DebyeHuckel {
          private:
            double shift;    // offset at cutoff distance
            double sqcutoff; // squared cutoff distance
          public:
            DebyeHuckelShift(InputMap&);                       //!< Construction from InputMap
            inline double operator() (const particle &a, const particle &b, double r2) const {
              if (r2>sqcutoff)
                return 0;
#ifdef FAU_APPROXMATH
              double rinv = invsqrtQuake(r2);
              return lB * a.charge * b.charge * ( exp_cawley(-k/rinv)*rinv - shift );
#else
              double r=sqrt(r2);
              return lB * a.charge * b.charge * ( exp(-k*r)/r - shift );
#endif
            }
        };

        /*!
         * \brief Combines two pair potentials
         * \details This combines two PairPotentialBases. The combined potential can subsequently
         * be used as a normal pair potential and even be combined with a third potential and
         * so forth. A number of typedefs such as Potential::CoulombHS is aliasing this.
         *
         * \code
         *   // mix two and three pair potentials
         *   using namespace Potential;
         *   typedef CombinedPairPotential< LennardJones, SquareWell > Tpairpot1;
         *   typedef CombinedPairPotential< Tpairpot1, Coulomb > Tpairpot2;
         *   Tpairpot2 mypairpot;
         *   std::cout << mypairpot.info();
         * \endcode
         *
         * \note tokT() functions are *ignored* and the two pair potentials must return
         *       energies in kT directly. The plan is to remove tokT() permanently.
         * \date Lund, 2012
         * \author Mikael Lund
         */
        template<class T1, class T2>
          class CombinedPairPotential : public PairPotentialBase {
            private:
              string _brief() { return first.brief() + " " + second.brief(); }
              //void _setScale(double s) {}
            public:
              T1 first;  //!< First pair potential of type T1
              T2 second; //!< Second pair potential of type T2
              CombinedPairPotential(InputMap &in) : first(in), second(in) {
                name=first.name+"+"+second.name;
              }
              CombinedPairPotential(InputMap &in, string pfx1, string pfx2) : first(in,pfx1), second(in,pfx2) {
                name=first.name+"+"+second.name;
              }
              inline double operator() (const particle &a, const particle &b, double r2) const FOVERRIDE {
                return first(a,b,r2) + second(a,b,r2);
              }
              inline Point force(const particle &a, const particle &b, double r2) FOVERRIDE {
                return first.force(a,b,r2) + second.force(a,b,r2);
              }
              string info(char w=20) { return first.info(w) + second.info(w); }
              void test(UnitTest &t) {
                first.test(t);
                second.test(t);
              }
          };

        /*!
         * \brief As CombinedPairPotential but with a spherical cut-off
         *
         * The following keywords are read from the InputMap:
         * \li \c pairpot_cutoff Spherical cut-off distance (Default: infinity)
         */
        template<class T1, class T2>
          class CombinedPairPotentialCutoff : public CombinedPairPotential<T1,T2> {
            private:
              typedef CombinedPairPotential<T1,T2> Tbase;
              double cutoff2; // squared, spherical cutoff
              void getin(InputMap &in) {
                cutoff2 = pow(in.get<double>("pairpot_cutoff",pc::infty),2);
              }
            public:
              CombinedPairPotentialCutoff(InputMap &in) : Tbase(in) { getin(in); }

              CombinedPairPotentialCutoff(InputMap &in, string pfx1, string pfx2) : Tbase(in,pfx1,pfx2) {
                getin(in);
              }

              inline double operator() (const particle &a, const particle &b, double r2) const {
                return (r2>cutoff2) ? 0 : Tbase::operator()(a,b,r2);
              }

              string info(char w=20) {
                std::ostringstream o;
                o << Tbase::info(w)
                  << textio::pad(textio::SUB,w,"Pair potential cutoff") << sqrt(cutoff2) << textio::_angstrom << endl;
                return o.str();
              }
          };

        class MultipoleEnergy {
          public:
            double lB;
            double ionion(double, double, double);
            double iondip(double, const Point&, double);
            double dipdip(const Point&, const Point&, double);
        };

        /*!
         * \brief Lennard-Jones potential with Lorentz-Berthelot mixing rule
         */
        typedef LennardJonesMixed<LorentzBerthelot> LennardJonesLB;

        /*!
         * \brief Combined Coulomb / HardSphere potential
         */
        typedef CombinedPairPotential<Coulomb, HardSphere> CoulombHS;

        /*!
         * \brief Combined Coulomb / LennardJones potential
         */
        typedef CombinedPairPotential<Coulomb, LennardJones> CoulombLJ;
    
        /*!
         * \brief Combined Coulomb / WeeksChandlerAndersen potential
         */
        typedef CombinedPairPotential<Coulomb, WeeksChandlerAndersen> CoulombWCA;
        
        /*!
         * \brief Combined Coulomb / LennardJonesTrunkShift potential
         */
        typedef CombinedPairPotential<Coulomb, LennardJonesTrunkShift> CoulombLJTS;

        /*!
         * \brief Combined Coulomb / LennardJones potential
         */
        typedef CombinedPairPotential<CoulombWolf, LennardJones> CoulombWolfLJ;

        /*!
         * \brief Combined DebyeHuckel / HardSphere potential
         */
        typedef CombinedPairPotential<DebyeHuckel, HardSphere> DebyeHuckelHS;

        /*!
         * \brief Combined DebyeHuckel / LennardJones potential
         */
        typedef CombinedPairPotential<DebyeHuckel, LennardJones> DebyeHuckelLJ;

        /*!
         * \brief Combined DebyeHuckel / R12Repulsion potential
         */
        typedef CombinedPairPotential<DebyeHuckel, R12Repulsion> DebyeHuckelr12;

    } //end of Potential namespace

  } //end of Faunus namespace
#endif
