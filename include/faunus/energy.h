#ifndef FAUNUS_ENERGY_H
#define FAUNUS_ENERGY_H

#ifndef SWIG
#include <faunus/common.h>
#include <faunus/point.h>
#include <faunus/group.h>
#include <faunus/inputfile.h>
#include <faunus/space.h>
#include <faunus/textio.h>
#include <faunus/potentials.h>
#include <faunus/auxiliary.h>
#include <faunus/analysis.h>
#endif

namespace Faunus {

  /**
   * @brief Classes/templates that calculate the energy of particles, groups, system.
   *
   * This namespace containes classes and templates for calculating
   * energies of the system. This can be pair interactions, external
   * potentials, long range corrections, constraints etc.
   * All energy calculation classes inherit from the base class
   * `Energy::Energybase` and thus have a common interface.
   */
  namespace Energy {

    /**
     *  @brief Base class for energy evaluation
     *
     *  This base class defines functions for evaluating interactions
     *  between particles, groups, external potentials etc.
     *  By default all energy functions return ZERO and derived classes are
     *  expected only to implement functions relevant for certain
     *  properties. I.e. a derived class for non-bonded interactions are
     *  not expected to implement `i_internal()`, for example.
     *
     *  @note All energy functions must return energies in units of kT.
     */
    template<class Tspace>
      class Energybase {
        protected:
          virtual std::string _info()=0;
          char w; //!< Width of info output
          Tspace* spc;
        public:
          typedef Tspace SpaceType;
          typedef typename Tspace::ParticleType Tparticle;
          typedef typename Tspace::GeometryType Tgeometry;
          typedef typename Tspace::p_vec Tpvec;

          string name;  //!< Short informative name

          inline virtual ~Energybase() {}

          inline Energybase() : w(25), spc(nullptr) {}

          inline virtual void setSpace(Tspace &s) { spc=&s; } 

          inline virtual Tspace& getSpace() {
            assert(spc!=nullptr);
            return *spc;
          }

          /** @brief Particle-particle energy */
          virtual double p2p(const Tparticle&, const Tparticle&)
          { return 0; }

          virtual Point f_p2p(const Tparticle&, const Tparticle&) // Particle-particle force
          { return Point(0,0,0); }

          virtual double all2p(const Tpvec&, const Tparticle&)  // Particle-Particle vector energy
          { return 0; }

          virtual double i2i(const Tpvec&, int, int)           // i'th particle with j'th particle
          { return 0; }

          virtual double i2g(const Tpvec&, Group &, int)       // i'th particle with group
          { return 0; }

          virtual double i2all(Tpvec&, int)              // i'th particle with all other particles
          { return 0; }

          virtual double i_external(const Tpvec&, int)         // internal energy of i'th particle
          { return 0; }

          virtual double i_internal(const Tpvec&, int)         // External energy of i'th particle
          { return 0; }

          virtual double p_external(const Tparticle&)   // External energy of particle
          { return 0; }

          /** @brief Total energy of i'th = i2all+i_external+i_internal */
          double i_total(Tpvec &p, int i)
          { return i2all(p,i) + i_external(p,i) + i_internal(p,i); }

          virtual double g2g(const Tpvec&, Group&, Group&)     // Group-Group energy
          { return 0; }

          virtual double g_external(const Tpvec&, Group&)      // External energy of group
          { return 0; }

          virtual double g_internal(const Tpvec&, Group&)      // Internal energy of group
          { return 0; }

          virtual double v2v(const Tpvec&, const Tpvec&)       // Particle vector-Particle vector energy
          { return 0; }

          virtual double external(const Tpvec&)                // External energy - pressure, for example.
          { return 0; }

          virtual void field(const Tpvec&, Eigen::MatrixXd&) //!< Calculate electric field on all particles
          { }

          inline virtual std::string info() {
            assert(!name.empty() && "Energy name cannot be empty");
            if (_info().empty())
              return string();
            return textio::header("Energy: " + name) + _info();
          }
      };

    /**
     * @brief Add two energy classes together
     */
    template<class T1, class T2>
      class CombinedEnergy : public Energybase<typename T1::SpaceType> {
        private:
          string _info() { return first.info()+second.info(); }
          typedef Energybase<typename T1::SpaceType> Tbase;
          typedef typename Tbase::Tparticle Tparticle;
          typedef typename Tbase::Tpvec Tpvec;
        public:
          T1 first;
          T2 second;

          CombinedEnergy(const T1 &a, const T2 &b) : first(a), second(b) {}

          string info() { return _info(); }

          void setSpace(typename T1::SpaceType &s) FOVERRIDE {
            first.setSpace(s);
            second.setSpace(s);
            Tbase::setSpace(s);
          } 

          double p2p(const Tparticle &a, const Tparticle &b) FOVERRIDE
          { return first.p2p(a,b)+second.p2p(a,b); }

          Point f_p2p(const Tparticle &a, const Tparticle &b) FOVERRIDE
          { return first.f_p2p(a,b)+second.f_p2p(a,b); }

          double all2p(const Tpvec &p, const Tparticle &a) FOVERRIDE
          { return first.all2p(p,a)+second.all2p(p,a); }

          double i2i(const Tpvec &p, int i, int j) FOVERRIDE
          { return first.i2i(p,i,j)+second.i2i(p,i,j); }

          double i2g(const Tpvec &p, Group &g, int i) FOVERRIDE
          { return first.i2g(p,g,i)+second.i2g(p,g,i); }

          double i2all(Tpvec &p , int i) FOVERRIDE
          { return first.i2all(p,i)+second.i2all(p,i); }

          double i_external(const Tpvec&p, int i) FOVERRIDE
          { return first.i_external(p,i)+second.i_external(p,i); }

          double i_internal(const Tpvec&p, int i) FOVERRIDE
          { return first.i_internal(p,i)+second.i_internal(p,i); }

          double g2g(const Tpvec&p, Group&g1, Group&g2) FOVERRIDE
          { return first.g2g(p,g1,g2)+second.g2g(p,g1,g2); }

          double g_external(const Tpvec&p, Group&g) FOVERRIDE
          { return first.g_external(p,g)+second.g_external(p,g); }

          double g_internal(const Tpvec&p, Group&g) FOVERRIDE
          { return first.g_internal(p,g)+second.g_internal(p,g); }

          double external(const Tpvec&p) FOVERRIDE
          { return first.external(p)+second.external(p); }

          double v2v(const Tpvec&p1, const Tpvec&p2) FOVERRIDE
          { return first.v2v(p1,p2)+second.v2v(p1,p2); }

          void field(const Tpvec&p, Eigen::MatrixXd&E) FOVERRIDE
          { first.field(p,E); second.field(p,E); }
      };

    /**
     * @brief Operator to conveniently add two energy classes together
     */
    template<class T1, class T2,
      class = typename
        std::enable_if<std::is_base_of<Energybase<typename T1::SpaceType>,T1>::value>::type,
      class = typename
        std::enable_if<std::is_base_of<Energybase<typename T1::SpaceType>,T2>::value>::type>
        CombinedEnergy<T1,T2>& operator+(const T1 &u1, const T2 &u2)
        { return *(new CombinedEnergy<T1,T2>(u1,u2)); }

    template<class Tgeometry> struct FunctorScalarDist {
      template<class Tparticle>
        inline double operator()(const Tgeometry &geo, const Tparticle &a, const Tparticle &b) const {
          return geo.sqdist(a,b);
        }
    };

    template<class Tgeometry> struct FunctorVectorDist {
      template<class Tparticle>
        Point operator()(const Tgeometry &geo, const Tparticle &a, const Tparticle &b) const {
          return geo.vdist(a,b);
        }
    };

    /**
     * @brief Energy class for non-bonded interactions.
     *
     * `Tpairpot` is expected to be a pair potential with the following
     * properties:
     * 
     * - `Tpairpot(const InputMap&)`
     * - `double Tpairpot::operator()(const particle&, const particle&, double sqdist))`
     *
     * For a list of implemented potentials, see the `Faunus::Potential`
     * namespace.
     */
    template<class Tspace, class Tpairpot>
      class Nonbonded : public Energybase<Tspace> {
        protected:
          string _info() { return pairpot.info(25); }
          typedef Energybase<Tspace> Tbase;
          typedef typename Tbase::Tparticle Tparticle;
          typedef typename Tbase::Tpvec Tpvec;

        public:
          typename Tspace::GeometryType geo;
          Tpairpot pairpot;

          Nonbonded(InputMap &in) : geo(in), pairpot(in) {
            static_assert(
                std::is_base_of<Potential::PairPotentialBase,Tpairpot>::value,
                "Tpairpot must be a pair potential" );
            Tbase::name="Nonbonded N" + textio::squared + " - " + pairpot.name;
          }

          void setSpace(Tspace &s) FOVERRIDE {
            geo=s.geo;
            Tbase::setSpace(s);
          } 

          //!< Particle-particle energy (kT)
          inline double p2p(const Tparticle &a, const Tparticle &b) {
            return pairpot( a,b,geo.sqdist(a,b));
          }

          //!< Particle-particle force (kT/Angstrom)
          inline Point f_p2p(const Tparticle &a, const Tparticle &b) FOVERRIDE {
            auto r=geo.vdist(a,b);
            return pairpot.force(a,b,r.squaredNorm(),r);
          }

          double all2p(const Tpvec &p, const Tparticle &a) {
            double u=0;
            for (auto &b : p)
              u+=pairpot(a,b,geo.sqdist(a,b));
            return u;
          }

          double i2i(const Tpvec &p, int i, int j) {
            return pairpot( p[i], p[j], geo.sqdist( p[i], p[j]) );
          }

          double i2g(const Tpvec &p, Group &g, int j) FOVERRIDE {
            double u=0;
            if ( !g.empty() ) {
              int len=g.back()+1;
              if ( g.find(j) ) {   //j is inside g - avoid self interaction
                for (int i=g.front(); i<j; i++)
                  u+=pairpot(p[i],p[j],geo.sqdist(p[i],p[j]));
                for (int i=j+1; i<len; i++)
                  u+=pairpot(p[i],p[j],geo.sqdist(p[i],p[j]));
              } else              //simple - j not in g
                for (int i=g.front(); i<len; i++)
                  u+=pairpot( p[i], p[j], geo.sqdist(p[i],p[j]));
            }
            return u;  
          }

          double i2all(Tpvec &p, int i) FOVERRIDE {
            assert(i>=0 && i<int(p.size()) && "index i outside particle vector");
            double u=0;
            int n=(int)p.size();
            for (int j=0; j!=i; ++j)
              u+=pairpot( p[i], p[j], geo.sqdist(p[i],p[j]) );
            for (int j=i+1; j<n; ++j)
              u+=pairpot( p[i], p[j], geo.sqdist(p[i],p[j]) );
            return u;
          }

          double g2g(const Tpvec &p, Group &g1, Group &g2) FOVERRIDE {
            double u=0;
            if (!g1.empty())
              if (!g2.empty()) {
                // IN CASE ONE GROUP IS A SUBGROUP OF THE OTHER
                if (g1.find(g2.front()))
                  if (g1.find(g2.back())) {  // g2 is a subgroup of g1
                    assert(g1.size()>=g2.size());
                    for (int i=g1.front(); i<g2.front(); i++)
                      for (auto j : g2)
                        u+=pairpot(p[i],p[j],geo.sqdist(p[i],p[j]));
                    for (int i=g2.back()+1; i<=g1.back(); i++)
                      for (auto j : g2)
                        u+=pairpot(p[i],p[j],geo.sqdist(p[i],p[j]));
                    return u;
                  }
                if (g2.find(g1.front()))
                  if (g2.find(g1.back())) {  // g1 is a subgroup of g2
                    assert(g2.size()>=g1.size());
                    for (int i=g2.front(); i<g1.front(); i++)
                      for (auto j : g1)
                        u+=pairpot(p[i],p[j],geo.sqdist(p[i],p[j]));
                    for (int i=g1.back()+1; i<=g2.back(); i++)
                      for (auto j : g1)
                        u+=pairpot(p[i],p[j],geo.sqdist(p[i],p[j]));
                    return u;
                  }

                // IN CASE BOTH GROUPS ARE INDEPENDENT (DEFAULT)
                int ilen=g1.back()+1, jlen=g2.back()+1;
#pragma omp parallel for reduction (+:u)
                for (int i=g1.front(); i<ilen; ++i)
                  for (int j=g2.front(); j<jlen; ++j)
                    u+=pairpot(p[i],p[j],geo.sqdist(p[i],p[j]));
              }
            return u;
          }

          double g_internal(const Tpvec &p, Group &g) FOVERRIDE { 
            double u=0;
            int b=g.back(), f=g.front();
            if (!g.empty())
              for (int i=f; i<b; ++i)
                for (int j=i+1; j<=b; ++j)
                  u+=pairpot(p[i],p[j],geo.sqdist(p[i],p[j]));
            return u;
          }

          double v2v(const Tpvec &p1, const Tpvec &p2) {
            double u=0;
            for (auto &i : p1)
              for (auto &j : p2)
                u+=p2p(i,j);
            return u;
          }
      };

    /**
     * @brief Energy class for non-bonded interactions.
     *
     * `Tpairpot` is expected to be a pair potential with the following properties:
     *
     * - `Tpairpot(const InputMap&)`
     * - `double Tpairpot::operator()(const particle&, const particle&, double sqdist))`
     *
     * For a list of implemented potentials, see the Faunus::Potential namespace.
     */
    template<class Tspace, class Tpairpot>
      class NonbondedVector : public Energybase<Tspace> {
        protected:
          string _info() { return pairpot.info(25); }
          typedef Energybase<Tspace> Tbase;
          typedef typename Tbase::Tparticle Tparticle;
          typedef typename Tbase::Tpvec Tpvec;

        public:
          typename Tspace::GeometryType geo;
          Tpairpot pairpot;

          NonbondedVector(InputMap &in) : geo(in), pairpot(in) {
            static_assert(
                std::is_base_of<Potential::PairPotentialBase,Tpairpot>::value,
                "Tpairpot must be a pair potential" );
            Tbase::name="Nonbonded N" + textio::squared + " - " + pairpot.name;
          }

          void setSpace(Tspace &s) FOVERRIDE {
            geo=s.geo;
            Tbase::setSpace(s);
          } 

          //!< Particle-particle energy (kT)
          inline double p2p(const Tparticle &a, const Tparticle &b) {
            return pairpot( a,b,geo.vdist(a,b));
          }

          //!< Particle-particle force (kT/Angstrom)
          inline Point f_p2p(const Tparticle &a, const Tparticle &b) FOVERRIDE {
            return pairpot.force( a,b,geo.sqdist(a,b),geo.vdist(a,b));
          }

          double all2p(const Tpvec &p, const Tparticle &a) {
            double u=0;
            for (auto &b : p)
              u+=pairpot(a,b,geo.vdist(a,b));
            return u;
          }

          double i2i(const Tpvec &p, int i, int j) {
            return pairpot( p[i], p[j], geo.vdist( p[i], p[j]) );
          }

          double i2g(const Tpvec &p, Group &g, int j) {
            double u=0;
            if ( !g.empty() ) {
              int len=g.back()+1;
              if ( g.find(j) ) {   //j is inside g - avoid self interaction
                for (int i=g.front(); i<j; i++)
                  u+=pairpot(p[i],p[j],geo.vdist(p[i],p[j]));
                for (int i=j+1; i<len; i++)
                  u+=pairpot(p[i],p[j],geo.vdist(p[i],p[j]));
              } else              //simple - j not in g
                for (int i=g.front(); i<len; i++)
                  u+=pairpot( p[i], p[j], geo.vdist(p[i],p[j]));
            }
            return u;  
          }

          double i2all(Tpvec &p, int i) {
            assert(i>=0 && i<int(p.size()) && "index i outside particle vector");
            double u=0;
            for (int j=0; j!=i; ++j)
              u+=pairpot( p[i], p[j], geo.vdist(p[i],p[j]) );
            int n=(int)p.size();
            for (int j=i+1; j<n; ++j)
              u+=pairpot( p[i], p[j], geo.vdist(p[i],p[j]) );
            return u;
          }

          double g2g(const Tpvec &p, Group &g1, Group &g2) FOVERRIDE {
            double u=0;
            if (!g1.empty())
              if (!g2.empty()) {
                // IN CASE ONE GROUP IS A SUBGROUP OF THE OTHER
                if (g1.find(g2.front()))
                  if (g1.find(g2.back())) {  // g2 is a subgroup of g1
                    assert(g1.size()>=g2.size());
                    for (int i=g1.front(); i<g2.front(); i++)
                      for (auto j : g2)
                        u+=pairpot(p[i],p[j],geo.vdist(p[i],p[j]));
                    for (int i=g2.back()+1; i<=g1.back(); i++)
                      for (auto j : g2)
                        u+=pairpot(p[i],p[j],geo.vdist(p[i],p[j]));
                    return u;
                  }
                if (g2.find(g1.front()))
                  if (g2.find(g1.back())) {  // g1 is a subgroup of g2
                    assert(g2.size()>=g1.size());
                    for (int i=g2.front(); i<g1.front(); i++)
                      for (auto j : g1)
                        u+=pairpot(p[i],p[j],geo.vdist(p[i],p[j]));
                    for (int i=g1.back()+1; i<=g2.back(); i++)
                      for (auto j : g1)
                        u+=pairpot(p[i],p[j],geo.vdist(p[i],p[j]));
                    return u;
                  }

                // IN CASE BOTH GROUPS ARE INDEPENDENT (DEFAULT)
                int ilen=g1.back()+1, jlen=g2.back()+1;
#pragma omp parallel for reduction (+:u) schedule (dynamic)
                for (int i=g1.front(); i<ilen; ++i)
                  for (int j=g2.front(); j<jlen; ++j)
                    u+=pairpot(p[i],p[j],geo.vdist(p[i],p[j]));
              }
            return u;
          }

          double g_internal(const Tpvec &p, Group &g) { 
            double u=0;
            const int b=g.back(), f=g.front();
            if (!g.empty())
              for (int i=f; i<b; ++i)
                for (int j=i+1; j<=b; ++j)
                  u+=pairpot(p[i],p[j],geo.vdist(p[i],p[j]));
            return u;
          }

          double v2v(const Tpvec &p1, const Tpvec &p2) {
            double u=0;
            for (auto &i : p1)
              for (auto &j : p2)
                u+=p2p(i,j);
            return u;
          }

          /**
           * Calculates the electric field on all particles
           * and stores (add) in the vector `E`.
           *
           * @param p Particle vector
           * @param E Holds field on each particle. Must have N columns.
           */
          void field(const Tpvec &p, Eigen::MatrixXd &E) FOVERRIDE {
            assert((int)p.size()==E.cols());
            size_t i=0;
            for (auto &pi : p) {
              for (auto &pj : p)
                if (&pi!=&pj)
                  E.col(i) = E.col(i) + pairpot.field(pj, geo.vdist(pi,pj));
              i++;
            }
          }
      };

    /**
     * @brief Cuts group-to-group interactions at specified mass-center separation
     *
     * For two molecular groups (`Group::isMolecular()==true`) this will invoke
     * a mass center cut-off and return zero if beyond. This is best used
     * in connection with cut-off based pair potentials as a simple alternative
     * to cell lists. In this case the `g2g` cutoff should be set to the
     * particle-particle cutoff plus the maximum radii of the two groups,
     * @f$ r_{cut}^{g2g} = r_{cut}^{p2p} + a_{g1} + a_{g2} @f$.
     *
     * Upon construction the `InputMap` is searched for the keyword
     * `g2g_cutoff` - the default value is infinity.
     */
    template<class Tspace, class Tpairpot>
      class NonbondedCutg2g : public Nonbonded<Tspace,Tpairpot> {
        private:
          typedef Nonbonded<Tspace,Tpairpot> base;
          double rcut2;

          Point getMassCenter(const typename base::Tpvec &p, const Group &g) {
            assert(&p==&base::spc->p || &p==&base::spc->trial);
            return (&p==&base::spc->p) ? g.cm : g.cm_trial;
          }

          bool cut(const typename base::Tpvec &p, const Group &g1, const Group &g2) {
            if (g1.isMolecular())
              if (g2.isMolecular()) {
                Point a = getMassCenter(p,g1);
                Point b = getMassCenter(p,g2);
                if (base::geo.sqdist(a,b)>rcut2)
                  return true;
              }
            return false;
          }
          
        public:
          bool noPairPotentialCutoff; //!< Set if range of pairpot is longer than rcut (default: false)

          NonbondedCutg2g(InputMap &in) : base(in) {
            noPairPotentialCutoff=false;
            rcut2 = pow(in.get<double>("g2g_cutoff",pc::infty), 2);
            base::name+=" (g2g cut=" + std::to_string(sqrt(rcut2))
              + textio::_angstrom + ")";
          }

          double g2g(const typename base::Tpvec &p, Group &g1, Group &g2) FOVERRIDE {
            return cut(p,g1,g2) ? 0 : base::g2g(p,g1,g2);
          }

          double i2g(const typename base::Tpvec &p, Group &g, int i) FOVERRIDE {
            auto gi=base::spc->findGroup(i);
            return cut(p,g,*gi) ? 0 : base::i2g(p,g,i);
          }

          /**
           * If no pair potential cutoff is applied, `i2all` may need
           * to re-calculate the full g2g interaction as the mass centers
           * may have moved. To activate this behavior set
           * `noPairPotentialCutoff=true`.
           */
          double i2all(typename base::Tpvec &p, int i) FOVERRIDE {
            if (noPairPotentialCutoff) {
              auto gi=base::spc->findGroup(i);
              assert(gi!=nullptr);
              double u=base::i2g(p,*gi,i); // i<->own group
              for (auto gj : base::spc->groupList())
                if (gj!=gi)
                  u+=g2g(p,*gi,*gj);
              return u;
            }
            else
            {
              double u=0;
              auto gi=base::spc->findGroup(i);
              for (auto gj : base::spc->groupList())
                if (!cut(p,*gi,*gj))
                  u+=base::i2g(p,*gj,i);
              return u;
            }
          }
          
          // unfinished
          void field(const typename base::Tpvec &p, Eigen::MatrixXd &E) FOVERRIDE {
            assert((int)p.size()==E.cols());
            assert(!"Validate this!");

            // double loop over all groups and test cutoff
            // todo: openmp pragma
            for (auto gi : base::spc->groupList())
              for (auto gj : base::spc->groupList())
                if (gi!=gj)
                  if (!cut(p,*gi,*gj))
                    for (int i : *gi)
                      for (int j : *gj)
                        E.col(i)+=base::pairpot.field(p[j],base::geo.vdist(p[i],p[j]));
              
              // now loop over all internal particles in groups
              for (auto g : base::spc->groupList())
                for (int i : *g)
                  for (int j : *g)
                    if (i!=j)
                      E.col(i)+=base::pairpot.field(p[j],base::geo.vdist(p[i],p[j]));
          }
      };

    /**
      @brief Treats groups as charged monopoles beyond cut-off
      */
    template<class Tspace, class Tpairpot, class Tmppot=Potential::DebyeHuckel>
      class NonbondedCutg2gMonopole : public NonbondedCutg2g<Tspace,Tpairpot> {
        private:
          double qscale,Q,R; // charge scaling
          Tmppot dh;
          typedef NonbondedCutg2g<Tspace,Tpairpot> base;
          string _info() {
            using namespace Faunus::textio;
            std::ostringstream o;
            o << pad(SUB,30,"Monopole radius") << R << _angstrom << endl
              << pad(SUB,30,"Monopole charge") << Q/qscale << endl
              << pad(SUB,30,"DH charge scaling") << qscale << endl;
            return o.str();
          }
        public:
          NonbondedCutg2gMonopole(InputMap &in) : base(in), dh(in) {
            base::name+="+MP";
            R=in.get<double>("monopole_radius", 0, "g2g cutoff monopole radius (angstrom)");
            double k=1/dh.debyeLength();
            qscale = std::sinh(k*R)/(k*R);
            Q = qscale * in.get<double>("monopole_charge", 0, "g2g cutoff monopole charge number (e)");
          }
          double g2g(const typename base::Tpvec &p, Group &g1, Group &g2) FOVERRIDE {
            if (g1.isMolecular())
              if (g2.isMolecular()) {
                Point a = base::getMassCenter(p,g1);
                Point b = base::getMassCenter(p,g2);
                double r2=base::geo.sqdist(a,b);
                if (r2>base::rcut2) {
                  typename base::Tparticle p1,p2;
                  p1.charge = Q; //qscale*g1.charge(p);
                  p2.charge = Q; //qscale*g2.charge(p);
                  return dh(p1,p2,r2);
                }
              }
            return base::base::g2g(p,g1,g2);
          }
      };

    /**
     * @brief Class for handling bond pairs
     *
     * Takes care of bonded interactions and can handle mixed bond types. If you create bond BETWEEN
     * groups, make sure to set the `CrossGroupBonds` to `true`.
     *
     * Example:
     *
     *     vector<particle> p(...);            // particle vector
     *     int i=10, j=11;                     // particle index
     *     Energy::Bonded b;
     *     b.add(i, j, Potential::Harmonic(0.1,5.0) );
     *     std::cout << b.info();
     *     double rij2 = ... ;                 // squared distance between i and j
     *     double u = b(i,j)( p[i], p[j], rij2 ); // i j bond energy in kT
     *
     * @date Lund, 2011-2012
     */
    template<class Tspace>
      class Bonded : public Energybase<Tspace>,
      protected pair_list<std::function<
                     double(const typename Tspace::ParticleType&,const typename Tspace::ParticleType&,double)> >
    {
      private:
        typedef typename Energybase<Tspace>::Tparticle Tparticle;
        typedef typename Energybase<Tspace>::Tpvec Tpvec;
        typedef pair_list<std::function<double(const Tparticle&,const Tparticle&,double)> > Tbase;
        typedef std::function<Point(const Tparticle&,const Tparticle&,double,const Point&)> Tforce;

        using Energybase<Tspace>::spc;
        std::map<typename Tbase::Tpair, Tforce> force_list;

        string _infolist;
        string _info() {
          using namespace Faunus::textio;
          std::ostringstream o;
          o << pad(SUB,30,"Look for group-group bonds:")
            << (CrossGroupBonds ? "yes (slow)" : "no (faster)") << endl << endl
            << indent(SUBSUB) << std::left
            << setw(7) << "i" << setw(7) << "j" << endl;
          return o.str() + _infolist;
        }

        template<class Tpairpot>
          struct ForceFunctionObject {
            Tpairpot pot;
            ForceFunctionObject(const Tpairpot &p) : pot(p) {}
            Point operator()(const Tparticle &a, const Tparticle &b, double r2, const Point &r) {
              return pot.force(a,b,r2,r);
            }
          };

      public:
        bool CrossGroupBonds; //!< Set to true if bonds cross groups (slower!). Default: false

        Bonded() {
          this->name="Bonded particles";
          CrossGroupBonds=false;
        }

        /**
         * @brief Bond energy i with j
         * @todo Optimize by using iterator directly as returned by find.
         */
        double i2i(const Tpvec &p, int i, int j) FOVERRIDE {
          assert(i!=j);
          auto f=this->list.find( opair<int>(i,j) );
          if (f!=this->list.end())
            return f->second( p[i], p[j], spc->geo.sqdist( p[i], p[j] ) );
          return 0;
        }

        /**
         * @note This will work only for particles contained inside
         * Space main particle vector.
         */
        Point f_p2p(const Tparticle &a, const Tparticle &b) FOVERRIDE {
          int i=&a-&spc->p[0]; // calc. index from addresses
          int j=&b-&spc->p[0];
          assert(i>=0 && j>=0);
          assert(i<(int)spc->p.size() && j<(int)spc->p.size());
          auto f=force_list.find( opair<int>(i,j) );
          if (f!=this->force_list.end()) {
            auto r = spc->geo.vdist(a,b);
            return f->second(a,b,r.squaredNorm(),r);
          }
          return Point(0,0,0);
        }

        //!< All bonds w. i'th particle 
        double i2all(Tpvec &p, int i) FOVERRIDE {
          assert( i>=0 && i<(int)p.size() ); //debug
          double u=0;
          auto eqr=this->mlist.equal_range(i);
          for (auto it=eqr.first; it!=eqr.second; ++it) {
            int j = it->second; // partner index
            u += this->list[opair<int>(i,j)](
                p[i], p[j], spc->geo.sqdist( p[i], p[j] ) );
          }
          return u;
        }

        double total(const Tpvec &p) {
          double u=0;
          for (auto &m : Tbase::list) {
            int i=m.first.first;
            int j=m.first.second;
            assert(i>=0 && i<(int)p.size() && j>=0 && j<(int)p.size()); //debug
            u += m.second( p[i], p[j], spc->geo.sqdist( p[i], p[j] ) );
          }
          return u;
        }

        /**
         * Group-to-group bonds are disabled by default as these are
         * rarely used and the current implementation is rather slow
         * for systems with *many* groups (but very general). To
         * activate `g2g()`, set `CrossGroupBonds=true`.
         *
         * @warning Untested!
         * @todo Possible optimization: swap groups so that `g1<g2`;
         */
        double g2g(const Tpvec &p, Group &g1, Group &g2) FOVERRIDE {
          double u=0;
          if (CrossGroupBonds)
            for (auto i : g1) {
              auto eqr=this->mlist.equal_range(i);
              for (auto it=eqr.first; it!=eqr.second; ++it) {
                int j = it->second; // partner index
                if (g2.find(j))
                  u+=this->list[opair<int>(i,j)](
                      p[i], p[j], spc->geo.sqdist( p[i], p[j] ) );
              }
            }
          return u;
        }

        /**
         * @brief Internal bonds in Group, only
         *
         * Accounts for bonds between particles within a group.
         * Bonds with particles
         * outside the group are skipped and should be accounted for
         * by the g2g() energy function.
         */
        double g_internal(const Tpvec &p, Group &g) FOVERRIDE {
          double u=0;
          if ( this->list.size() > pow(g.size(),2) ) {
            // small group compared to bond list
            auto end=this->list.end();
            for (auto i=g.front(); i<g.back(); i++)
              for (auto j=i+1; j<=g.back(); j++) {
                auto it = this->list.find(opair<int>(i,j));
                if (it!=end)
                  u+=it->second( p[i],p[j],spc->geo.sqdist( p[i],p[j] ) );
              }
          } else {
            // big group compared to bond list
            for (auto &m : this->list) {
              int i=m.first.first, j=m.first.second;
              if (g.find(i))
                if (g.find(j))
                  u += m.second(p[i],p[j],spc->geo.sqdist( p[i],p[j] ));
            }
          }
          return u;
        }

        template<class Tpairpot>
          void add(int i, int j, Tpairpot pot) {
            std::ostringstream o;
            o << textio::indent(textio::SUBSUB) << std::left << setw(7) << i
              << setw(7) << j << pot.brief() + "\n";
            _infolist += o.str();
            pot.name.clear();   // potentially save a
            pot.prefix.clear(); // little bit of memory
            Tbase::add(i,j,pot);// create and add functor to pair list
            force_list[ typename Tbase::Tpair(i,j) ]
              = ForceFunctionObject<decltype(pot)>(pot);
          }
    };

    /**
     * @brief Dummy energy class that sums missed energy changes to avoid energy drifts
     *
     * This energy function is designed to be used with Move::Movebase classes
     * that returns energy changes not detectable in the energy drift checkup
     * routines. The idea is simply to sum the energy change discrepancy
     * and treat
     * this as an external potential.
     */
    template<class Tspace>
      class EnergyRest : public Energybase<Tspace> {
        private:
          double usum;
          string _info(){ return ""; }
        public:
          EnergyRest() {
            usum=0;
            this->name="Dummy energy (drift calc. aid)";
          }

          //!< Add energy change disrepancy, dU = U(metropolis) - U(as in drift calculation)
          void add(double du) { usum+=du; }

          //!< Dumme rest treated as external potential to whole system
          double external(const typename Energybase<Tspace>::Tpvec &p) FOVERRIDE {
            return usum;
          }
      };

    /**
     * @brief Energy from external pressure for use in the NPT-ensemble.
     *
     * @details This will count the number of particles in the system and
     * calculate the energy contribution to the pressure at a given
     * volume and external pressure.
     *
     * The system energy is
     * \f$\beta u = \beta pV - \ln V - N\ln V\f$
     * where the first two terms are returned by `external()` while the last
     * term is obtained by summing `g_external()` over molecular
     * and or atomic groups.
     * If applied on an atomic group, `N` will be set to the number of
     * atoms in the group, while for a molecular group `N=1`.
     *
     * @date Lund, 2011
     */
    template<class Tspace>
      class ExternalPressure : public Energybase<Tspace> {
        private:
          typedef typename Energybase<Tspace>::Tpvec Tpvec;
          double P; //!< Pressure, p/kT
          string _info() {
            std::ostringstream o;
            o << textio::pad(textio::SUB,15,"Pressure")
              << P*1e30/pc::Nav << " mM = "
              << P*pc::kB*pc::T()*1e30 << " Pa = "
              << P*pc::kB*pc::T()*1e30/0.980665e5 << " atm\n";
            return o.str();
          }
        public:
          ExternalPressure(InputMap &in) {
            this->name="External Pressure";
            P=in.get<double>("npt_P",0,
                "NPT external pressure P/kT (mM)")/1e30*pc::Nav;
            assert(P>=0 && "Specify non-negative pressure");
          }

          /** @brief External energy working on system. pV/kT-lnV */
          double external(const Tpvec&p) FOVERRIDE {
            double V=this->getSpace().geo.getVolume();
            assert(V>1e-9 && "Volume must be non-zero!");
            return P*V - log(V); 
          }

          double g_external(const Tpvec &p, Group &g) FOVERRIDE {
            int N=g.numMolecules();
            double V=this->getSpace().geo.getVolume();
            return -N*log(V);
          }
      };

    /**
     * @brief Sum additive external potential on particles.
     * @tparam Texpot External potential typically derived from
     *         `Potential::ExternalPotentialBase`
     */
    template<typename Tspace, typename Texpot>
      class ExternalPotential : public Energybase<Tspace> {
        private:
          typedef Energybase<Tspace> base;
          string _info() { return expot.info(); }
        public:
          Texpot expot;
          ExternalPotential(InputMap &in) : expot(in) {
            base::name="External Potential ("+expot.name+")";
          }
          double p_external(const typename base::Tparticle &p) FOVERRIDE {
            return expot(p); 
          }
          double i_external(const typename base::Tpvec &p, int i) FOVERRIDE {
            return p_external( p[i] );
          }
          double g_external(const typename base::Tpvec &p, Group &g) FOVERRIDE {
            double u=0;
            for (auto i : g)
              u+=p_external(p[i]);
            return u;
          }
      };

    /**
     * @brief Constrain two group mass centra within a certain distance interval [mindist:maxdist]
     * @author Mikael Lund
     * @date Lund, 2012
     * @todo Prettify output
     *
     * This energy class will constrain the mass center separation between selected groups to a certain
     * interval. This can be useful to sample rare events and the constraint is implemented as an external
     * group energy that return infinity if the mass center separation are outside the defined range.
     * An arbitrary number of group pairs can be added with the addPair() command, although one would
     * rarely want to have more than one.
     * In the following example,
     * the distance between `mygroup1` and `mygroup2` are constrained to the range `[10:50]` angstrom:
     * @code
     * Energy::Hamiltonian pot;
     * auto nonbonded = pot.create( Energy::Nonbonded<Tpairpot,Tgeometry>(mcp) );
     * auto constrain = pot.create( Energy::MassCenterConstrain(pot.getGeometry()) );
     * constrain->addPair( mygroup1, mygroup2, 10, 50); 
     * @endcode
     */
    template<typename Tspace>
      class MassCenterConstrain : public Energybase<Tspace> {
        private:
          string _info();
          struct data {
            double mindist, maxdist;
          };
          std::map<opair<Faunus::Group*>, data> gmap;
        public:
          MassCenterConstrain(Geometry::Geometrybase&);      //!< Constructor
          void addPair(Group&, Group&, double, double);      //!< Add constraint between two groups
          double g_external(const p_vec&, Group&) FOVERRIDE; //!< Constrain treated as external potential
      };

    /**
     * @brief Calculates the total system energy
     *
     * For a given particle vector, space, and energy class we try to
     * calculate the total energy taking into account inter- and
     * intra-molecular interactions as well as external potentials.
     * While this may not work for all systems it may be a useful
     * first guess.
     * This is the default energy routine for `Move::ParallelTempering`
     * and may also be used for checking energy drifts.
     */
    template<class Tspace, class Tenergy, class Tpvec>
      double systemEnergy(Tspace &spc, Tenergy &pot, const Tpvec &p) {
        pot.setSpace(spc); // ensure pot geometry is in sync with spc
        double u = pot.external(p);
        for (auto g : spc.groupList())
          u += pot.g_external(p, *g) + pot.g_internal(p, *g);
        for (size_t i=0; i<spc.groupList().size()-1; i++)
          for (size_t j=i+1; j<spc.groupList().size(); j++)
            u += pot.g2g(p, *spc.groupList()[i], *spc.groupList()[j]);
        return u;
      }

    /* typedefs */

    // template aliasing requires gcc 4.7+
    // template<class Tpairpot, class Tgeometry>
    //   using PairWise = Nonbonded<Tpairpot, Tgeometry>;

  }//Energy namespace
}//Faunus namespace
#endif
