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
#endif
#ifdef ENABLE_MPI
#include <faunus/mpi.h>
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

          /** @brief Determines if given particle vector is the trial vector */
          template<class Tpvec>
            bool isTrial(const Tpvec &p) const { return (&p==&spc->trial); }

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

          virtual double p2p(const Tparticle&, const Tparticle&) // Particle-particle energy
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

          /* @brief Total energy of i'th = i2all+i_external+i_internal */
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
            pairpot.setSpace(s);
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
        private:
          string _info() { return pairpot.info(25); }
          typedef Energybase<Tspace> Tbase;
          typedef typename Tbase::Tparticle Tparticle;
          typedef typename Tbase::Tpvec Tpvec;
          bool groupBasedField;

        public:
          typename Tspace::GeometryType geo;
          Tpairpot pairpot;

          NonbondedVector(InputMap &in) : geo(in), pairpot(in) {
            static_assert(
                std::is_base_of<Potential::PairPotentialBase,Tpairpot>::value,
                "Tpairpot must be a pair potential" );
            Tbase::name="Nonbonded N" + textio::squared + " - " + pairpot.name;
            groupBasedField = in.get<bool>("pol_g2g",false,"Field will exclude own group");
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
            std::swap(p[0],p[i]);
            double u=0;
#pragma omp parallel for reduction (+:u)
            for (int j=1; j<(int)p.size(); ++j)
              u+=pairpot( p[0], p[j], geo.vdist(p[0],p[j]) );
            std::swap(p[0],p[i]);
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

            // Include field only from external molecules
            if (groupBasedField) {
              for (auto gi : Tbase::spc->groupList()) // loop over all groups
                for (auto gj : Tbase::spc->groupList()) // loop over all group
                  if (gi!=gj) // discard identical groups (addresss comparison)
                    for (auto i : *gi)  // loop over particle index in 1st group
                      for (auto j : *gj) // loop ove
                        E.col(i) = E.col(i) + pairpot.field(p[j],geo.vdist(p[i],p[j]));
            } else {
              size_t i=0;
              for (auto &pi : p) {
                for (auto &pj : p)
                  if (&pi!=&pj)
                    E.col(i) = E.col(i) + pairpot.field(pj,geo.vdist(pi,pj));
                i++;
              }
            }
          }
      };

    /**
     * @brief Nonbonded with early rejection for infinite energies
     *
     * Useful for potentials with a hard sphere part
     */
    template<class Tspace, class Tpairpot, class Tnonbonded=Energy::Nonbonded<Tspace,Tpairpot> >
      class NonbondedEarlyReject : public Tnonbonded {
        private:
          typedef Tnonbonded base;
        public:
          NonbondedEarlyReject(InputMap &in) : base(in) {
            base::name+=" (early reject)";
          }

          double g2g(const typename base::Tpvec &p, Group &g1, Group &g2) FOVERRIDE {
            double u=0;
            for (auto i : g1)
              for (auto j : g2) {
                auto _u = base::pairpot(p[i],p[j],base::geo.sqdist(p[i],p[j]));
                if (std::isinf(_u))
                  return INFINITY;
                u+=_u;
              }
            return u;
          }

          double g_internal(const typename base::Tpvec &p, Group &g) FOVERRIDE { 
            double u=0;
            int b=g.back(), f=g.front();
            if (!g.empty())
              for (int i=f; i<b; ++i)
                for (int j=i+1; j<=b; ++j) {
                  auto _u = base::pairpot(p[i],p[j],base::geo.sqdist(p[i],p[j]));
                  if (std::isinf(_u))
                    return INFINITY;
                  u+=_u;
                }
            return u;
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
    template<class Tspace, class Tpairpot, class Tnonbonded=Energy::Nonbonded<Tspace,Tpairpot> >
      class NonbondedCutg2g : public Tnonbonded {
        private:
          typedef Tnonbonded base;
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
            << std::boolalpha << CrossGroupBonds << endl
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
          int i=spc->findIndex(a);
          int j=spc->findIndex(b);
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

        /* Commented to allow for gcc 4.6 compatibility
           auto getBondList() -> decltype(this->list) {
           return Tbase::getBondList();
           }
           */

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

          /** @brief Field on all particles due to external potential */
          void field(const typename base::Tpvec &p, Eigen::MatrixXd &E) FOVERRIDE {
            assert((int)p.size()==E.cols());
            for (size_t i=0; i<p.size(); i++)
              E.col(i) += expot.field(p[i]);
          }
      };

    /**
     * @brief Restrain a group mass center within a certain z-axis interval [bin_min:bin_max]
     * @author Joao Henriques
     * @date Lund, 2013     
     * 
     * This energy class will restrain the mass center of a given group within a certain window/bin along the z-axis.
     * Mainly for use with free energy vs. surface distance simulations, where the the surface is too attractive and prevents
     * correct sampling. Can also be used to restrain the mass center of a group to a subset of the simulation cell volume in 
     * other type of systems.
     *
     * The InputMap parameters are:
     *
     * Key                | Description
     * :----------------- | :---------------------------
     * `bin_min`          | Lower limit (always positive, i.e. from 0 to `cuboid_zlen`), [angstrom]
     * `bin_max`          | Higher limit (from `bin_min` to `cuboid_zlen`), [angstrom]
     *
     */
    template<class Tspace>
      class MassCenterRestrain : public Energybase<Tspace> {
        private:
          typedef Energybase<Tspace> base;
          typedef typename Tspace::p_vec Tpvec;
          string _info() {
            std::ostringstream o;
            o << pad(textio::SUB,25,"Bin limits (z-axis)") << "[" << min << ":" << max << "]" << textio::_angstrom << endl;
            return o.str();
          };
          double min, max;
          Group* gPtr;
        public:
          MassCenterRestrain(InputMap &in) {
            min = in.get<double>("bin_min", 0);
            max = in.get<double>("bin_max", pc::infty);
            base::name = "Mass Center Restrain";
            gPtr = nullptr;
          }
          void add(Group &g) {
            gPtr = &g;
          }
          double g_external(const Tpvec &p, Group &g) {
            if (&g != gPtr)
              return 0;
            double boxlenz = base::spc->geo.len.z();
            double mc = Geometry::massCenter(base::spc->geo, p, g).z() + (boxlenz / 2);
            if (mc >= min && mc <= max)
              return 0;
            else
              return pc::infty;
          }
          double i_external(const Tpvec &p, int i) {
            auto gi = base::spc->findGroup(i);
            return g_external(p, *gi);
          }
      };

    /**
     * @brief Constrain two group mass centre separations to a distance interval [mindist:maxdist]
     *
     * This energy class will constrain the mass center separation between selected groups to a certain
     * interval. This can be useful to sample rare events and the constraint is implemented as an external
     * group energy that return infinity if the mass center separation are outside the defined range.
     * An arbitrary number of group pairs can be added with the addPair() command, although one would
     * rarely want to have more than one.
     * In the following example, the distance between `mygroup1` and `mygroup2`
     * are constrained to the range `[10:50]` angstrom:
     *
     * @code
     * InputMap mcp;
     * Energy::MassCenterConstrain<Tspace> pot(mcp);
     * pot.addPair( mygroup1, mygroup2, 10, 50); 
     * @endcode
     *
     * The `addPair` function can be called without distance interval in
     * which case the default window is used. This is read during from
     * the `InputMap` using the keywords,
     *
     * Keyword           | Description
     * :---------------- | :----------------------------------------
     * `cmconstrain_min` | Minimum mass center separation (angstrom)
     * `cmconstrain_max` | Maximum mass center separation (angstrom)
     *
     * @date Lund, 2012
     * @todo Prettify output
     */
    template<typename Tspace>
      class MassCenterConstrain : public Energybase<Tspace> {
        private:
          string _info() {
            using namespace Faunus::textio;
            std::ostringstream o;
            o << indent(SUB) << "The following groups have mass center constraints:\n";
            for (auto &m : gmap)
              o << indent(SUBSUB) << m.first.first->name << " " << m.first.second->name
                << " " << m.second.mindist << "-" << m.second.maxdist << _angstrom << endl;
            return o.str();
          }

          struct data {
            double mindist, maxdist;
          };

          data defaultWindow;

          std::map<opair<Faunus::Group*>, data> gmap;

        public:

          MassCenterConstrain(InputMap &in) {
            this->name="Group Mass Center Distance Constraints";
            defaultWindow.mindist = in("cmconstrain_min", 0.0);
            defaultWindow.maxdist = in("cmconstrain_max", 1.0e20);

            assert(defaultWindow.mindist<defaultWindow.maxdist);
            assert(defaultWindow.mindist>0);
            assert(defaultWindow.maxdist>0);
          }

          /** @brief Add constraint between two groups */
          void addPair(Group &a, Group &b, double mindist, double maxdist) {
            data d = {mindist, maxdist};
            opair<Group*> p(&a, &b);
            gmap[p] = d;
          }

          /** @brief Add constraint between two groups - distances read from user input */
          void addPair(Group &a, Group &b) {
            addPair(a,b,defaultWindow.mindist,defaultWindow.maxdist);
          }

          /** @brief Constrain treated as external potential */
          double g_external(const p_vec &p, Group &g1) {
            for (auto m : gmap)              // scan through pair map
              if (m.first.find(&g1)) {       // and look for group g1
                Group *g2ptr;                // pointer to constrained partner
                if (&g1 == m.first.first)
                  g2ptr = m.first.second;
                else
                  g2ptr = m.first.first;
                Point cma = Geometry::massCenter(this->getSpace().geo, p, g1);
                Point cmb = Geometry::massCenter(this->getSpace().geo, p, *g2ptr);
                double r2 = this->getSpace().geo.sqdist(cma,cmb);
                double min = m.second.mindist;
                double max = m.second.maxdist;
                if (r2<min*min || r2>max*max) 
                  return pc::infty;
              }
            return 0;
          }
      };

    /**
     * @brief Hydrophobic interactions using SASA
     *
     * This will add a square-well attraction between solvent
     * exposed, hydrophobic sites where the strength is given
     * by a surface tension argument, 
     *
     * @f[ U = \sum_i^{N_{contact}} \gamma a_i @f]
     *
     * where @f$a@f$ is the SASA value for a residue and `i` runs
     * over all particles in contact with at least one other
     * hydrophobic and exposed site(s). Note that the value of the surface
     * tension, @f$\gamma@f$, is not comparable to a true, macroscopic
     * tension.
     *
     * Upon construction the following keywords are read from `InputMap`:
     *
     *  Keyword                | Comment
     *  :--------------------- | :--------------
     *  `sasahydro_sasafile`   | SASA file - one column (angstrom^2)
     *  `sasahydro_duplicate`  | read SASA file n times (default: 1)
     *  `sasahydro_tension`    | surface tension (default: 0 dyne/cm)
     *  `sasahydro_threshold`  | surface distance threshold (default: 3.0 angstrom)
     *  `sasahydro_uofr`       | set to "yes" if U(r) should be sampled (default: "no")
     *  `sasahydro_dr`         | distance resolution of U(r) (default: 0.5 angstrom)
     *
     *  The `sasafile` should be a single column file with SASA values for each particle in
     *  the system (angstrom^2). Alternatively one can load the file several times using
     *  the `duplicate` option. To generate a SASA file of a protein, use the vmd-sasa.tcl
     *  VMD script, found in the `scripts` folder.
     *
     *  If `sasahydro_uofr` is specified, the hydrophobic interaction
     *  energy is averaged as a function of mass center separation between groups and
     *  saved to disk upon calling the destructor.
     *
     * @todo
     * Currently only `g2g` is implemented and it is assumed that a hydrophobic site on
     * one group is in contact with sites on *maximum one* other group. For large macro 
     * molecules this is hardly a problem; if it is, an energy drift should show.
     *
     * @author Kurut / Lund
     * @date Lund, 2014
     * @note Experimental
     */
    template<class Tspace>
      class HydrophobicSASA : public Energybase<Tspace> {
        private:
          typedef Energybase<Tspace> base;
          typedef typename Tspace::p_vec Tpvec;
          typedef std::map<int, Average<double> > Tuofr;

          std::map<string,Tuofr> uofr; // sasa energy vs. group-2-group distance

          vector<bool> v;      // bool for all particles; true=active
          vector<double> sasa; // sasa for all particles
          double threshold;    // surface-surface distance threshold
          double tension, tension_dyne;
          int duplicate;       // how many times to load sasa file
          string file;         // sasa file

          bool sample_uofr;    // set to true if we should sample U_sasa(r)
          double dr;           // U(r) resolution in r

          /** @brief Fraction of hydrophobic area */
          double fracHydrophobic() const {
            double h=0, noh=0;
            if (sasa.size()==base::spc->p.size()) {
              for (size_t i=0; i<sasa.size(); i++)
                if (base::spc->p[i].hydrophobic)
                  h+=sasa[i];
                else
                  noh+=sasa[i];
              return h/(h+noh);
            }
            return 0;
          }

          string _info() {
            char w=25;
            using namespace textio;
            std::ostringstream o;
            o << pad(SUB,w,"SASA file") << file << " (duplicated " << duplicate << " times)\n"
              << pad(SUB,w,"SASA vector size") << sasa.size()
              << " (particle vector = " << base::spc->p.size() << ")\n"
              << pad(SUB,w,"Hydrophobic SASA") << fracHydrophobic()*100 << percent + "\n"
              << pad(SUB,w,"Threshold") << threshold << _angstrom + "\n"
              << pad(SUB,w,"Surface tension") << tension_dyne << " dyne/cm = " << tension
              << kT+"/"+angstrom+squared+"\n"
              << pad(SUB,w+1,"Sample "+beta+"U(r)") << std::boolalpha << sample_uofr << "\n"
              << pad(SUB,w+1,beta+"U(r) resolution") << dr << _angstrom + "\n";
            return o.str();
          }

          void loadSASA(const string &file, unsigned int n) {
            double s;
            std::vector<double> t;
            std::ifstream f(file);
            while (f>>s)
              t.push_back(s);
            sasa.clear();
            while (n-->0)
              sasa.insert(sasa.end(), t.begin(), t.end());
            if (sasa.empty())
              std::cerr << "Warning: No SASA data loaded from "+file+"\n";
          }

          /** @brief Save U(r) to disk if sampled */
          void save(const string &pfx=textio::prefix) {
            for (auto &i : uofr) {
              std::ofstream f(pfx+"Usasa_"+i.first);
              if (f)
                for (auto &j : i.second)
                  f << j.first*dr << " " << j.second << "\n";
            }
          }

        public:
          HydrophobicSASA(InputMap &in) {
            string pfx="sasahydro_";
            base::name = "Hydrophobic SASA";
            threshold = in.get<double>(pfx+"threshold", 3.0);
            tension_dyne = in.get<double>(pfx+"tension", 0);
            tension = tension_dyne * 1e-23 / (pc::kB * pc::T());  // dyne/cm converted to kT/A^2, 1dyne/cm = 0.001 J/m^2
            file = in.get(pfx+"sasafile", string());
            duplicate = in.get(pfx+"duplicate", 0);
            sample_uofr = in.get<bool>(pfx+"uofr", false);
            dr = in.get<double>(pfx+"dr", 0.5);
            loadSASA(file, duplicate);
          }

          ~HydrophobicSASA() { save(); }

          /** @brief Group-to-group energy */
          double g2g(const Tpvec &p, Group &g1, Group &g2) FOVERRIDE {
            double dsasa=0;
            if (sasa.size()==p.size())
              if (g1.isMolecular())
                if (g2.isMolecular()) {
                  v.resize(p.size());
                  std::fill(v.begin(), v.end(), true);
                  for (auto i : g1)
                    for (auto j : g2)
                      if (p[i].hydrophobic)
                        if (p[j].hydrophobic)
                          if (sasa[i]>1e-3)
                            if (sasa[j]>1e-3)
                              if (v[i] || v[j]) {
                                double r2=base::spc->geo.sqdist(p[i],p[j]);
                                if (r2<pow(threshold+p[i].radius+p[j].radius,2)) {
                                  if (v[i])
                                    dsasa += sasa[i];
                                  if (v[j])
                                    dsasa += sasa[j];
                                  v[i]=v[j]=false;
                                }
                              }
                  // analyze
                  if (sample_uofr && !base::isTrial(p)) {
                    double r = base::spc->geo.dist(g1.cm,g2.cm);
                    auto k = std::minmax(g1.name,g2.name);
                    uofr[k.first+"-"+k.second][to_bin(r,dr)] += -tension*dsasa;
                  }
                }

            return -tension*dsasa; // in kT
          }
      };

#ifdef ENABLE_MPI
    template<typename Tcoord=double>
      class PenaltyFunction1D : public Table2D<Tcoord,double> {
        private:
          int _cnt;
          int _Nupdate;
          double _du;
          double _du_sum;
          int _size; // maximum number of keys in the map
          float _lo1, _hi1;
          typedef Faunus::MPI::FloatTransmitter::floatp floatp;
          typedef Table2D<Tcoord,double> Tbase;
          typedef Table2D<Tcoord,double> Thist;
          typedef std::pair<Tcoord,Tcoord> Tpair;
          Thist hist;
          Faunus::MPI::MPIController *mpiPtr; 
          Faunus::MPI::FloatTransmitter ft;
        public:
          /**
           * @brief Constructor
           *
           * @param mpi MPI controller
           * @param Nupdate Number of moves between updates of the penalty function
           * @param size Total number of points in the penalty function (default 2000)
           * @param bw1 Bin width of 1st reaction coordinate (default 1)
           * @param bw1 Bin width of 2nd reaction coordinate (default 1)
           * @param lo1 Lower limit of 1st reaction coordinate (default -20)
           * @param hi1 Lower limit of 1st reaction coordinate (default 20)
           * @param lo2 Lower limit of 2nd reaction coordinate (default -20)
           * @param hi2 Lower limit of 2nd reaction coordinate (default 20)
           */
          PenaltyFunction1D(Faunus::MPI::MPIController &mpi, int Nupdate, size_t size, Tcoord res1, Tcoord res2,
              float lo1, float hi1, float lo2, float hi2)
            : Tbase(res1, Tbase::XYDATA), hist(res1, Thist::XYDATA), mpiPtr(&mpi) {
              Tbase::name = "penalty1D";
              _size = size*2;
              _cnt = 0;
              _du = 0;
              _Nupdate = Nupdate;
              _du_sum = 0;
              _lo1 = lo1;
              _hi1 = hi1;
            }
          /** @brief Check if coordinate is within user-defined range */
          bool isInrange(Tcoord &coord) {
            return (coord>=_lo1 && coord<=_hi1);
          }
          /**
           * @brief Merge histograms obtained from parallel processes with different seeds
           *
           * @details Slave processes send histograms to the master. The master computes the 
           * sum over all histograms and sends it back to the slaves.
           */
          void exchange() {
            if (!mpiPtr->isMaster()) {
              std::vector<floatp> sendBuf = hist.hist2buf(_size);
              std::vector<floatp> recvBuf = ft.swapf(*mpiPtr, sendBuf, mpiPtr->rankMaster());
              hist.buf2hist(recvBuf);
            }
            if (mpiPtr->isMaster()) {
              std::vector<floatp> sendBuf = hist.hist2buf(_size);
              std::vector<floatp> recvBuf(_size);
              for (int i=0; i<mpiPtr->nproc(); ++i) {
                if (i==mpiPtr->rankMaster()) continue;
                ft.recvf(*mpiPtr, i, recvBuf);
                ft.waitrecv();
                sendBuf.insert(sendBuf.end(), recvBuf.begin(), recvBuf.end());
              }
              hist.buf2hist(sendBuf);
              sendBuf = hist.hist2buf(_size);
              for (int i=1; i<mpiPtr->nproc(); ++i) {
                ft.sendf(*mpiPtr, sendBuf, i);
                ft.waitsend();
              }
              hist.buf2hist(sendBuf);
            }
          }
          /** @brief Update histogram of single processes and penalty function 
           * 
           *  @details Penalty function is updated using histogram obtained
           *  by merging contributions from all processes.
           */
          double update(Tcoord &coord) {
            _cnt++;
            hist(coord)++; // increment internal histogram
            _du = 0;
            if ((_cnt%_Nupdate)==0) { // if Nupdate'th time
              exchange();
              for (auto &m : hist.getMap()) { // update penalty function
                Tbase::operator()(m.first) += std::log(m.second);
                double log2 = std::log(2); // half bin width compensation
                if (m.first==_lo1 || m.first==_hi1) 
                  Tbase::operator()(m.first) += log2;
              }
              double min = Tbase::min()->second;
              for (auto &m : Tbase::getMap()) {
                Tbase::operator()(m.first) -= min;
              }
              assert(hist(coord)!=0 && "hist_size (>= max # of points in histogram) is set too small.");
              _du = std::log(hist(coord)) - min; // prevents energy drift
              _du_sum += _du;
              hist.clear();
            }
            return _du;
          }
          /**
           * @brief Update histogram of single processes using waste-recycling method.
           *
           * @details If energy change > 0, histogram receives contributions from both
           * accepted and rejected configurations.
           * Penalty function is updated using histogram obtained
           * by merging contributions from all processes.
           *
           * [More info](http://dx.doi.org/10.1007/3-540-35273-2_4)
           */
          double update(Tpair &coord, double &weight, bool &rejection) {
            _cnt++;
            if (weight < 1) {
              hist(coord.first) += (1-weight);
              hist(coord.second) += weight;
            }
            else hist(coord.second)++;
            _du = 0;
            if ((_cnt%_Nupdate)==0) { // if Nupdate'th time
              exchange();
              for (auto &m : hist.getMap()) { // update penalty function
                Tbase::operator()(m.first) += std::log(m.second);
                double log2 = std::log(2); // half bin width compensation
                if (m.first==_lo1 || m.first==_hi1) 
                  Tbase::operator()(m.first) += log2;
              }
              double min = Tbase::min()->second;
              for (auto &m : Tbase::getMap()) {
                Tbase::operator()(m.first) -= min;
              }
              if (!rejection) _du = std::log(hist.find(coord.second)) - min;
              else if (1-weight>1e-20) _du = std::log(hist.find(coord.first)) - min;
              assert(hist(coord.second)!=0 && "hist_size (>= max # of points in histogram) is set too small.");
              _du_sum += _du;
              hist.clear();
            }
            return _du;
          }
          /** @brief Save table to disk */
          void save(const string &filename) {
            Tbase::save(filename+"penalty");
            hist.save(filename+"histo");
          }
          /** 
           * @brief Translate penalty by a reference value and save it to disk. 
           *
           * @details The reference value is obtained by averaging over a region of the
           * reaction coordinate space.
           */
          void save_final(const string &filename, 
              Tcoord a, Tcoord b, Tcoord c, Tcoord d) {
            double ref_value = - Tbase::ave(a,b);
            Tbase::save(filename+"penalty_fin",1.,ref_value);
          }

          /** @brief Load table to disk */
          void load(const string &filename) {
            Tbase::load(filename+"penalty");
            hist.load(filename+"histo");
          }
          void test(UnitTest &t) {
            if (_cnt/_Nupdate>=1) {
              auto it_max = Tbase::max();
              auto it_min = Tbase::min();
              t("penalty2D_range",it_max->second-it_min->second,0.1);
            }
          }
          string info() {
            using namespace Faunus::textio;
            std::ostringstream o;
            o << header("1D penalty function");
            if (_cnt/_Nupdate>=1) {
              char w=25;
              auto it_max = Tbase::max();
              auto it_min = Tbase::min();
              o << textio::pad(SUB,w, "Number of updates") << _cnt/_Nupdate << endl
                << textio::pad(SUB,w, "Point of lowest penalty") << it_min->first 
                << " " << angstrom << endl
                << textio::pad(SUB,w, "Point of highest penalty") << it_max->first
                << " " << angstrom << endl
                << textio::pad(SUB,w, "Penalty range") << it_max->second-it_min->second << kT << endl
                << textio::pad(SUB,w, "Drift compensation") << _du_sum << kT << endl;
            }
            return o.str();
          }
      };

    template<typename Tcoord=double>
      class PenaltyFunction2D : public Table3D<Tcoord,double> {
        private:
          int _cnt;
          int _Nupdate;
          double _du;
          double _du_sum;
          int _size; // maximum number of keys in the map
          float _lo1, _hi1, _lo2, _hi2;
          typedef Faunus::MPI::FloatTransmitter::floatp floatp;
          typedef Table3D<Tcoord,double> Tbase;
          typedef Table3D<Tcoord,double> Thist;
          typedef std::pair<Tcoord,Tcoord> Tpair;
          Thist hist;
          Faunus::MPI::MPIController *mpiPtr; 
          Faunus::MPI::FloatTransmitter ft;
        public:
          /**
           * @brief Constructor
           *
           * @param mpi MPI controller
           * @param Nupdate Number of moves between updates of the penalty function
           * @param size Total number of points in the penalty function (default 2000)
           * @param bw1 Bin width of 1st reaction coordinate (default 1)
           * @param bw1 Bin width of 2nd reaction coordinate (default 1)
           * @param lo1 Lower limit of 1st reaction coordinate (default -20)
           * @param hi1 Lower limit of 1st reaction coordinate (default 20)
           * @param lo2 Lower limit of 2nd reaction coordinate (default -20)
           * @param hi2 Lower limit of 2nd reaction coordinate (default 20)
           */
          PenaltyFunction2D(Faunus::MPI::MPIController &mpi, int Nupdate, size_t size, Tcoord res1, Tcoord res2,
              float lo1, float hi1, float lo2, float hi2)
            : Tbase(res1, res2, Tbase::XYDATA), hist(res1, res2, Thist::XYDATA), mpiPtr(&mpi) {
              Tbase::name = "penalty2D";
              _size = size*3;
              _cnt = 0;
              _du = 0;
              _Nupdate = Nupdate;
              _du_sum = 0;
              _lo1 = lo1;
              _hi1 = hi1;
              _lo2 = lo2;
              _hi2 = hi2;
            }
          /** @brief Check if coordinates are within user-defined ranges */
          bool isInrange(Tpair &coord) {
            return (coord.first>=_lo1 && coord.first<=_hi1 && coord.second>=_lo2 && coord.second<=_hi2);
          }
          /**
           * @brief Merge histograms obtained from parallel processes with different seeds
           *
           * @details Slave processes send histograms to the master. The master computes the 
           * sum over all histograms and sends it back to the slaves.
           */
          void exchange() {
            if (!mpiPtr->isMaster()) {
              std::vector<floatp> sendBuf = hist.hist2buf(_size);
              std::vector<floatp> recvBuf = ft.swapf(*mpiPtr, sendBuf, mpiPtr->rankMaster());
              hist.buf2hist(recvBuf);
            }
            if (mpiPtr->isMaster()) {
              std::vector<floatp> sendBuf = hist.hist2buf(_size);
              std::vector<floatp> recvBuf(_size);
              for (int i=0; i<mpiPtr->nproc(); ++i) {
                if (i==mpiPtr->rankMaster()) continue;
                ft.recvf(*mpiPtr, i, recvBuf);
                ft.waitrecv();
                sendBuf.insert(sendBuf.end(), recvBuf.begin(), recvBuf.end());
              }
              hist.buf2hist(sendBuf);
              sendBuf = hist.hist2buf(_size);
              for (int i=1; i<mpiPtr->nproc(); ++i) {
                ft.sendf(*mpiPtr, sendBuf, i);
                ft.waitsend();
              }
              hist.buf2hist(sendBuf);
            }
          }
          /** @brief Update histogram of single processes and penalty function
           *                               
           *  @details Penalty function is updated using histogram obtained 
           *  by merging contributions from all processes.                               
           */ 
          double update(Tpair &coord) {
            _cnt++;
            hist(coord.first, coord.second)++; // increment internal histogram
            _du = 0;
            if ((_cnt%_Nupdate)==0) { // if Nupdate'th time
              exchange();
              for (auto &m : hist.getMap()) { // update penalty function
                Tbase::operator()(m.first.first, m.first.second) += std::log(m.second);
                double log2 = std::log(2); // half bin width compensation
                if (m.first.first==_lo1 || m.first.first==_hi1) {
                  Tbase::operator()(m.first.first, m.first.second) += log2;
                  if (m.first.second==_lo2 || m.first.second==_hi2)
                    Tbase::operator()(m.first.first, m.first.second) += log2;
                }
                else if (m.first.second==_lo2 || m.first.second==_hi2)
                  Tbase::operator()(m.first.first, m.first.second) += log2;
              }
              double min = Tbase::min()->second;
              for (auto &m : Tbase::getMap()) {
                Tbase::operator()(m.first.first, m.first.second) -= min;
              }
              assert(hist(coord.first, coord.second)!=0 && "hist_size (>= max # of points in histogram) is set too small.");
              if (hist(coord.first, coord.second)==0) cout << mpiPtr->rank() << " hist(coord1, coord2)=0\n";
              _du = std::log(hist(coord.first, coord.second)) - min; // prevents energy drift
              _du_sum += _du;
              hist.clear();
            }
            return _du;
          }
          /**
           * @brief Update histogram of single processes using waste-recycling method.
           *
           * @details If energy change > 0, histogram receives contributions from both
           * accepted and rejected configurations.
           * Penalty function is updated using histogram obtained
           * by merging contributions from all processes.
           *
           * [More info](http://dx.doi.org/10.1007/3-540-35273-2_4)
           */
          double update(std::pair<Tpair,Tpair> &coord, double &weight, bool &rejection) {
            _cnt++;
            if (weight < 1) {
              hist(coord.first.first, coord.first.second) += (1-weight);
              hist(coord.second.first, coord.second.second) += weight;
            }
            else hist(coord.second.first, coord.second.second)++;
            _du = 0;
            if ((_cnt%_Nupdate)==0) { // if Nupdate'th time
              exchange();
              for (auto &m : hist.getMap()) { // update penalty function
                Tbase::operator()(m.first.first, m.first.second) += std::log(m.second);
                double log2 = std::log(2); // half bin width compensation
                if (m.first.first==_lo1 || m.first.first==_hi1) {
                  Tbase::operator()(m.first.first, m.first.second) += log2;
                  if (m.first.second==_lo2 || m.first.second==_hi2)
                    Tbase::operator()(m.first.first, m.first.second) += log2;
                }
                else if (m.first.second==_lo2 || m.first.second==_hi2)
                  Tbase::operator()(m.first.first, m.first.second) += log2;
              }
              double min = Tbase::min()->second;
              for (auto &m : Tbase::getMap()) {
                Tbase::operator()(m.first.first, m.first.second) -= min;
              }
              if (!rejection) _du = std::log(hist.find(coord.second)) - min;
              else if (1-weight>1e-20) _du = std::log(hist.find(coord.first)) - min;
              assert(hist(coord.second.first, coord.second.second)!=0 
                  && "hist_size (>= max # of points in histogram) is set too small.");
              _du_sum += _du;
              hist.clear();
            }
            return _du;
          }

          /** @brief Save table to disk */
          void save(const string &filename) {
            Tbase::save(filename+"penalty");
            hist.save(filename+"histo");
          }
          /** 
           * @brief Translate penalty by a reference value and save it to disk. 
           *
           * @details The reference value is obtained by averaging over a region of the
           * reaction coordinate space.
           */
          void save_final(const string &filename,
              Tcoord a, Tcoord b, Tcoord c, Tcoord d) {
            double ref_value = - Tbase::ave(a,b,c,d);
            Tbase::save(filename+"penalty_fin",1.,ref_value);
          }
          /** @brief Load table to disk */
          void load(const string &filename) {
            Tbase::load(filename+"penalty");
            hist.load(filename+"histo");
          }
          void test(UnitTest &t) {
            if (_cnt/_Nupdate>=1) {
              auto it_max = Tbase::max();
              auto it_min = Tbase::min();
              t("penalty2D_range",it_max->second-it_min->second,0.1);
            }
          }
          string info() {
            using namespace Faunus::textio;
            std::ostringstream o;
            o << header("2D penalty function");
            if (_cnt/_Nupdate>=1) {
              char w=25;
              auto it_max = Tbase::max();
              auto it_min = Tbase::min();
              o << textio::pad(SUB,w, "Number of updates") << _cnt/_Nupdate << endl
                << textio::pad(SUB,w, "Point of lowest penalty") << "(" << it_min->first.first 
                << ", " << it_min->first.second << ") " << angstrom << endl
                << textio::pad(SUB,w, "Point of highest penalty") << "(" << it_max->first.first 
                << ", " << it_max->first.second << ") " << angstrom << endl
                << textio::pad(SUB,w, "Penalty range") << it_max->second-it_min->second << kT << endl
                << textio::pad(SUB,w, "Drift compensation") << _du_sum << kT << endl;
            }
            return o.str();
          }
      };

    /**
      @brief General class for penalty functions along a coordinate
      @date Malmo, 2011

      This class stores a penalty function, f(x), along a given coordinate, x,
      of type `Tcoordinate` which could be a distance, angle, volume etc.
      Initially f(x) is zero for all x.
      Each time the system visits x the update(x) function should be called
      so as to add the penalty energy, du. In the energy evaluation, the
      coordinate x should be associated with the extra energy f(x).

      This will eventually ensure uniform sampling. Example:

      ~~~
      PenaltyFunction<double> f(0.1,1000,6.0); // 0.1 kT penalty
      Point masscenter;           // some 3D coordinate...
      ...
      f.update(masscenter.z);     // update penalty energy for z component
      double u = f(masscenter.z); // get accumulated penalty at coordinate (kT)
      f.save("penalty.dat");      // save to disk
      ~~~

      In the above example, the penalty energy will be scaled by 0.5 if the
      sampling along the coordinate is less than 6 kT between the least and
      most likely position.
      This threshold check is carried out every 1000th call to `update()`.
      Note also that when the penalty energy is scaled, so is the threshold
      (also by a factor of 0.5).
      */
    template<typename Tcoord=float>
      class PenaltyFunction : public Table2D<Tcoord,double> {
        private:
          unsigned long long _cnt;
          int _Ncheck;
          double _kTthreshold;
          typedef Table2D<Tcoord,double> Tbase;
          typedef Table2D<Tcoord,unsigned long long int> Thist;
          Thist hist;
          Tcoord _du; //!< penalty energy
          std::string _log;
        public:
          /**
           * @brief Constructor
           * @param penalty Penalty energy for each update (kT)
           * @param Ncheck Check histogram every Nscale'th step
           *        (put large number for no scaling, default)
           * @param kTthreshold Half penalty energy once this
           *        threshold in distribution has been reached
           * @param res Resolution of the penalty function (default 0.1)
           */
          PenaltyFunction(double penalty, int Ncheck=1e9, double kTthreshold=5, Tcoord res=0.1)
            : Tbase(res, Tbase::XYDATA), hist(res, Thist::HISTOGRAM) {
              Tbase::name="penalty";
              _cnt=0;
              _Ncheck=Ncheck;
              _kTthreshold=kTthreshold;
              _du=penalty;
              assert(Ncheck>0);
              _log="#   initial penalty energy = "+std::to_string(_du)+"\n";
            }
          /** @brief Update penalty for coordinate */
          double update(Tcoord coordinate) {
            _cnt++;
            Tbase::operator()(coordinate)+=_du;  // penalize coordinate
            hist(coordinate)++;                  // increment internal histogram
            if ((_cnt%_Ncheck)==0) {             // if Ncheck'th time
              double deltakT=log( hist(hist.maxy()) / double(hist(hist.miny())) );
              assert(deltakT>0);
              std::ostringstream o;
              o << "#   n=" << _cnt << " dkT=" << deltakT;
              if (deltakT<_kTthreshold) {   // if histogram diff. is smaller than threshold
                _kTthreshold*=0.5;          // ...downscale threshold
                scale(0.5);                 // ...and penalty energy
                o << " update: du=" << _du << " threshold=" << _kTthreshold;
              }
              _log += o.str() + "\n";       // save info to log
            }
            return _du;
          }
          /*! \brief Manually scale penalty energy */
          void scale(double s) { _du*=s; }

          /*! \brief Save table to disk */
          void save(const string &filename) {
            Tbase::save(filename);
            hist.save(filename+".dist");
          }

          string info() {
            return "# Penalty function log:\n" + _log;
          }
      };

    /**
     * @brief General energy class for handling penalty function (PF) methods.
     *
     * @details User-defined reaction coordinate(s) are provided by an external function, f.
     * f is a member of a class given as a template argument.
     * The return type of f can be either a double or a pair of doubles.
     * The dimensionality of the penalty function is inferred from the return value of f.
     * The InputMap class is scanned for the following keys:
     *
     * Key              | Description
     * :--------------- | :-----------------------------
     * `penalty_update` | number of moves between updates of PF
     * `penalty_size`   | total number of points in the PF
     * `penalty_bw1`    | bin width of 1st coordinate
     * `penalty_bw2`    | bin width of 2nd coordinate
     * `penalty_lo1`    | lower limit of 1st coordinate
     * `penalty_hi1`    | upper limit of 1st coordinate
     * `penalty_lo2`    | lower limit of 2nd coordinate
     * `penalty_hi2`    | upper limit of 2nd coordinate
     */
    template<class Tspace, class Tfunction>
       class PenaltyEnergy : public Energybase<Tspace> {
         private:
           string _info() { return "Energy from Penalty Function\n"; }
           Tspace* spcPtr;
           string pfx = "penalty";
           Tfunction f;
           typedef decltype(f(spcPtr->p)) Treturn; // type of coordinate
           typedef Energybase<Tspace> Tbase;
           typedef typename Energy::PenaltyFunction1D<double> Tone;
           typedef typename Energy::PenaltyFunction2D<double> Ttwo;
           typedef typename std::conditional<std::is_same<double,Treturn>::value,Tone,Ttwo>::type Tpf;
           Tpf pf;
         public:
           std::pair<Treturn,Treturn> coordpair; // current and trial coordinates
           PenaltyEnergy(Faunus::MPI::MPIController &mpi, InputMap &in) :
             pf(mpi,in.get<int>(pfx+"_update",1e5),
                 in.get<size_t>(pfx+"_size",2000),
                 in.get<float>(pfx+"_bw1",0.1),
                 in.get<float>(pfx+"_bw2",0.1),
                 in.get<float>(pfx+"_lo1",-20),
                 in.get<float>(pfx+"_hi1",20),
                 in.get<float>(pfx+"_lo2",-20),
                 in.get<float>(pfx+"_hi2",20)
               ) {}
           string info() { return pf.info(); }
           void test(UnitTest &t) { pf.test(t); }
           void load(const string &filename) { pf.load(filename); }
           void save(const string &filename) { pf.save(filename); }
           void save_final(const string &filename, double a, double b, double c=0, double d=0) { 
             pf.save_final(filename, a, b, c, d); 
           }
           double update(Treturn &c) {
             return pf.update(c);
           }
           double update(std::pair<Treturn,Treturn> &c, double &w, bool &r) {
             return pf.update(c, w, r);
           }
           std::map<Treturn,double> getMap() { return pf.getMap(); }
           double find(Treturn c) { return pf.find(c); }
           double external(const p_vec &p) {
             Treturn coor = f(p);
             double du;
             if (!Tbase::isTrial(p)) coordpair.first=coor; // trial coordinate is stored
             else coordpair.second=coor; // current coordinate is stored
             if (!pf.isInrange(coor)) du = std::numeric_limits<double>::max();
             else du = pf.find(coor);
             return du;
           }
       };
#endif

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
