#include<faunus/faunus.h>

namespace Faunus {
  namespace Potential {

    // Base class for manybody potentials
    class ManybodyBase {
      protected:
        typedef vector<int> Tivec;
        Tivec v; // particle index
        string name;
      public:
        virtual string brief() const {
          std::ostringstream o;
          o << name+". Index: ";
          for (auto i : v)
            o << i << " ";
          return o.str();
        };
        virtual void setIndex(const Tivec &index) {
          v = index;
        }
    };

    class Angular : public ManybodyBase {
      private:
        double c1,c2; // pot parameters
      public:
        Angular(InputMap &in) { name = "angle";}

        /** @brief Energy function (kT) */
        template<class Tgeo, class Tpvec>
          double operator()(Tgeo &g, const Tpvec &p) {
            assert(v.size()==3);
            assert((int)p.size()>*std::max_element(v.begin(), v.end()));
            auto ab = g.vdist(p[v[0]], p[v[1]]);
            auto ac = g.vdist(p[v[0]], p[v[2]]);
            auto bc = g.vdist(p[v[1]], p[v[2]]);
            return c1 * ab.dot(bc) / ac.norm(); // for example...
          }
    };

  }//namespace

  namespace Energy {

    /**
     * @brief Energy class for manybody interactions such as dihedrals and angular potentials
     */
    template<class Tspace>
      class Manybody : public Energybase<Tspace> {
        private:
          string _infosum;
        protected:
          string _info() {
            return _infosum;
          }
          typedef Energybase<Tspace> Tbase;
          typedef typename Tbase::Tpvec Tpvec;

          typedef std::function<double(typename Tbase::Tgeometry&,const Tpvec&)> EnergyFunct;
          vector< vector<EnergyFunct> > list; 

        public:
          Manybody(Tspace &spc) {
            Tbase::setSpace(spc);
            list.resize(spc.p.size());
          }

          /**
           * @brief Associate a manybody potential with particle i
           */
          template<class Tmanybodypot>
          void add(size_t i, const Tmanybodypot &f) {
            assert(i>=0 && i<list.size());
            list[i].push_back(f);
            _infosum+=f.brief()+"\n";
          }

          double i2all(Tpvec &p, int i) FOVERRIDE {
            assert(p.size()==list.size());
            double sum=0;
            for (auto &u : list.at(i))
              sum += u(Tbase::spc->geo,p);
            return sum;
          }
      };

  }//namespace

}//namespace

using namespace Faunus;                   // use Faunus namespace
typedef Space<Geometry::Cuboid> Tspace;   // Type of simulation space
typedef Potential::CoulombLJ Tpair;       // and pair potential

int main() {
  ::atom.includefile("../../examples/minimal.json");     // load atom properties
  InputMap in("../../examples/minimal.input");           // open parameter file for user input
  Energy::Nonbonded<Tspace,Tpair> pot(in);// Hamiltonian, non-bonded only
  Tspace spc(in);                         // Simulation space, particles etc.
  Group salt;                             // Group for salt particles
  salt.addParticles(spc,in);              // Add according to user input
  Move::AtomicTranslation<Tspace> mv(in,pot,spc);// particle move class
  mv.setGroup(salt);                      // move class acts on salt group
  //mv.move(1e5);                           // move salt randomly 100000 times
  cout << spc.info() + pot.info() + mv.info(); // final information

  Energy::Manybody<Tspace> three(spc);
  Potential::Angular ang(in);
  ang.setIndex({0,1,2});
  three.add(0, ang);
  ang.setIndex({0,3,4});
  three.add(0, ang);
  cout << three.i2all(spc.p, 0) << endl;
}
