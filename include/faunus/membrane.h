#ifndef FAU_MEMBRANE
#define FAU_MEMBRANE

#include <faunus/faunus.h>

namespace Faunus {
  template<class Tgeometry=Geometry::Cuboid>
    class DesernoMembrane {
      private:
        class _energy : public Energy::Energybase {
          private:
            string _info() { return fene.brief()+harm.brief(); }
            Potential::FENE fene;
            Potential::Harmonic harm;
          public:
            GroupArray lipids;
            Tgeometry geometry;
            _energy(InputMap &in) : fene(in), harm(in), lipids(3), geometry(in) {
              double sigma=in.get<double>("lipid_sigma", 10);
              double epsilon = in.get<double>("lipid_epsilon", 1);
              double headtail_k=0.5*10*epsilon/(sigma*sigma);
              double headtail_req=4*sigma;
              double fene_k=30*epsilon/(sigma*sigma);
              double fene_rmax=1.5*sigma;
              harm = Potential::Harmonic(headtail_k,headtail_req);
              fene = Potential::FENE(fene_k,fene_rmax);
            }
            double g_internal(const p_vec &p, Group &g) { return i_internal(p, g.front()); }
            double i_internal(const p_vec &p, int i) {
              double u=0;
              if (lipids.find(i)) {
                int k=i-1;              // i=middle?
                if ( (i+1)%3==0) k=i-2; // i=last?
                else if (i%3==0) k=i;   // i=first
                u+=fene(p[k], p[k+1], geometry.sqdist(p[k],p[k+1]));
                u+=fene(p[k+1], p[k+2], geometry.sqdist(p[k+1],p[k+2]));
                u+=harm(p[k], p[k+2], geometry.sqdist(p[k],p[k+2]));
              }
              return u;
            }
        };
      public:
        vector<GroupMolecular> pol;
        DesernoMembrane(InputMap &in, Energy::Hamiltonian &pot, Space &spc) {
          auto ePtr = pot.create( _energy(in) );
          // Add lipids
          int N=in.get("lipid_N",0);
          pol.resize(N);
          string polyfile = in.get<string>("lipid_file", "");
          FormatAAM aam;                                         // AAM structure file I/O
          for (auto &g : pol) {                                  // load polymers
            aam.load(polyfile);
            Geometry::FindSpace f;
            f.find(*spc.geo, spc.p, aam.particles() );           // find empty spot in particle vector
            g = spc.insert( aam.particles() );                   // insert into space
            g.name="lipid";
            spc.enroll(g);
            ePtr->lipids.add(g);

            // dope head groups w. charge
            double rho = 1 / in.get<double>("lipid_areapercharge",0);
            cout << "# Surface charge density: " << 1/rho << " A^2/charge\n";
            double area=2 * ePtr->geometry.len.x * ePtr->geometry.len.y;
            double p = rho*area/N;
            if (slp_global.randOne()<p) {
              spc.p[ g.front() ].charge=-1;
              spc.trial[ g.front() ].charge=-1;
            }
          }
        }
    };
}//namespace
#endif
