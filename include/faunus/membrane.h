#ifndef FAU_MEMBRANE
#define FAU_MEMBRANE

#include <faunus/faunus.h>

namespace Faunus {

  namespace Potential {
    // Mixes any pair potential with CosAttract, only between "TL" type particles.
    template<class Tpairpot>
      class CosAttractCombi : public CombinedPairPotential<Tpairpot,CosAttract> {
        private:
          typedef CombinedPairPotential<Tpairpot,CosAttract> Tbase;
        public:
          particle::Tid id; // T2 potential added to particle pairs of this type only!
          CosAttractCombi(InputMap &in) : Tbase(in), id(atom["TL"].id) {}
          inline double operator() (const particle &a, const particle &b, double r2) const {
            double u=Tbase::first(a,b,r2);
            if (a.id==id)
              if (b.id==id)
                u+=Tbase::second(a,b,r2);
            return u;
          }
      };
  }

  /*!
   * \brief Coarse grained three-bead lipid membrane
   * \details This will insert a coarse grained lipid membrane into a current simulation system.
   * The following is handled upon constuction:
   * \li Insertion of lipid bilayed in the XY plane (Cuboid class required)
   * \li Addition of bonded interactions to any existing Hamiltonian (FENE+Harmonic bonds)
   * \li Initialization of Week-Chandler-Andersen parameters
   */
  template<class Tgeometry=Geometry::Cuboid>
    class DesernoMembrane {
      private:
        class _energy : public Energy::Energybase {
          private:
            string _info() {
              using namespace textio;
              std::ostringstream o;
              o << pad(SUB, 10, "Bonds:") << fene.brief()+" - " + harm.brief() + "\n"
                << pad(SUB, 10, "Sigma") << sigma << _angstrom << endl
                << pad(SUB, 10, "Epsilon") << epsilon << kT << endl;
              return o.str();
            }
            Potential::FENE fene;
            Potential::Harmonic harm;
            GroupArray* lipidsPtr; //pointer to group with all lipid molecules
          public:
            double sigma, epsilon;

            PointParticle::Tid idHead, idTail; // particle id for tail and head groups of the lipids

            Tgeometry geometry;

            _energy(
                InputMap &in,
                Potential::WeeksChandlerAndersen &wca,
                Potential::CosAttract &cos,
                GroupArray& lipids) : fene(in), harm(in), lipidsPtr(&lipids), geometry(in) {

              idHead=atom["HD"].id;
              idTail=atom["TL"].id;
              atom["HD"].dp=4.0;
              atom["TL"].dp=3.0;

              name="Deserno 3 bead lipid hamiltonian";
              sigma   = in.get<double>("lipid_sigma", 10);
              epsilon = in.get<double>("lipid_epsilon", 1);
              double headtail_k=0.5*10*epsilon/(sigma*sigma);
              double headtail_req=4*sigma;
              double fene_k=30*epsilon/(sigma*sigma);
              double fene_rmax=1.5*sigma;

              harm = Potential::Harmonic(headtail_k,headtail_req);
              fene = Potential::FENE(fene_k,fene_rmax);
              cos=Potential::CosAttract(in);

              wca.customSigma(idHead, idHead, 0.95*sigma);            
              wca.customSigma(idHead, idTail, 0.95*sigma);            
              wca.customSigma(idTail, idTail, sigma);                 
              wca.customEpsilon(idHead, idHead, epsilon);             
              wca.customEpsilon(idHead, idTail, epsilon);             
              wca.customEpsilon(idTail, idTail, epsilon);   
            }

            double g_internal(const p_vec &p, Group &g) {
              double u=0;
              if (&g==lipidsPtr)  // is g our full molecular array of lipids?
                for (int i=g.front(); i<=g.back(); i+=3)
                  u+=i_internal(p,i);
              else if (lipidsPtr->find(g.front())) // if not, perhaps a single lipid?
                u=i_internal(p, g.front());
              return u;
            }

            double i_internal(const p_vec &p, int i) {
              double u=0;
              if (lipidsPtr->find(i)) {
                int k=i-1;              // i=middle?
                if ( ((i+1)%3)==0) k=i-2; // i=last?
                else if ((i%3)==0) k=i; // i=first
                u+=fene(p[k], p[k+1], geometry.sqdist(p[k],p[k+1]));
                u+=fene(p[k+1], p[k+2], geometry.sqdist(p[k+1],p[k+2]));
                u+=harm(p[k], p[k+2], geometry.sqdist(p[k],p[k+2]));
                assert( p[k].id==idHead );
              }
              return u;
            }
        }; // end of energy class

        _energy* ePtr;

      public:

        GroupArray lipids;  // Group containing all lipid molecules

        DesernoMembrane(
            InputMap &in,
            Energy::Hamiltonian &pot,
            Space &spc,
            Potential::WeeksChandlerAndersen &wca,
            Potential::CosAttract &cos) : lipids(3) {

          // Extend existing Hamiltonian w. bonded interactions
          auto ePtr = pot.create( _energy(in, wca, cos, lipids) );

          // Create a lipid molecule
          p_vec p(3);
          p[0] = atom[ ePtr->idHead ];
          p[0].z=-2.5*ePtr->sigma;
          for (int i=1; i<3; i++) {
            p[i] = atom[ ePtr->idTail ];
            p[i].z = p[i-1].z + ePtr->sigma;
          }

          // Create a bilayer in XY plane and insert info space
          Point u;
          double qfrac=in.get<double>("lipid_chargefraction", 0);          // head group ionization?
          int N=in.get<int>("lipid_N",0);                                  // # of lipids
          while (N-->0) {
            u.x=slp_global.randHalf() * ePtr->geometry.len.x;              // random xy position
            u.y=slp_global.randHalf() * ePtr->geometry.len.y;
            ePtr->geometry.boundary(u);                                    // respect periodicity
            (slp_global.randOne()>0.5) ? u.z=1 : u.z=-1;                   // 50% inverted lipids
            (slp_global.randOne()<qfrac) ? p[0].charge=-1 : p[0].charge=0; // dope w. charge
            for (auto &i : p) {
              i.x=u.x;
              i.y=u.y;
              i.z*=u.z; // flip 50% of the peptides (see above)
            }
            GroupMolecular g = spc.insert(p); // insert lipid into simulation space (no overlap check!)
            lipids.add(g);                    // add inserted lipid to "lipids" group
          }
          spc.enroll( lipids );               // Space needs to know about ALL groups
          cout << "Capacity = " << spc.p.capacity() << endl;
        }
    };
}//namespace
#endif
