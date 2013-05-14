#include <faunus/common.h>
#include <faunus/species.h>
#include <faunus/inputfile.h>
#include <faunus/space.h>
#include <faunus/energy.h>
#include <faunus/titrate.h>
#include <faunus/json.h>

namespace Faunus {

  namespace Energy {

 
 
  }//namespace Energy

  namespace Move {

    /**
     * @brief As SwapMove but Minimizes Short Ranged interactions
     *        within a molecule upon swapping
     *
     * Before calculating dU of an attempted swap move, radii on
     * particles within the SAME group are set to minus radius of
     * the swapped particle and hydrophobicity is set to false.
     * This to minimize large interactions in molecules with overlapping
     * particles - i.e LJ will be zero. It can also be used to avoid
     * internal hydrophobic interactions in rigid groups upon swapping
     * between hydrophobic and non-hydrophobic species.
     * Group information is found in Space::g and to avoid energy drifts by
     * ignoring hydrophobic interactions internally in groups the
     * Energy::EnergyRest is used to collect the missing contribution to dU.
     */
    class SwapMoveMSR : public SwapMove {
    private:
      std::map<int, double> radiusbak;    // backup for radii
      std::map<int, bool> hydrophobicbak; // backup for hydrophobic state
      Energy::EnergyRest potrest; // dummy energy class missing contributions to dU in energy drift calc.
      double _energyChange();
      void modify();
      void restore();
    public:
      SwapMoveMSR(InputMap&, Energy::Hamiltonian&, Space&, string="swapmv_");
    };


    SwapMoveMSR::SwapMoveMSR(
        InputMap &in, Energy::Hamiltonian &ham, Space &spc, string pfx) : SwapMove(in,ham,spc,pfx) {
      title+=" (min. shortrange)";
      ham.add( potrest );
    }

    void SwapMoveMSR::modify() {
      radiusbak.clear();
      hydrophobicbak.clear();
      for (auto g : spc->groupList() )   // loop over all groups
        if (g->find(ipart)) {  //   is ipart part of a group?
          for (auto i : *g)    //     if so, loop over that group
            if (i!=ipart) {    //       and ignore ipart
              assert( abs(spc->p[i].radius-spc->trial[i].radius)<1e-9);
              assert( spc->p[i].hydrophobic==spc->trial[i].hydrophobic);

              //radiusbak[i]         = spc->p[i].radius;
              //spc->p[i].radius     = -spc->p[ipart].radius;
              //spc->trial[i].radius = -spc->p[ipart].radius;

              hydrophobicbak[i]         = spc->p[i].hydrophobic;
              spc->p[i].hydrophobic     = false;
              spc->trial[i].hydrophobic = false;
            }
          return; // a particle can be part of a single group, only
        }
    }

    void SwapMoveMSR::restore() {
      for (auto &m : radiusbak) {
        spc->p[m.first].radius = m.second;
        spc->trial[m.first].radius = m.second;
      }
      for (auto &m : hydrophobicbak) {
        spc->p[m.first].hydrophobic = m.second;
        spc->trial[m.first].hydrophobic = m.second;
      }
    }

    double SwapMoveMSR::_energyChange() {
      double du_orig = SwapMove::_energyChange();
      modify();
      double du = SwapMove::_energyChange();
      restore();
      potrest.add( du-du_orig );
      return du;
    }

  }//Move namespace
}//Faunus namespace

