#pragma once

#include "mpi.h"
#include "energy.h"
#include "reactioncoordinate.h"

namespace Faunus {
namespace Energy {

/**
 * `udelta` is the total change of updating the energy function. If
 * not handled this will appear as an energy drift (which it is!). To
 * avoid this, this term is added to the energy but since it's the
 * same in both the trial and old state energies it will not affect
 * MC move acceptance.
 */
class Penalty : public Energybase {
  protected:
    typedef typename Tspace::Tgroup Tgroup;
    typedef typename Tspace::Tpvec Tpvec;
    typedef typename std::shared_ptr<ReactionCoordinate::ReactionCoordinateBase> Tcoord;

    Tspace &spc;
    bool overwrite_penalty = true; // overwrites the input penalty function
    bool nodrift;                  // avoid energy drift when upgrading penalty function
    bool quiet;                    // hold mund
    size_t dim = 0;                // number of reaction coordinate
    size_t cnt = 0;                // number of calls to `sync()`
    size_t nupdate;                // update frequency [steps]
    size_t samplings;
    size_t nconv = 0;
    double udelta = 0; // total energy change of updating penalty function
    double scale;      // scaling factor for f0
    double f0;         // penalty increment
    std::string file, hisfile;
    std::vector<Tcoord> rcvec; // vector of reaction coordinate functions (length = 1 or 2)
    std::vector<double> coord; // latest reaction coordinate (length = 1 or 2)

    Table<int> histo;      // sampling along reaction coordinates
    Table<double> penalty; // penalty function

  public:
    Penalty(const json &j, Tspace &spc);
    virtual ~Penalty();
    void to_json(json &j) const override;
    double energy(Change &change) override;

    /*
     * @todo: If this is called before `energy()`, the coord
     * is never calculated and causes undefined behavior
     */
    virtual void update(const std::vector<double> &c);

    void sync(Energybase *basePtr, Change &) override; // @todo: this doubles the MPI communication
};

#ifdef ENABLE_MPI
struct PenaltyMPI : public Penalty {
    Eigen::VectorXi weights; // array w. mininum histogram counts
    Eigen::VectorXd buffer;  // receive buffer for penalty functions

    PenaltyMPI(const json &j, Tspace &spc);
    void update(const std::vector<double> &c) override; //!< Average penalty function across all nodes
};    //!< Penalty function with MPI exchange
#endif

} // end of Energy namespace
} // end of Faunus namespace
