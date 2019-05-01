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
    using Base::cnt;
    using Base::coord;
    using Base::f0;
    using Base::file;
    using Base::hisfile;
    using Base::histo;
    using Base::nconv;
    using Base::penalty;
    using Base::quiet;
    using Base::samplings;
    using Base::scale;
    using Base::udelta;

    Eigen::VectorXi weights; // array w. mininum histogram counts
    Eigen::VectorXd buffer;  // receive buffer for penalty functions

    PenaltyMPI(const json &j, Tspace &spc) : Base(j, spc) {
        weights.resize(MPI::mpi.nproc());
        buffer.resize(penalty.size() * MPI::mpi.nproc()); // recieve buffer for penalty func
    }

    void update(const std::vector<double> &c) override {
        using namespace Faunus::MPI;
        double uold = penalty[c];
        if (++cnt % this->nupdate == 0 and f0 > 0) {

            int min = histo.minCoeff(); // if min>0 --> all RC's visited
            MPI_Barrier(mpi.comm);      // wait for all walkers to reach here
            MPI_Allgather(&min, 1, MPI_INT, weights.data(), 1, MPI_INT, mpi.comm);

            // if at least one walker has sampled full RC space at least `samplings` times
            if (weights.maxCoeff() > samplings) { // change to minCoeff()?
                // master collects penalty from all slaves
                MPI_Gather(penalty.data(), penalty.size(), MPI_DOUBLE, buffer.data(), penalty.size(), MPI_DOUBLE, 0,
                           mpi.comm);

                // master performs the average
                if (mpi.isMaster()) {
                    penalty.setZero();
                    for (int i = 0; i < mpi.nproc(); i++)
                        penalty += Eigen::Map<Eigen::MatrixXd>(buffer.data() + i * penalty.size(), penalty.rows(),
                                                               penalty.cols());
                    penalty = (penalty.array() - penalty.minCoeff()) / double(mpi.nproc());
                }

                // master sends the averaged penalty function to all slaves
                MPI_Bcast(penalty.data(), penalty.size(), MPI_DOUBLE, 0, mpi.comm);
                nconv += 1;

                // at this point, *all* penalty functions shall be identical

                // save penalty function to disk
                if (mpi.isMaster()) {
                    std::ofstream f(file + ".walkersync" + std::to_string(nconv));
                    if (f) {
                        f.precision(16);
                        f << "# " << f0 << " " << samplings << " " << nconv << "\n" << penalty.array() << endl;
                    }
                }

                // save histogram to disk
                std::ofstream f(MPI::prefix + hisfile + ".walkersync" + std::to_string(nconv));
                if (f) {
                    f << histo << endl;
                    f.close();
                }

                // print information to console
                if (min > 0 and not quiet) {
                    cout << "Barriers/kT: penalty = " << penalty.maxCoeff()
                         << " histogram = " << std::log(double(histo.maxCoeff()) / histo.minCoeff()) << endl;
                }

                histo.setZero();
                f0 = f0 * scale; // reduce penalty energy
                samplings = std::ceil(samplings / scale);
            }
        }
        coord = c;
        histo[coord]++;
        penalty[coord] += f0;
        udelta += penalty[coord] - uold;
    } //!< Average penalty function across all nodes
};    //!< Penalty function with MPI exchange
#endif

} // end of Energy namespace
} // end of Faunus namespace
