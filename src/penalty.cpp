#include "penalty.h"

namespace Faunus {
namespace Energy {

Penalty::Penalty(const json &j, Tspace &spc) : spc(spc) {
    using namespace ReactionCoordinate;
    name = "penalty";
    overwrite_penalty = j.value("overwrite", true);
    f0 = j.at("f0").get<double>();
    scale = j.at("scale").get<double>();
    quiet = j.value("quiet", true);
    nupdate = j.at("update").get<double>();
    samplings = j.value("samplings", 1);
    nodrift = j.value("nodrift", true);
    file = j.at("file").get<std::string>();
    hisfile = j.value("histogram", "penalty-histogram.dat");
    std::vector<double> binwidth, min, max;

    if (scale < 0 or scale > 1)
        throw std::runtime_error("`scale` must be in the interval [0:1]");

    for (auto &i : j.at("coords")) {
        if (i.is_object() and i.size() == 1) {
            auto rc = ReactionCoordinate::createReactionCoordinate(i, spc); // throws if unknown rc given
            assert(rc);                                                     // after this rc *must* be defined
            if (rc->min >= rc->max || rc->binwidth <= 0)
                throw std::runtime_error("min<max and binwidth>0 required for penalty reaction coordinate");
            rcvec.push_back(rc);
            binwidth.push_back(rc->binwidth);
            min.push_back(rc->min);
            max.push_back(rc->max);
        }
    }
    dim = binwidth.size();
    if (dim < 1 or dim > 2)
        throw std::runtime_error("exactly one or two coordinates required");

    coord.resize(rcvec.size(), 0);
    histo.reInitializer(binwidth, min, max);
    penalty.reInitializer(binwidth, min, max);

    std::ifstream f(MPI::prefix + file);
    if (f) {
        cout << "Loading penalty function '" << MPI::prefix + file << "'" << endl;
        std::string hash;
        f >> hash >> f0 >> samplings >> nconv;
        for (int row = 0; row < penalty.rows(); row++)
            for (int col = 0; col < penalty.cols(); col++)
                if (not f.eof())
                    f >> penalty(row, col);
                else
                    throw std::runtime_error("penalty file dimension mismatch");
    }
}
Penalty::~Penalty() {
    if (overwrite_penalty) {
        std::ofstream f(MPI::prefix + file);
        if (f) {
            f.precision(16);
            f << "# " << f0 << " " << samplings << " " << nconv << "\n" << penalty.array() - penalty.minCoeff() << "\n";
            f.close();
        }
    }

    std::ofstream f2(MPI::prefix + hisfile);
    if (f2)
        f2 << histo << "\n";
    // add function to save to numpy-friendly file...
}
void Penalty::to_json(json &j) const {
    j["file"] = file;
    j["scale"] = scale;
    j["update"] = nupdate;
    j["nodrift"] = nodrift;
    j["histogram"] = hisfile;
    j["f0_final"] = f0;
    j["overwrite"] = overwrite_penalty;
    auto &_j = j["coords"] = json::array();
    for (auto rc : rcvec)
        _j.push_back(*rc); // `ReactionCoordinateBase` --> `json`
}
double Penalty::energy(Change &change) {
    assert(rcvec.size() <= coord.size());
    double u = 0;
    coord.resize(rcvec.size());
    if (change) {
        for (size_t i = 0; i < rcvec.size(); i++) {
            coord.at(i) = rcvec[i]->operator()();
            if (not rcvec[i]->inRange(coord[i]))
                return pc::infty;
        }
        penalty.to_index(coord);
        u = penalty[coord];
    }
    // reaching here, `coord` always reflects
    // the current reaction coordinate
    return (nodrift) ? u - udelta : u;
}
void Penalty::update(const std::vector<double> &c) {
    if (++cnt % nupdate == 0 and f0 > 0) {
        bool b = histo.minCoeff() >= (int)samplings;
        if (b) {
            double min = penalty.minCoeff(); // define minimun penalty energy
            penalty = penalty.array() - min; // ...to zero
            if (not quiet)
                cout << "Barriers/kT: penalty = " << penalty.maxCoeff()
                     << " histogram = " << std::log(double(histo.maxCoeff()) / histo.minCoeff()) << endl;
            f0 = f0 * scale; // reduce penalty energy
            samplings = std::ceil(samplings / scale);
            histo.setZero();
            udelta += -min;
        }
    }
    coord = c;
    histo[coord]++;
    penalty[coord] += f0;
    udelta += f0;
}
void Penalty::sync(Energybase *basePtr, Change &) {
    // this function is called when a move is accepted
    // or rejected, as well as when initializing the system
    auto other = dynamic_cast<decltype(this)>(basePtr);
    assert(other);
    update(other->coord);
    other->update(other->coord); // this is to keep cnt and samplings in sync

    // some assertions...
    assert(samplings == other->samplings);
    assert(coord == other->coord);
    assert(cnt == other->cnt);
    assert(udelta == other->udelta);
}

#ifdef ENABLE_MPI

PenaltyMPI::PenaltyMPI(const json &j, Tspace &spc) : Penalty(j, spc) {
    weights.resize(MPI::mpi.nproc());
    buffer.resize(penalty.size() * MPI::mpi.nproc()); // recieve buffer for penalty func
}

void PenaltyMPI::update(const std::vector<double> &c) {
    using namespace Faunus::MPI;
    double uold = penalty[c];
    if (++cnt % this->nupdate == 0 and f0 > 0) {

        int min = histo.minCoeff(); // if min>0 --> all RC's visited
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
                    penalty +=
                        Eigen::Map<Eigen::MatrixXd>(buffer.data() + i * penalty.size(), penalty.rows(), penalty.cols());
                penalty = (penalty.array() - penalty.minCoeff()) / double(mpi.nproc());
            }

            // master sends the averaged penalty function to all slaves
            MPI_Bcast(penalty.data(), penalty.size(), MPI_DOUBLE, 0, mpi.comm);
            nconv += 1;

            // at this point, *all* penalty functions shall be identical

            /* save penalty function to disk
            if (mpi.isMaster()) {
                std::ofstream f(file + ".walkersync" + std::to_string(nconv));
                if (f) {
                    f.precision(16);
                    f << "# " << f0 << " " << samplings << " " << nconv << "\n" << penalty.array() << endl;
                }
            }*/

            /* save histogram to disk
            std::ofstream f(MPI::prefix + hisfile + ".walkersync" + std::to_string(nconv));
            if (f) {
                f << histo << endl;
                f.close();
            }*/

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
}

#endif
} // namespace Energy
} // namespace Faunus
