#pragma once

#include "core.h"
#include "energy.h"
#include "average.h"
#include "potentials.h"
#include "mpi.h"

namespace Faunus {
namespace Move {

class Movebase {
  private:
    virtual void _move(Change &) = 0;                          //!< Perform move and modify change object
    virtual void _accept(Change &);                            //!< Call after move is accepted
    virtual void _reject(Change &);                            //!< Call after move is rejected
    virtual void _to_json(json &j) const = 0;                  //!< Extra info for report if needed
    virtual void _from_json(const json &j) = 0;                //!< Extra info for report if needed
    TimeRelativeOfTotal<std::chrono::microseconds> timer;      //!< Timer for whole move
    TimeRelativeOfTotal<std::chrono::microseconds> timer_move; //!< Timer for _move() only
  protected:
    unsigned long cnt = 0;
    unsigned long accepted = 0;
    unsigned long rejected = 0;

  public:
    static Random slump; //!< Shared for all moves
    std::string name;    //!< Name of move
    std::string cite;    //!< Reference
    int repeat = 1;      //!< How many times the move should be repeated per sweep

    void from_json(const json &j);
    void to_json(json &j) const; //!< JSON report w. statistics, output etc.
    void move(Change &change);   //!< Perform move and modify given change object
    void accept(Change &c);
    void reject(Change &c);
    virtual double bias(Change &, double uold,
                        double unew); //!< adds extra energy change not captured by the Hamiltonian
    inline virtual ~Movebase() = default;
};

void from_json(const json &j, Movebase &m); //!< Configure any move via json
void to_json(json &j, const Movebase &m);

/**
 * @brief Swap the charge of a single atom
 */
template <typename Tspace> class AtomicSwapCharge : public Movebase {
  private:
    typedef typename Tspace::Tpvec Tpvec;
    typedef typename Tspace::Tparticle Tparticle;
    Tspace &spc; // Space to operate on
    int molid = -1;
    double ln10 = log(10);
    double pKa, pH;
    Average<double> msqd; // mean squared displacement
    double _sqd, _bias;   // squared displament
    std::string molname;  // name of molecule to operate on
    Change::data cdata;

    void _to_json(json &j) const override {
        j = {{"pH", pH},
             {"pka", pKa},
             {"molid", molid},
             {u8::rootof + u8::bracket("r" + u8::squared), std::sqrt(msqd.avg())},
             {"molecule", molname}};
        _roundjson(j, 3);
    }

    void _from_json(const json &j) override {
        assert(!molecules.empty());
        try {
            molname = j.at("molecule");
            auto it = findName(molecules, molname);
            if (it == molecules.end())
                throw std::runtime_error("unknown molecule '" + molname + "'");
            molid = it->id();
            pH = j.at("pH").get<double>();
            pKa = j.at("pKa").get<double>();
            if (repeat < 0) {
                auto v = spc.findMolecules(molid);
                repeat = std::distance(v.begin(), v.end()); // repeat for each molecule...
                if (repeat > 0)
                    repeat = repeat * v.front().size(); // ...and for each atom
            }
        } catch (std::exception &e) {
            std::cerr << name << ": " << e.what();
            throw;
        }
    } //!< Configure via json object

    typename Tpvec::iterator randomAtom() {
        assert(molid >= 0);
        auto mollist = spc.findMolecules(molid); // all `molid` groups
        if (size(mollist) > 0) {
            auto git = slump.sample(mollist.begin(), mollist.end()); // random molecule iterator
            if (!git->empty()) {
                auto p = slump.sample(git->begin(), git->end());         // random particle iterator
                cdata.index = Faunus::distance(spc.groups.begin(), git); // integer *index* of moved group
                cdata.atoms[0] = std::distance(git->begin(), p);         // index of particle rel. to group
                return p;
            }
        }
        return spc.p.end();
    }

    void _move(Change &change) override {
        auto p = randomAtom();
        if (p != spc.p.end()) {
            auto &g = spc.groups[cdata.index];
            double oldcharge = p->charge;
            p->charge = fabs(oldcharge - 1);
            _sqd = fabs(oldcharge - 1) - oldcharge;
            change.groups.push_back(cdata);   // add to list of moved groups
            _bias = _sqd * (pH - pKa) * ln10; // one may add bias here...
        }
    }

    double bias(Change &, double, double) override {
        return _bias;
    } //!< adds extra energy change not captured by the Hamiltonian

    void _accept(Change &) override { msqd += _sqd; }
    void _reject(Change &) override { msqd += 0; }

  public:
    AtomicSwapCharge(Tspace &spc) : spc(spc) {
        name = "swapcharge";
        repeat = -1; // meaning repeat N times
        cdata.atoms.resize(1);
        cdata.internal = true;
    }
};

/**
 * @brief Translate and rotate a molecular group
 */
class AtomicTranslateRotate : public Movebase {
  protected:
    typedef typename Tspace::Tpvec Tpvec;
    Tspace &spc; // Space to operate on
    int molid = -1;
    Point dir = {1, 1, 1};
    Average<double> msqd; // mean squared displacement
    double _sqd;          // squared displament
    std::string molname;  // name of molecule to operate on
    Change::data cdata;

    void _to_json(json &j) const override;
    void _from_json(const json &j) override; //!< Configure via json object
    std::vector<Particle>::iterator randomAtom();

    /**
     * @brief translates a single particle.
     */
    virtual void translateParticle(typename Tpvec::iterator p, double dp);
    void _move(Change &change) override;
    void _accept(Change &) override;
    void _reject(Change &) override;

  public:
    AtomicTranslateRotate(Tspace &spc);
};

/**
 * @brief Translate and rotate an atom on a 2D hypersphere-surface
 * @todo under construction
 */
/*
template<typename Tspace>
   class Atomic2dTranslateRotate : public AtomicTranslateRotate {
       protected:
           typedef AtomicTranslateRotate base;
           using base::spc;

           void translateParticle(typename base::Tpvec::iterator p, double dp) override {
               auto &g = spc.groups[base::cdata.index];
               Point oldpos = p->pos;

               Point rtp = xyz2rtp(p->pos); // Get the spherical coordinates of the particle
               double slump_theta = dp * (base::slump() - 0.5); // Get random theta-move
               double slump_phi = dp * (base::slump() - 0.5);   // Get random phi-move

               double scalefactor_theta = spc.geo.getRadius() * sin(rtp.z()); // Scale-factor for theta
               double scalefactor_phi = spc.geo.getRadius();                  // Scale-factor for phi

               Point theta_dir = Point(-sin(rtp.y()), cos(rtp.y()), 0); // Unit-vector in theta-direction
               Point phi_dir = Point(cos(rtp.y()) * cos(rtp.z()), sin(rtp.y()) * cos(rtp.z()),
                                     -sin(rtp.z())); // Unit-vector in phi-direction
               Point xyz = oldpos + scalefactor_theta * theta_dir * slump_theta +
                           scalefactor_phi * phi_dir * slump_phi; // New position
               p->pos = spc.geo.getRadius() * xyz / xyz.norm();   // Convert to cartesian coordinates

               spc.geo.boundary(p->pos);
               base::_sqd = spc.geo.sqdist(oldpos, p->pos); // squared displacement
               if (not g.atomic) {                          // recalc mass-center for non-molecular groups
                   g.cm = Geometry::massCenter(g.begin(), g.end(), spc.geo.getBoundaryFunc(), -g.cm);
#ifndef NDEBUG
                   Point cmbak = g.cm;                             // backup mass center
                   g.translate(-cmbak, spc.geo.getBoundaryFunc()); // translate to {0,0,0}
                   double should_be_zero = spc.geo.sqdist({0, 0, 0}, Geometry::massCenter(g.begin(),
g.end())); if (should_be_zero > 1e-6) throw std::runtime_error("atomic move too large"); else g.translate(cmbak,
spc.geo.getBoundaryFunc()); #endif
               }
           }

       public:
           Atomic2dTranslateRotate(Tspace &spc) : base(spc) {
               base::name = "transrot 2d";
           }
   };*/

/**
 * @brief Translate and rotate a molecular group
 */
template <typename Tspace> class TranslateRotate : public Movebase {
  protected:
    typedef typename Tspace::Tpvec Tpvec;
    Tspace &spc; // Space to operate on
    int molid = -1;
    double dptrans = 0;
    double dprot = 0;
    Point dir = {1, 1, 1};
    double _sqd;          // squared displacement
    Average<double> msqd; // mean squared displacement

    void _to_json(json &j) const override {
        j = {{"dir", dir},
             {"dp", dptrans},
             {"dprot", dprot},
             {"molid", molid},
             {u8::rootof + u8::bracket("r" + u8::squared), std::sqrt(msqd.avg())},
             {"molecule", molecules[molid].name}};
        _roundjson(j, 3);
    }

    void _from_json(const json &j) override {
        assert(!molecules.empty());
        try {
            std::string molname = j.at("molecule");
            auto it = findName(molecules, molname);
            if (it == molecules.end())
                throw std::runtime_error("unknown molecule '" + molname + "'");
            molid = it->id();
            dir = j.value("dir", Point(1, 1, 1));
            dprot = j.at("dprot");
            dptrans = j.at("dp");
            if (repeat < 0) {
                auto v = spc.findMolecules(molid);
                repeat = std::distance(v.begin(), v.end());
            }
        } catch (std::exception &e) {
            throw std::runtime_error(name + ": " + e.what());
        }
    } //!< Configure via json object

    void _move(Change &change) override {
        assert(molid >= 0);
        assert(!spc.groups.empty());
        assert(spc.geo.getVolume() > 0);

        // pick random group from the system matching molecule type
        // TODO: This can be slow -- implement look-up-table in Space
        auto mollist = spc.findMolecules(molid, Tspace::ACTIVE); // list of molecules w. 'molid'
        if (size(mollist) > 0) {
            auto it = slump.sample(mollist.begin(), mollist.end());
            if (not it->empty()) {
                assert(it->id == molid);

                if (dptrans > 0) { // translate
                    Point oldcm = it->cm;
                    Point dp = ranunit(slump, dir) * dptrans * slump();

                    it->translate(dp, spc.geo.getBoundaryFunc());
                    _sqd = spc.geo.sqdist(oldcm, it->cm); // squared displacement
                }

                if (dprot > 0) { // rotate
                    Point u = ranunit(slump);
                    double angle = dprot * (slump() - 0.5);
                    Eigen::Quaterniond Q(Eigen::AngleAxisd(angle, u));
                    it->rotate(Q, spc.geo.getBoundaryFunc());
                }

                if (dptrans > 0 || dprot > 0) { // define changes
                    Change::data d;
                    d.index = Faunus::distance(spc.groups.begin(), it); // integer *index* of moved group
                    d.all = true;                                       // *all* atoms in group were moved
                    change.groups.push_back(d);                         // add to list of moved groups
                }
                assert(spc.geo.sqdist(it->cm, Geometry::massCenter(it->begin(), it->end(), spc.geo.getBoundaryFunc(),
                                                                   -it->cm)) < 1e-6);
            }
        }
    }

    void _accept(Change &) override { msqd += _sqd; }
    void _reject(Change &) override { msqd += 0; }

  public:
    TranslateRotate(Tspace &spc) : spc(spc) {
        name = "moltransrot";
        repeat = -1; // meaning repeat N times
    }
};

/**
 * @brief Move that will swap conformation of a molecule
 *
 * This will swap between different molecular conformations
 * as defined in `MoleculeData` with `traj` and `weight`.
 * If defined, the weight
 * distribution is respected, otherwise all conformations
 * have equal intrinsic weight. Upon insertion, the new conformation
 * is randomly oriented and placed on top of the mass-center of
 * an exising molecule. That is, there is no mass center movement.
 *
 * @todo Add feature to align molecule on top of an exiting one
 * @todo Expand `_info()` to show number of conformations
 * @warning Weighted distributions untested and not verified for correctness
 * @date Malmo, November 2016
 */
template <class Tspace> class ConformationSwap : public Movebase {
  private:
    typedef typename Tspace::Tpvec Tpvec;
    typedef MoleculeData Tmoldata;
    RandomInserter inserter;
    Tspace &spc; // Space to operate on
    int molid = -1;
    int newconfid = -1;

    void _to_json(json &j) const override {
        j = {{"molid", molid}, {"molecule", molecules[molid].name}};
        _roundjson(j, 3);
    }

    void _from_json(const json &j) override {
        assert(!molecules.empty());
        try {
            std::string molname = j.at("molecule");
            auto it = findName(molecules, molname);
            if (it == molecules.end())
                throw std::runtime_error("unknown molecule '" + molname + "'");
            molid = it->id();
            if (molecules[molid].conformations.size() < 2)
                throw std::runtime_error("minimum two conformations required");
            if (repeat < 0) {
                auto v = spc.findMolecules(molid);
                repeat = std::distance(v.begin(), v.end());
            }
        } catch (std::exception &e) {
            throw std::runtime_error(name + ": " + e.what());
        }
    } //!< Configure via json object

    void _move(Change &change) override {
        assert(molid >= 0);
        assert(change.empty());

        auto mollist = spc.findMolecules(molid, Tspace::ACTIVE); // list of molecules w. 'molid'
        if (size(mollist) > 0) {
            auto g = slump.sample(mollist.begin(), mollist.end());
            if (not g->empty()) {
                inserter.offset = g->cm;

                // Get a new conformation that should be properly wrapped around the boundaries
                // (if applicable) and have the same mass-center as "g->cm".
                Tpvec p = inserter(spc.geo, spc.p, molecules[molid]);
                if (p.size() not_eq g->size())
                    throw std::runtime_error(name + ": conformation atom count mismatch");

                newconfid = molecules[molid].conformations.index;

                std::copy(p.begin(), p.end(), g->begin()); // override w. new conformation
#ifndef NDEBUG
                // this move shouldn't move mass centers, so let's check if this is true:
                Point newcm = Geometry::massCenter(p.begin(), p.end(), spc.geo.getBoundaryFunc(), -g->cm);
                if ((newcm - g->cm).norm() > 1e-6)
                    throw std::runtime_error(name + ": unexpected mass center movement");
#endif
                Change::data d;
                d.index = Faunus::distance(spc.groups.begin(), g); // integer *index* of moved group
                d.all = true;                                      // *all* atoms in group were moved
                d.internal = false;                                // we *don't* want to calculate the internal energy
                change.groups.push_back(d);                        // add to list of moved groups
            }
        }
    }

    void _accept(Change &change) override {
        assert(change.groups.size() == 1);
        spc.groups[change.groups.front().index].confid = newconfid;
    }

  public:
    ConformationSwap(Tspace &spc) : spc(spc) {
        name = "conformationswap";
        repeat = -1; // meaning repeat n times
        inserter.dir = {0, 0, 0};
        inserter.rotate = true;
        inserter.allowoverlap = true;
    }

}; // end of conformation swap move

/**
 * @brief Sketch for MD move
 */
template <typename Tspace> class ForceMove : public Movebase {
  private:
    typedef typename Tspace::Tpvec Tpvec;
    void _to_json(json &) const override{};
    void _from_json(const json &) override{};
    std::vector<Point> forces, velocities;

  public:
    ForceMove() {
        // resize forces and velocities to mathc spc.p
    }
}; // end of forcemove

#ifdef DOCTEST_LIBRARY_INCLUDED
TEST_CASE("[Faunus] TranslateRotate") {
    typedef typename Space::Tpvec Tpvec;

    CHECK(!atoms.empty());     // set in a previous test
    CHECK(!molecules.empty()); // set in a previous test

    Tspace spc;
    TranslateRotate<Tspace> mv(spc);
    json j = R"( {"molecule":"B", "dp":1.0, "dprot":0.5, "dir":[0,1,0], "repeat":2 })"_json;
    mv.from_json(j);

    j = json(mv).at(mv.name);
    CHECK(j.at("molecule") == "B");
    CHECK(j.at("dir") == Point(0, 1, 0));
    CHECK(j.at("dp") == 1.0);
    CHECK(j.at("repeat") == 2);
    CHECK(j.at("dprot") == 0.5);
}
#endif

template <typename Tspace> class VolumeMove : public Movebase {
  private:
    const std::map<std::string, Geometry::VolumeMethod> methods = {
        {"xy", Geometry::XY}, {"isotropic", Geometry::ISOTROPIC}, {"isochoric", Geometry::ISOCHORIC}};
    typename decltype(methods)::const_iterator method;
    typedef typename Tspace::Tpvec Tpvec;
    Tspace &spc;
    Average<double> msqd, Vavg; // mean squared displacement
    double dV = 0, deltaV = 0, Vnew = 0, Vold = 0;

    void _to_json(json &j) const override {
        using namespace u8;
        if (cnt > 0) {
            j = {{"dV", dV},
                 {"method", method->first},
                 {bracket("V"), Vavg.avg()},
                 {rootof + bracket(Delta + "V" + squared), std::sqrt(msqd.avg())},
                 {cuberoot + rootof + bracket(Delta + "V" + squared), std::cbrt(std::sqrt(msqd.avg()))}};
            _roundjson(j, 3);
        }
    }

    void _from_json(const json &j) override {
        try {
            method = methods.find(j.value("method", "isotropic"));
            if (method == methods.end())
                std::runtime_error("unknown volume change method");
            dV = j.at("dV");
        } catch (std::exception &e) {
            throw std::runtime_error(e.what());
        }
    }

    void _move(Change &change) override {
        if (dV > 0) {
            change.dV = true;
            change.all = true;
            Vold = spc.geo.getVolume();
            Vnew = std::exp(std::log(Vold) + (slump() - 0.5) * dV);
            deltaV = Vnew - Vold;
            spc.scaleVolume(Vnew, method->second);
        } else
            deltaV = 0;
    }

    void _accept(Change &) override {
        msqd += deltaV * deltaV;
        Vavg += spc.geo.getVolume();
    }
    void _reject(Change &) override {
        msqd += 0;
        Vavg += spc.geo.getVolume();
    }

  public:
    VolumeMove(Tspace &spc) : spc(spc) {
        name = "volume";
        repeat = 1;
    }
}; // end of VolumeMove

/**
 * @brief Displaces charge on a single atom
 */
template <typename Tspace> class ChargeMove : public Movebase {
  private:
    typedef typename Tspace::Tpvec Tpvec;
    Tspace &spc;          // Space to operate on
    Average<double> msqd; // mean squared displacement
    double dq = 0, deltaq = 0;
    int atomIndex;
    Change::data cdata;

    void _to_json(json &j) const override {
        using namespace u8;
        j = {{"index", atomIndex},
             {"dq", dq},
             {rootof + bracket(Delta + "q" + squared), std::sqrt(msqd.avg())},
             {cuberoot + rootof + bracket(Delta + "q" + squared), std::cbrt(std::sqrt(msqd.avg()))}};
        _roundjson(j, 3);
    }

    void _from_json(const json &j) override {
        dq = j.at("dq").get<double>();
        atomIndex = j.at("index").get<int>();
        auto git = spc.findGroupContaining(spc.p[atomIndex]);                    // group containing atomIndex
        cdata.index = std::distance(spc.groups.begin(), git);                    // integer *index* of moved group
        cdata.atoms[0] = std::distance(git->begin(), spc.p.begin() + atomIndex); // index of particle rel. to group
    }

    void _move(Change &change) override {
        if (dq > 0) {
            auto &p = spc.p[atomIndex]; // refence to particle
            double qold = p.charge;
            p.charge += dq * (slump() - 0.5);
            deltaq = p.charge - qold;
            change.groups.push_back(cdata); // add to list of moved groups
        } else
            deltaq = 0;
    }

    void _accept(Change &) override { msqd += deltaq * deltaq; }
    void _reject(Change &) override { msqd += 0; }

  public:
    ChargeMove(Tspace &spc) : spc(spc) {
        name = "charge";
        repeat = 1;
        cdata.internal = true; // the group is internally changed
        cdata.atoms.resize(1); // we change exactly one atom
    }
};

/**
 * @brief QuadrantJump translates a molecule to another quadrant
 * considering as the origin the center of the box or the center of mass
 * of a range of atomic indexes specified by "index": [start:stop].
 */
template <typename Tspace> class QuadrantJump : public Movebase {
  private:
    typedef typename Tspace::Tpvec Tpvec;
    typedef typename Tspace::Tparticle Tparticle;
    Tspace &spc; // Space to operate on
    int molid = -1;
    Point dir = {1, 1, 1};
    std::vector<size_t> index;
    double _sqd;          // squared displacement
    Average<double> msqd; // mean squared displacement

    void _to_json(json &j) const override {
        j = {{"dir", dir},
             {"molid", molid},
             {u8::rootof + u8::bracket("r" + u8::squared), std::sqrt(msqd.avg())},
             {"molecule", molecules[molid].name}};
        _roundjson(j, 3);
    }

    void _from_json(const json &j) override {
        assert(!molecules.empty());
        try {
            std::string molname = j.at("molecule");
            auto it = findName(molecules, molname);
            if (it == molecules.end())
                throw std::runtime_error("unknown molecule '" + molname + "'");
            molid = it->id();
            dir = j.value("dir", Point(1, 1, 1));
            index = j.value("index", decltype(index)());
            if (repeat < 0) {
                auto v = spc.findMolecules(molid);
                repeat = std::distance(v.begin(), v.end());
            }
        } catch (std::exception &e) {
            throw std::runtime_error(name + ": " + e.what());
        }
    } //!< Configure via json object

    void _move(Change &change) override {
        assert(molid >= 0);
        assert(!spc.groups.empty());
        assert(spc.geo.getVolume() > 0);

        // pick random group from the system matching molecule type
        // TODO: This can be slow -- implement look-up-table in Space
        auto mollist = spc.findMolecules(molid, Tspace::ACTIVE); // list of molecules w. 'molid'
        if (size(mollist) > 0) {
            auto it = slump.sample(mollist.begin(), mollist.end());
            if (not it->empty()) {
                assert(it->id == molid);
                Point oldcm = it->cm;
                if (index.size() == 2) {
                    auto cm_O = Geometry::massCenter(spc.p.begin() + index[0], spc.p.begin() + index[1] + 1,
                                                     spc.geo.getBoundaryFunc());
                    it->translate(-2 * spc.geo.vdist(oldcm, cm_O).cwiseProduct(dir.cast<double>()),
                                  spc.geo.getBoundaryFunc());
                } else {
                    it->translate(-2 * oldcm.cwiseProduct(dir.cast<double>()), spc.geo.getBoundaryFunc());
                }
                _sqd = spc.geo.sqdist(oldcm, it->cm); // squared displacement
                Change::data d;
                d.index = Faunus::distance(spc.groups.begin(), it); // integer *index* of moved group
                d.all = true;                                       // *all* atoms in group were moved
                change.groups.push_back(d);                         // add to list of moved groups

                assert(spc.geo.sqdist(it->cm, Geometry::massCenter(it->begin(), it->end(), spc.geo.getBoundaryFunc(),
                                                                   -it->cm)) < 1e-9);
            }
        } else
            std::cerr << name << ": no molecules found" << std::endl;
    }

    void _accept(Change &) override { msqd += _sqd; }
    void _reject(Change &) override { msqd += 0; }

  public:
    QuadrantJump(Tspace &spc) : spc(spc) {
        name = "quadrantjump";
        repeat = -1; // meaning repeat N times
    }
};

/**
 * @brief Molecular cluster move
 *
 * @todo fix so that it works w. GC molecules (index are calculating before simulation)
 */
template <typename Tspace> class Cluster : public Movebase {
  private:
    typedef typename Tspace::Tpvec Tpvec;
    typedef typename Tspace::Tgroup Tgroup;
    Tspace &spc;
    Average<double> msqd, msqd_angle, N;
    double dptrans = 0, dprot = 0, angle = 0, _bias = 0;
    size_t bias_rejected = 0;
    bool rotate; // true if cluster should be rotated
    Point dir = {1, 1, 1}, dp;
    std::vector<std::string> names; // names of molecules to be considered
    std::vector<int> ids;           // molecule id's of molecules to be considered
    std::set<int> satellites;       // subset of molecules to cluster, but NOT act as nuclei (cluster centers)
    std::vector<size_t> index;      // index of all possible molecules to be considered
    std::map<size_t, size_t> clusterSizeDistribution; // distribution of cluster sizes
    PairMatrix<double, true> thresholdsq;

    virtual double clusterProbability(const Tgroup &g1, const Tgroup &g2) const {
        if (spc.geo.sqdist(g1.cm, g2.cm) <= thresholdsq(g1.id, g2.id))
            return 1.0;
        return 0.0;
    }

    void _to_json(json &j) const override {
        using namespace u8;
        j = {{"dir", dir},
             {"dp", dptrans},
             {"dprot", dprot},
             {rootof + bracket("r" + squared), std::sqrt(msqd.avg())},
             {rootof + bracket(theta + squared) + "/" + degrees, std::sqrt(msqd_angle.avg()) / 1.0_deg},
             {bracket("N"), N.avg()},
             {"bias rejection rate", double(bias_rejected) / cnt},
             {"clusterdistribution", clusterSizeDistribution}};
        _roundjson(j, 3);

        // print threshold matrix
        auto &_j = j["threshold"];
        for (auto i : ids)
            for (auto j : ids)
                if (i >= j) {
                    auto str = Faunus::molecules[i].name + " " + Faunus::molecules[j].name;
                    _j[str] = std::sqrt(thresholdsq(i, j));
                    _roundjson(_j[str], 3);
                }

        // print satellite molecules
        if (not satellites.empty()) {
            auto &_j = j["satellites"];
            _j = json::array();
            for (auto id : satellites)
                _j.push_back(Faunus::molecules[id].name);
        }
    }

    void _from_json(const json &j) override {
        assertKeys(j, {"dp", "dprot", "dir", "threshold", "molecules", "repeat", "satellites"});
        dptrans = j.at("dp");
        dir = j.value("dir", Point(1, 1, 1));
        dprot = j.at("dprot");
        names = j.at("molecules").get<decltype(names)>(); // molecule names
        ids = names2ids(molecules, names);                // names --> molids
        index.clear();
        for (auto &g : spc.groups)            // loop over all groups
            if (not g.atomic)                 // only molecular groups
                if (g.size() == g.capacity()) // only active particles
                    if (std::find(ids.begin(), ids.end(), g.id) != ids.end())
                        index.push_back(&g - &spc.groups.front());
        if (repeat < 0)
            repeat = index.size();

        // read satellite ids (molecules NOT to be considered as cluster centers)
        auto satnames = j.value("satellites", std::vector<std::string>()); // molecule names
        auto vec = names2ids(Faunus::molecules, satnames);                 // names --> molids
        satellites = std::set<int>(vec.begin(), vec.end());

        for (auto id : satellites)
            if (std::find(ids.begin(), ids.end(), id) == ids.end())
                throw std::runtime_error("satellite molecules must be defined in `molecules`");

        // read cluster thresholds
        if (j.count("threshold") == 1) {
            auto &_j = j.at("threshold");
            // threshold is given as a single number
            if (_j.is_number()) {
                for (auto i : ids)
                    for (auto j : ids)
                        if (i >= j)
                            thresholdsq.set(i, j, std::pow(_j.get<double>(), 2));
            }
            // threshold is given as pairs of clustering molecules
            else if (_j.is_object()) {
                for (auto it = _j.begin(); it != _j.end(); ++it) {
                    auto v = words2vec<std::string>(it.key());
                    if (v.size() == 2) {
                        auto it1 = findName(Faunus::molecules, v[0]);
                        auto it2 = findName(Faunus::molecules, v[1]);
                        if (it1 == Faunus::molecules.end() or it2 == Faunus::molecules.end())
                            throw std::runtime_error("unknown molecule(s): ["s + v[0] + " " + v[1] + "]");
                        thresholdsq.set(it1->id(), it2->id(), std::pow(it.value().get<double>(), 2));
                    } else
                        throw std::runtime_error("threshold requires exactly two space-separated molecules");
                }
            } else
                throw std::runtime_error("threshold must be a number or object");
        }
    }

    /**
     * @param spc Space
     * @param first Index of initial molecule (randomly selected)
     * @param index w. all molecules clustered around first (first included)
     */
    void findCluster(Tspace &spc, size_t first, std::set<size_t> &cluster) {
        assert(first < spc.p.size());
        std::set<size_t> pool(index.begin(), index.end());
        assert(pool.count(first) > 0);

        cluster.clear();
        cluster.insert(first);
        pool.erase(first);

        size_t n;
        do { // find cluster (not very clever...)
        start:
            n = cluster.size();
            for (size_t i : cluster)
                if (not spc.groups.at(i).empty()) // check if group is inactive
                    for (size_t j : pool)
                        if (i != j)
                            if (not spc.groups.at(j).empty()) { // check if group is inactive
                                // probability to cluster
                                double P = clusterProbability(spc.groups.at(i), spc.groups.at(j));
                                if (Movebase::slump() <= P) {
                                    cluster.insert(j);
                                    pool.erase(j);
                                    goto start; // wow, first goto ever!
                                }
                            }
        } while (cluster.size() != n);

        // check if cluster is too large
        double max = spc.geo.getLength().minCoeff() / 2;
        for (auto i : cluster)
            for (auto j : cluster)
                if (j > i)
                    if (spc.geo.sqdist(spc.groups.at(i).cm, spc.groups.at(j).cm) >= max * max)
                        rotate = false; // skip rotation if cluster larger than half the box length
    }

    void _move(Change &change) override {
        _bias = 0;
        rotate = true;
        if (not index.empty()) {
            std::set<size_t> cluster; // all group index in cluster

            // find "nuclei" or cluster center and exclude any molecule id listed as "satellite".
            size_t first;
            do {
                first = *slump.sample(index.begin(), index.end()); // random molecule (nuclei)
            } while (satellites.count(spc.groups[first].id) != 0);

            findCluster(spc, first, cluster); // find cluster around first

            N += cluster.size();                       // average cluster size
            clusterSizeDistribution[cluster.size()]++; // update cluster size distribution
            Change::data d;
            d.all = true;
            dp = ranunit(slump, dir) * dptrans * slump();

            if (rotate)
                angle = dprot * (slump() - 0.5);
            else
                angle = 0;

            auto boundary = spc.geo.getBoundaryFunc();

            // lambda function to calculate cluster COM
            auto clusterCOM = [&]() {
                double m, sum_m = 0;
                Point cm(0, 0, 0);
                Point O = spc.groups[*cluster.begin()].cm;
                for (auto i : cluster) {
                    auto &g = spc.groups[i];
                    Point t = g.cm - O;
                    boundary(t);
                    m = g.mass();
                    cm += m * t;
                    sum_m += m;
                }
                cm = cm / sum_m + O;
                boundary(cm);
                return cm;
            };

            Point COM = clusterCOM(); // org. cluster center
            Eigen::Quaterniond Q;
            Q = Eigen::AngleAxisd(angle, ranunit(slump)); // quaternion

            for (auto i : cluster) { // loop over molecules in cluster
                auto &g = spc.groups[i];
                if (rotate) {
                    Geometry::rotate(g.begin(), g.end(), Q, boundary, -COM);
                    g.cm = g.cm - COM;
                    boundary(g.cm);
                    g.cm = Q * g.cm + COM;
                    boundary(g.cm);
                }
                g.translate(dp, boundary);
                d.index = i;
                change.groups.push_back(d);
            }

            change.moved2moved = false; // do not calc. internal cluster energy

            // Reject if cluster composition changes during move
            // Note: this only works for the binary 0/1 probability function
            // currently implemented in `findCluster()`.

            std::set<size_t> aftercluster;         // all group index in cluster _after_move
            findCluster(spc, first, aftercluster); // find cluster around first
            if (aftercluster == cluster)
                _bias = 0;
            else {
                _bias = pc::infty; // bias is infinite --> reject
                bias_rejected++;   // count how many time we reject due to bias
            }
#ifndef NDEBUG
            // check if cluster mass center movement matches displacement
            if (_bias == 0) {
                Point newCOM = clusterCOM();          // org. cluster center
                Point d = spc.geo.vdist(COM, newCOM); // distance between new and old COM
                double _zero = (d + dp).norm();       // |d+dp| should ideally be zero...
                assert(std::fabs(_zero) < 1e-6 && "cluster likely too large");
            }
#endif
        }
    }

    double bias(Change &, double, double) override {
        return _bias;
    } //!< adds extra energy change not captured by the Hamiltonian

    void _reject(Change &) override {
        msqd += 0;
        msqd_angle += 0;
    }

    void _accept(Change &) override {
        msqd += dp.squaredNorm();
        msqd_angle += angle * angle;
    }

  public:
    Cluster(Tspace &spc) : spc(spc) {
        cite = "doi:10/cj9gnn";
        name = "cluster";
        repeat = -1; // meaning repeat N times
    }
};

/**
 * @brief An abstract base class for rotational movements of a polymer chain
 */
class ChainRotationMovebase : public Movebase {
  protected:
    std::string molname;
    size_t molid;
    double dprot;             //!< maximal angle of rotation, ±0.5*dprot
    double sqdispl;           //!< center-of-mass displacement squared
    Average<double> msqdispl; //!< center-of-mass mean squared displacement
    bool permit_move = true;
    bool allow_small_box = false;
    int small_box_encountered = 0; //!< number of skipped moves due to too small container
  private:
    virtual size_t select_segment() = 0;           //!< selects a chain segment and return the number of atoms in it
    virtual void rotate_segment(double angle) = 0; //!< rotates the selected chain segment
    virtual void store_change(Change &change) = 0; //!< stores changes made to atoms
    void _move(Change &change) override;
    void _accept(Change &change) override;
    void _reject(Change &change) override;

  protected:
    void _from_json(const json &j) override;
    void _to_json(json &j) const override;

  public:
    double bias(Change &, double uold, double unew) override;
};

/**
 * @brief An abstract class that rotates a selected segment of a polymer chain in the given simulation box.
 */
template <typename Tspace> class ChainRotationMove : public ChainRotationMovebase {
    using Tbase = ChainRotationMovebase;

  protected:
    Tspace &spc;
    typedef typename Tspace::Tpvec Tpvec;
    typename Tspace::Tgvec::iterator molecule_iter;
    //! Indices of atoms in the spc.p vector that mark the origin and the direction of the axis of rotation.
    std::array<size_t, 2> axis_ndx;
    //! Indices of atoms in the spc.p vector that shall be rotated.
    std::vector<size_t> segment_ndx;

  public:
    explicit ChainRotationMove(Tspace &spc) : spc(spc) { repeat = -1; }

  protected:
    void _from_json(const json &j) override {
        Tbase::_from_json(j);
        auto moliter = findName(molecules, molname);
        if (moliter == molecules.end())
            throw std::runtime_error("unknown molecule '" + molname + "'");
        molid = moliter->id();
    }

  private:
    /**
     * @brief Rotates the chain segment around the axes by the given angle.
     * @param angle
     */
    void rotate_segment(double angle) override {
        if (!segment_ndx.empty()) {
            auto &chain = *molecule_iter;
            auto old_cm = chain.cm;
            // Uses an implementation from the old Pivot class. The translation of the chain might be unnecessary.
            auto shift_pos = spc.p[axis_ndx[0]].pos;
            // chain.unwrap(spc.geo.getDistanceFunc()); // remove pbc
            chain.translate(-shift_pos, spc.geo.getBoundaryFunc());
            auto origin_pos = spc.p[axis_ndx[0]].pos; // != shift_pos because of chain.translate
            auto axis_pos = spc.geo.vdist(origin_pos, spc.p[axis_ndx[1]].pos).normalized();
            Eigen::Quaterniond Q(Eigen::AngleAxisd(angle, axis_pos));
            auto M = Q.toRotationMatrix();
            for (auto i : segment_ndx) {
                spc.p[i].rotate(Q, M);                                       // internal rot.
                spc.p[i].pos = Q * (spc.p[i].pos - origin_pos) + origin_pos; // positional rot.
            }
            chain.cm = Geometry::massCenter(chain.begin(), chain.end());
            chain.translate(shift_pos, spc.geo.getBoundaryFunc());
            // chain.wrap(spc.geo.getBoundaryFunc()); // re-apply pbc
            if (box_big_enough()) {
                sqdispl = spc.geo.sqdist(chain.cm, old_cm); // CM movement
            }
        }
    }

    /**
     * Stores changes of atoms after the move attempt.
     * @param change
     */
    void store_change(Change &change) override {
        if (!segment_ndx.empty()) {
            auto &chain = *molecule_iter;
            auto offset = std::distance(spc.p.begin(), chain.begin());
            Change::data change_data;
            for (int i : segment_ndx) {
                change_data.atoms.push_back(i - offset); // `atoms` index are relative to chain
            }
            change_data.index = Faunus::distance(spc.groups.begin(), &chain); // integer *index* of moved group
            change_data.all = false;
            change_data.internal = true;          // trigger internal interactions
            change.groups.push_back(change_data); // add to list of moved groups
        }
    }

    /** In periodic systems (cuboid, slit, etc.) a chain rotational move can cause the molecule to be larger
     *  than half the box length which we catch here.
     *  @throws std::runtime_error
     */
    bool box_big_enough() {
        auto &chain = *molecule_iter;
        auto cm_pbc = Geometry::massCenter(chain.begin(), chain.end(), spc.geo.getBoundaryFunc(), -chain.cm);
        double cm_diff = spc.geo.sqdist(chain.cm, cm_pbc);
        if (cm_diff > 1e-6) {
            small_box_encountered++;
            permit_move = false;
            if (!allow_small_box) {
                throw std::runtime_error("Container too small for the molecule '" + molname + "'");
            }
            return false;
        }
        return true;
    }
};

/**
 * @brief Performs a crankshaft move of a random segment in a polymer chain.
 *
 * Two random atoms are select from a random polymer chain to be the joints of a crankshaft.
 * The chain segment between the joints is then rotated by a random angle around the axis
 * determined by the joints. The extend of the angle is limited by dprot.
 */
template <typename Tspace> class CrankshaftMove : public ChainRotationMove<Tspace> {
    using Tbase = ChainRotationMove<Tspace>;

  public:
    explicit CrankshaftMove(Tspace &spc) : ChainRotationMove<Tspace>(spc) { this->name = "crankshaft"; }

  protected:
    void _from_json(const json &j) override {
        Tbase::_from_json(j);
        if (this->repeat < 0) {
            // set the number of repetitions to the length of the chain (minus 2) times the number of the chains
            auto moliter = this->spc.findMolecules(this->molid);
            auto &molecule = *moliter.begin();
            this->repeat = std::distance(moliter.begin(), moliter.end()) * (molecule.size() - 2);
        }
    }

  private:
    /** Randomly selects two atoms as joints in a random chain. The joints then determine the axis of rotation
     *  of the chain segment between the joints.
     *  The member vectors containing atoms' indices of the axis and the segment are populated accordingly.
     *  Returns the segment size as atom count.
     *  A non-branched chain is assumed having atom indices in a dense sequence.
     */
    size_t select_segment() override {
        size_t segment_size = 0;
        this->segment_ndx.clear();
        this->molecule_iter = this->spc.randomMolecule(this->molid, this->slump); // a random chain
        if (this->molecule_iter != this->spc.groups.end()) {
            auto &molecule = *this->molecule_iter;
            if (molecule.size() > 2) { // must have at least three atoms
                auto joint0 = this->slump.sample(molecule.begin(), molecule.end());
                auto joint1 = this->slump.sample(molecule.begin(), molecule.end());
                if (joint0 != molecule.end() && joint1 != molecule.end()) {
                    auto joint_distance = std::distance(joint0, joint1);
                    if (joint_distance < 0) {
                        joint_distance *= -1;
                        std::swap(joint0, joint1);
                    }
                    if (joint_distance > 1) { // at least one atom between the joints
                        auto joint0_ndx = std::distance(this->spc.p.begin(), joint0);
                        auto joint1_ndx = std::distance(this->spc.p.begin(), joint1);
                        if (joint0_ndx < 0 || joint1_ndx < 0) {
                            throw std::range_error("A negative index of the atom encountered.");
                        }
                        this->axis_ndx = {(size_t)joint0_ndx, (size_t)joint1_ndx}; // joints create the axis
                        for (size_t i = joint0_ndx + 1; i < joint1_ndx; i++)
                            this->segment_ndx.push_back(i); // add segment's atom indices
                        segment_size = this->segment_ndx.size();
                    }
                }
            }
        }
        return segment_size;
    }
};

/**
 * @brief Performs a pivot move of a random tail part of a polymer chain.
 *
 * A random harmonic bond is selected from a random polymer chain. The bond determines the axes the rotation.
 * A part of the chain either before or after the bond (the selection has an equal probability )
 * then constitutes a segment which is rotated by a random angle. The extend of the angle is limited by dprot.
 */
template <typename Tspace> class PivotMove : public ChainRotationMove<Tspace> {
    using Tbase = ChainRotationMove<Tspace>;

  private:
    std::vector<std::shared_ptr<Potential::BondData>> bonds;

  public:
    explicit PivotMove(Tspace &spc) : ChainRotationMove<Tspace>(spc) { this->name = "pivot"; }

  protected:
    void _from_json(const json &j) override {
        Tbase::_from_json(j);
        bonds = Potential::filterBonds(molecules[this->molid].bonds, Potential::BondData::HARMONIC);

        if (this->repeat < 0) {
            // set the number of repetitions to the length of the chain (minus 2) times the number of the chains
            auto moliter = this->spc.findMolecules(this->molid);
            this->repeat = std::distance(moliter.begin(), moliter.end());
            if (this->repeat > 0) {
                this->repeat *= bonds.size();
            }
        }
    }

    /** Selects a random harmonic bond of a random polymer chain which atoms then create an axis of rotation.
     *  Atoms between the randomly selected chain's end and the bond atom compose a segment to be rotated.
     *  The member vectors containing atoms' indices of the axis and the segment are populated accordingly.
     *  Returns the segment size as atom count.
     *  A non-branched chain is assumed having atom indices in a dense sequence.
     */
    size_t select_segment() override {
        size_t segment_size = 0;
        this->segment_ndx.clear();
        this->molecule_iter = this->spc.randomMolecule(this->molid, this->slump);
        if (this->molecule_iter != this->spc.groups.end()) {
            auto &chain = *this->molecule_iter;
            if (chain.size() > 2) {                                         // must have at least three atoms
                auto bond = this->slump.sample(bonds.begin(), bonds.end()); // a random harmonic bond
                if (bond != bonds.end()) {
                    auto chain_offset = std::distance(this->spc.p.begin(), chain.begin());
                    auto atom0_ndx = (*bond)->index.at(0) + chain_offset;
                    auto atom1_ndx = (*bond)->index.at(1) + chain_offset;
                    if (atom0_ndx < 0 || atom1_ndx < 0) {
                        throw std::range_error("A negative index of the atom occured.");
                    }
                    if (this->slump() > 0.5) {
                        for (size_t i = atom1_ndx + 1; i < chain_offset + chain.size(); i++)
                            this->segment_ndx.push_back(i);
                        std::swap(atom0_ndx, atom1_ndx); // the second atom is the origin
                    } else {
                        for (size_t i = chain_offset; i < atom0_ndx; i++)
                            this->segment_ndx.push_back(i);
                    }
                    this->axis_ndx = std::array<size_t, 2>{(size_t)atom0_ndx, (size_t)atom1_ndx};
                    segment_size = this->segment_ndx.size();
                }
            }
        }
        return segment_size;
    }
};

#ifdef ENABLE_MPI
/**
 * @brief Class for parallel tempering (aka replica exchange) using MPI
 *
 * Although not completely correct, the recommended way of performing a temper move
 * is to do `N` Monte Carlo passes with regular moves and then do a tempering move.
 * This is because the MPI nodes must be in sync and if you have a system where
 * the random number generator calls are influenced by the Hamiltonian we could
 * end up in a deadlock.
 *
 * @date Lund 2012, 2018
 */
class ParallelTempering : public Movebase {
  private:
    typedef typename Tspace::Tpvec Tpvec;

    Tspace &spc; // Space to operate on
    MPI::MPIController &mpi;

    int partner;                   //!< Exchange replica (partner)
    enum extradata { VOLUME = 0 }; //!< Structure of extra data to send
    std::map<std::string, Average<double>> accmap;

    MPI::FloatTransmitter ft;           //!< Class for transmitting floats over MPI
    MPI::ParticleTransmitter<Tpvec> pt; //!< Class for transmitting particles over MPI

    void findPartner(); //!< Find replica to exchange with
    bool goodPartner(); //!< Is partner valid?
    void _to_json(json &j) const override;
    void _move(Change &change) override;
    double exchangeEnergy(double mydu); //!< Exchange energy with partner
    double bias(Change &, double uold, double unew) override;
    std::string id(); //!< Unique string to identify set of partners
    void _accept(Change &) override;
    void _reject(Change &) override;
    void _from_json(const json &j) override;

  public:
    ParallelTempering(Tspace &spc, MPI::MPIController &mpi);
};
#endif

class Propagator : public BasePointerVector<Movebase> {
  private:
    int _repeat;
    std::discrete_distribution<> dist;
    std::vector<double> w; // list of weights for each move

    void addWeight(double weight = 1);

  public:
    using BasePointerVector<Movebase>::vec;
    Propagator() = default;

    Propagator(const json &j, Tspace &spc, MPI::MPIController &mpi);

    int repeat() { return _repeat; }

    auto sample() {
        if (!vec.empty()) {
            assert(w.size() == vec.size());
            return vec.begin() + dist(Move::Movebase::slump.engine);
        }
        return vec.end();
    } //!< Pick move from a weighted, random distribution
};

} // namespace Move

class MCSimulation {
  private:
    typedef typename Tspace::Tpvec Tpvec;

    std::string lastMoveName; //!< name of latest move

    bool metropolis(double du) const; //!< Metropolis criterion (true=accept)

    struct State {
        Tspace spc;
        Energy::Hamiltonian pot;
        State(const json &j);

        void sync(State &other, Change &change);
    }; //!< Contains everything to describe a state

    State state1, // old state (accepted)
        state2;   // new state (trial)
    double uinit = 0, dusum = 0;
    Average<double> uavg;

    void init();

  public:
    Move::Propagator moves;

    auto &pot() { return state1.pot; }
    auto &space() { return state1.spc; }
    const auto &pot() const { return state1.pot; }
    const auto &space() const { return state1.spc; }
    const auto &geometry() const { return state1.spc.geo; }
    const auto &particles() const { return state1.spc.p; }

    MCSimulation(const json &j, MPI::MPIController &mpi);
    double drift(); //!< Calculates the relative energy drift from initial configuration

    /* currently unused -- see Analysis::SaveState.
                    void store(json &j) const {
                        j = state1.spc;
                        j["random-move"] = Move::Movebase::slump;
                        j["random-global"] = Faunus::random;
                    } // store system to json object
    */
    void restore(const json &j); //!< restore system from previously store json object
    void move();
    void to_json(json &j);
};

void to_json(json &j, MCSimulation &mc);

/**
 * @brief Ideal energy contribution of a speciation move
 * This funciton calculates the contribution to the energy change arising from the
 * change in concentration of reactant and products in the current and in the trial state.
 *
 * @f[
 *     \beta \Delta U = - \sum \ln ( N_o!/N_n! V^{N_n - N_o} )
 * @f]
 *
 * where the sum runs over all products and reactants.
 *
 * @todo
 * - use exception message to suggest how to fix the problem
 */
double IdealTerm(Tspace &spc_n, Tspace &spc_o, const Change &change);

} // namespace Faunus
