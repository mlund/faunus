#include "chainmove.h"
#include "iteratorsupport.h"

namespace Faunus {
namespace Move {

void ChainRotationMovebase::_from_json(const json &j) {
    molname = j.at("molecule");
    dprot = j.at("dprot");
    allow_small_box = j.value("skiplarge", true); // todo rename the json attribute and make false default
}

void ChainRotationMovebase::_to_json(json &j) const {
    using namespace u8;
    j = {{"molecule", molname},
         {"dprot", dprot},
         {u8::rootof + u8::bracket("r_cm" + u8::squared), std::sqrt(msqdispl.avg())}};
    if (small_box_encountered > 0) {
        j["skipped"] = double(small_box_encountered) / cnt; // todo rename the json attribute
    }
    _roundjson(j, 3);
}

void ChainRotationMovebase::_move(Change &change) {
    permit_move = true;
    sqdispl = 0;
    if (std::fabs(dprot) > 1e-9) {
        if (select_segment() > 0) {
            double angle = dprot * (slump() - 0.5);
            rotate_segment(angle);
            store_change(change);
        }
    }
}

void ChainRotationMovebase::_accept(Change &) { msqdispl += sqdispl; }
void ChainRotationMovebase::_reject(Change &) { msqdispl += 0; }
double ChainRotationMovebase::bias(Change &, double, double) { return permit_move ? 0 : pc::infty; }

ChainRotationMove::ChainRotationMove(Space &spc) : spc(spc) { repeat = -1; }
void ChainRotationMove::_from_json(const json &j) {
    Tbase::_from_json(j);
    auto moliter = findName(molecules, molname);
    if (moliter == molecules.end())
        throw std::runtime_error("unknown molecule '" + molname + "'");
    molid = moliter->id();
}
void ChainRotationMove::rotate_segment(double angle) {
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
void ChainRotationMove::store_change(Change &change) {
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
bool ChainRotationMove::box_big_enough() {
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
void CrankshaftMove::_from_json(const json &j) {
    Tbase::_from_json(j);
    if (this->repeat < 0) {
        // set the number of repetitions to the length of the chain (minus 2) times the number of the chains
        auto moliter = this->spc.findMolecules(this->molid);
        auto &molecule = *moliter.begin();
        this->repeat = std::distance(moliter.begin(), moliter.end()) * (molecule.size() - 2);
    }
}
size_t CrankshaftMove::select_segment() {
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
PivotMove::PivotMove(Space &spc) : ChainRotationMove(spc) { this->name = "pivot"; }
void PivotMove::_from_json(const json &j) {
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
size_t PivotMove::select_segment() {
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

} // end of namespace Move
} // namespace Faunus
