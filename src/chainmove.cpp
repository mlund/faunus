#include "chainmove.h"
#include "aux/iteratorsupport.h"
#include "bonds.h"

namespace Faunus {
namespace Move {

void ChainRotationMoveBase::_from_json(const json &j) {
    molname = j.at("molecule");
    dprot = j.at("dprot");
    allow_small_box = j.value("skiplarge", true); // todo rename the json attribute and make false default
}

void ChainRotationMoveBase::_to_json(json &j) const {
    using namespace u8;
    j = {{"molecule", molname},
         {"dprot", dprot},
         {u8::rootof + u8::bracket("r_cm" + u8::squared), std::sqrt(msqdispl.avg())}};
    if (small_box_encountered > 0) {
        j["skipped"] = double(small_box_encountered) / number_of_attempted_moves; // todo rename the json attribute
    }
    roundJSON(j, 3);
}

void ChainRotationMoveBase::_move(Change &change) {
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

void ChainRotationMoveBase::_accept(Change &) { msqdispl += sqdispl; }
void ChainRotationMoveBase::_reject(Change &) { msqdispl += 0; }
double ChainRotationMoveBase::bias(Change &, double, double) { return permit_move ? 0 : pc::infty; }

ChainRotationMoveBase::ChainRotationMoveBase(Space &spc, std::string name, std::string cite)
    : MoveBase(spc, name, cite) {}

ChainRotationMove::ChainRotationMove(Space &spc, std::string name, std::string cite)
    : ChainRotationMoveBase(spc, name, cite) {
    repeat = -1;
}

void ChainRotationMove::_from_json(const json &j) {
    TBase::_from_json(j);
    molid = findMoleculeByName(molname).id();
}

void ChainRotationMove::rotate_segment(double angle) {
    if (!segment_ndx.empty()) {
        auto &chain = *molecule_iter;
        auto old_cm = chain.mass_center;
        // Uses an implementation from the old Pivot class. The translation of the chain might be unnecessary.
        auto shift_pos = spc.particles.at(axis_ndx[0]).pos;
        // chain.unwrap(spc.geo.getDistanceFunc()); // remove pbc
        chain.translate(-shift_pos, spc.geometry.getBoundaryFunc());
        auto origin_pos = spc.particles.at(axis_ndx[0]).pos; // != shift_pos because of chain.translate
        auto axis_pos = spc.geometry.vdist(origin_pos, spc.particles.at(axis_ndx[1]).pos).normalized();
        Eigen::Quaterniond Q(Eigen::AngleAxisd(angle, axis_pos));
        auto M = Q.toRotationMatrix();
        for (auto index : segment_ndx) {
            auto& particle = spc.particles.at(index);
            particle.rotate(Q, M);                                       // internal rot.
            particle.pos = Q * (particle.pos - origin_pos) + origin_pos; // positional rot.
        }
        chain.mass_center = Geometry::massCenter(chain.begin(), chain.end());
        chain.translate(shift_pos, spc.geometry.getBoundaryFunc());
        // chain.wrap(spc.geo.getBoundaryFunc()); // re-apply pbc
        if (box_big_enough()) {
            sqdispl = spc.geometry.sqdist(chain.mass_center, old_cm); // CM movement
        }
    }
}
void ChainRotationMove::store_change(Change &change) {
    if (!segment_ndx.empty()) {
        auto &chain = *molecule_iter;
        auto offset = std::distance(spc.particles.begin(), chain.begin());
        Change::GroupChange change_data;
        for (int i : segment_ndx) {
            change_data.relative_atom_indices.push_back(i - offset); // `atoms` index are relative to chain
        }
        change_data.group_index = Faunus::distance(spc.groups.begin(), &chain); // integer *index* of moved group
        change_data.all = false;
        change_data.internal = true;          // trigger internal interactions
        change.groups.push_back(change_data); // add to list of moved groups
    }
}
bool ChainRotationMove::box_big_enough() {
    auto &chain = *molecule_iter;
    auto cm_pbc = Geometry::massCenter(chain.begin(), chain.end(), spc.geometry.getBoundaryFunc(), -chain.mass_center);
    double cm_diff = spc.geometry.sqdist(chain.mass_center, cm_pbc);
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

CrankshaftMove::CrankshaftMove(Space &spc, std::string name, std::string cite) : ChainRotationMove(spc, name, cite) {}

CrankshaftMove::CrankshaftMove(Space &spc) : CrankshaftMove(spc, "crankshaft", "") {}

void CrankshaftMove::_from_json(const json &j) {
    TBase::_from_json(j);
    // maximum number of bonds between the joints of a crankshaft
    joint_max = j.value("joint_max", std::numeric_limits<decltype(joint_max)>::max());
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
            assert(joint0 >= molecule.begin() && joint0 < molecule.end());
            size_t range_begin = std::min(static_cast<size_t>(std::distance(molecule.begin(), joint0)), joint_max);
            size_t range_end = std::min(static_cast<size_t>(std::distance(joint0, molecule.end()) - 1), joint_max);
            auto joint1 = joint0 + this->slump.range(-range_begin, range_end);
            assert(joint1 >= molecule.begin() && joint1 < molecule.end());
            auto joint_distance = std::distance(joint0, joint1);
            if (joint_distance < 0) {
                joint_distance *= -1;
                std::swap(joint0, joint1);
            }
            if (joint_distance > 1) { // at least one atom between the joints
                auto joint0_ndx = std::distance(this->spc.particles.begin(), joint0);
                auto joint1_ndx = std::distance(this->spc.particles.begin(), joint1);
                if (joint0_ndx < 0 || joint1_ndx < 0) {
                    throw std::range_error("A negative index of the atom encountered.");
                }
                // joints create the axis
                this->axis_ndx = {static_cast<size_t>(joint0_ndx), static_cast<size_t>(joint1_ndx)};
                for (size_t i = joint0_ndx + 1; i < joint1_ndx; i++)
                    this->segment_ndx.push_back(i); // add segment's atom indices
                segment_size = this->segment_ndx.size();
            }
        }
    }
    return segment_size;
}

PivotMove::PivotMove(Space &spc, std::string name, std::string cite) : ChainRotationMove(spc, name, cite) {}

PivotMove::PivotMove(Space &spc) : PivotMove(spc, "pivot", "") {}

void PivotMove::_from_json(const json &j) {
    TBase::_from_json(j);
    bonds = molecules[this->molid].bonds.find<Potential::HarmonicBond>();
    if (bonds.empty()) {
        throw ConfigurationError("no harmonic bonds found for pivot move");
    }
    if (repeat < 0) {
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
                auto chain_offset = std::distance(this->spc.particles.begin(), chain.begin());
                auto atom0_ndx = (*bond)->indices.at(0) + chain_offset;
                auto atom1_ndx = (*bond)->indices.at(1) + chain_offset;
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
