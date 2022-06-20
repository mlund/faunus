#pragma once
#include "energy.h"
#include "nonbonded/accumulate.h"
#include "nonbonded/sum_by_group_policy.h"

namespace Faunus {

namespace Potential {
class PairPotentialBase;
}

namespace Energy {

/**
 * @brief Computes pair quantity difference for a systen perturbation. Such quantity can be energy using nonponded
 * pair potential
 * .
 * @tparam TPolicy  a pairing policy
 */
template <typename TPolicy> class GroupPairing {
    const Space& spc;
    TPolicy pairing;

  protected:
    /**
     * @brief Computes pair quantity difference if only a single group has changed.
     *
     * @tparam TAccumulator
     * @param pair_accumulator  accumulator of interacting pairs of particles
     * @param change
     */
    template <RequireEnergyAccumulator TAccumulator>
    void accumulateGroup(TAccumulator& pair_accumulator, const Change& change) {
        const auto& change_data = change.groups.at(0);
        const auto& group = spc.groups.at(change_data.group_index);
        if (change_data.relative_atom_indices.size() == 1) {
            // faster algorithm if only a single particle moves
            pairing.group2all(pair_accumulator, group, change_data.relative_atom_indices[0]);
            if (change_data.internal) {
                pairing.groupInternal(pair_accumulator, group, change_data.relative_atom_indices[0]);
            }
        } else {
            const bool change_all = change_data.relative_atom_indices.empty(); // all particles or only their subset?
            if (change_all) {
                pairing.group2all(pair_accumulator, group);
                if (change_data.internal) {
                    pairing.groupInternal(pair_accumulator, group);
                }
            } else {
                pairing.group2all(pair_accumulator, group, change_data.relative_atom_indices);
                if (change_data.internal) {
                    pairing.groupInternal(pair_accumulator, group, change_data.relative_atom_indices);
                }
            }
        }
    }

    /**
     * @brief Computes pair quantity difference if the number of particles has changed.
     *
     * Particles have to be explicitly enumerated in the atom indices of the changed group. Implicit addition of atoms
     * with a group is not supported yet. Note that we do not have to care about missing (removed) particles at all.
     * They are taken into account in the original (old) space where they are present.
     *
     * @param pair_accumulator  accumulator of interacting pairs of particles
     * @param change
     */
    template <RequireEnergyAccumulator TAccumulator>
    void accumulateSpeciation(TAccumulator& pair_accumulator, const Change& change) {
        assert(change.matter_change);
        const auto& moved = change.touchedGroupIndex(); // index of moved groups
        const auto fixed =
            indexComplement(spc.groups.size(), moved) | ranges::to<std::vector>; // index of static groups
        auto filter_active = [](int size) { return ranges::views::filter([size](const auto i) { return i < size; }); };

        // loop over all changed groups
        for (auto change_group1_it = change.groups.begin(); change_group1_it < change.groups.end();
             ++change_group1_it) {
            const auto& group1 = spc.groups.at(change_group1_it->group_index);
            // filter only active particles
            const auto index1 =
                change_group1_it->relative_atom_indices | filter_active(group1.size()) | ranges::to<std::vector>;
            if (!index1.empty()) {
                // particles added into the group: compute (changed group) <-> (static group)
                pairing.group2groups(pair_accumulator, group1, fixed, index1);
            }
            // loop over successor changed groups (hence avoid double counting group1×group2 and group2×group1)
            for (auto change_group2_it = std::next(change_group1_it); change_group2_it < change.groups.end();
                 ++change_group2_it) {
                const auto& group2 = spc.groups.at(change_group2_it->group_index);
                const auto index2 =
                    change_group2_it->relative_atom_indices | filter_active(group2.size()) | ranges::to<std::vector>;
                if (!index1.empty() || !index2.empty()) {
                    // particles added into one or other group: compute (changed group) <-> (changed group)
                    pairing.group2group(pair_accumulator, group1, group2, index1, index2);
                }
            }
            if (!index1.empty() && !molecules.at(group1.id).rigid) {
                // compute internal energy in the changed group
                if (change_group1_it->all) {
                    pairing.groupInternal(pair_accumulator, group1);
                } else {
                    pairing.groupInternal(pair_accumulator, group1, index1);
                };
            }
        }
    }

  public:
    /**
     * @brief Computes pair quantity difference from changed particles.
     *
     * The internal energy contribution, i.e., the contribution from the intra group interactions, is added
     * only if a single group is changed or if all changed.
     *
     * @param change
     * @param pair_accumulator  accumulator of interacting pairs of particles
     */
    template <RequireEnergyAccumulator TAccumulator>
    void accumulate(TAccumulator& pair_accumulator, const Change& change) {
        assert(std::is_sorted(change.groups.begin(), change.groups.end()));
        if (change.everything) {
            pairing.all(pair_accumulator);
        } else if (change.volume_change) {
            // sum all interaction energies except the internal energies of incompressible molecules
            pairing.all(pair_accumulator, [](auto& group) { return group.isAtomic() || group.traits().compressible; });
        } else if (!change.matter_change) {
            if (change.groups.size() == 1) {
                // if only a single group changes use faster algorithm and optionally add the internal energy
                accumulateGroup(pair_accumulator, change);
            } else {
                // if multiple groups move, no internal energies are computed
                const auto& moved = change.touchedGroupIndex(); // index of moved groups
                pairing.groups2all(pair_accumulator, moved);
            }
        } else { // change.dN
            accumulateSpeciation(pair_accumulator, change);
        }
    }

    explicit GroupPairing(Space& spc)
        : spc(spc)
        , pairing(spc) {}

    void from_json(const json& j) { pairing.from_json(j); }

    void to_json(json& j) const { pairing.to_json(j); }

    // FIXME a temporal fix for non-refactorized NonbondedCached
    template <typename Accumulator>
    void group2group(Accumulator& pair_accumulator, const Space::GroupType& group1, const Space::GroupType& group2) {
        pairing.group2group(std::forward<Accumulator&>(pair_accumulator), std::forward<const Space::GroupType&>(group1),
                            std::forward<const Space::GroupType&>(group2));
    }
};

class NonbondedBase : public Energybase {
  public:
    /**
     * Helper function to retrive particle<->particle energy. Usually used by analysis functions.
     */
    virtual double particleParticleEnergy(const Particle& particle1, const Particle& particle2) = 0;
    /**
     * Helper function to retrive group<->group energy. Usually used by analysis functions.
     */
    virtual double groupGroupEnergy(const Group& group1, const Group& group2) = 0;
};

/**
 * @brief Computes change in the non-bonded energy, assuming pairwise additive energy terms.
 *
 * @tparam TPairEnergy  a functor to compute non-bonded energy between two particles
 * @tparam TPairingPolicy  pairing policy to effectively sum up the pairwise additive non-bonded energy
 */
template <RequirePairEnergy TPairEnergy, typename TPairingPolicy> class Nonbonded : public NonbondedBase {
  protected:
    const Space& spc;        //!< space to operate on
    TPairEnergy pair_energy; //!< a functor to compute non-bonded energy between two particles, see PairEnergy
    TPairingPolicy pairing;  //!< pairing policy to effectively sum up the pairwise additive non-bonded energy
    std::shared_ptr<EnergyAccumulatorBase>
        energy_accumulator; //!< energy accumulator used for storing and summing pair-wise energies

  public:
    Nonbonded(const json& j, Space& spc, BasePointerVector<Energybase>& pot)
        : spc(spc)
        , pair_energy(spc, pot)
        , pairing(spc) {
        name = "nonbonded";
        from_json(j);
        energy_accumulator = createEnergyAccumulator(j, pair_energy, 0.0);
        energy_accumulator->reserve(spc.numParticles()); // attempt to reduce memory fragmentation
    }

    double particleParticleEnergy(const Particle& particle1, const Particle& particle2) override {
        return pair_energy(particle1, particle2);
    }

    void updateState(const Change& change) override;

    void sync(Energybase* other_energy, const Change& change) override;

    double groupGroupEnergy(const Group& group1, const Group& group2) override {
        InstantEnergyAccumulator<TPairEnergy> accumulator(pair_energy);
        pairing.group2group(accumulator, group1, group2);
        return static_cast<double>(accumulator);
    }

    void from_json(const json& j) {
        pair_energy.from_json(j);
        pairing.from_json(j);
    }

    void to_json(json& j) const override {
        pair_energy.to_json(j);
        pairing.to_json(j);
        energy_accumulator->to_json(j);
    }

    double energy(const Change& change) override {
        energy_accumulator->clear();
        // down-cast to avoid slow, virtual function calls:
        if (auto ptr = std::dynamic_pointer_cast<InstantEnergyAccumulator<TPairEnergy>>(energy_accumulator)) {
            pairing.accumulate(*ptr, change);
        } else if (auto ptr = std::dynamic_pointer_cast<DelayedEnergyAccumulator<TPairEnergy>>(energy_accumulator)) {
            pairing.accumulate(*ptr, change);
        } else {
            pairing.accumulate(*energy_accumulator, change);
        }
        return static_cast<double>(*energy_accumulator);
    }

    /**
     * @brief Calculates the force on all particles.
     *
     * @todo A stub. Change to reflect only active particle, see Space::activeParticles().
     */
    void force(std::vector<Point>& forces) override {
        // just a temporary hack; perhaps better to allow PairForce instead of the PairEnergy template
        assert(forces.size() == spc.particles.size() && "the forces size must match the particle size");
        for (size_t i = 0; i < spc.particles.size() - 1; ++i) {
            for (size_t j = i + 1; j < spc.particles.size(); ++j) {
                const Point f = pair_energy.force(spc.particles[i], spc.particles[j]);
                forces[i] += f;
                forces[j] -= f;
            }
        }
    }
};
template <RequirePairEnergy TPairEnergy, typename TPairingPolicy>
void Nonbonded<TPairEnergy, TPairingPolicy>::updateState(const Change& change) {
    energy_accumulator->updateState(change);
}

/**
 * If the energy accumulator depends on the system state (e.g. cell lists) then we need to sync this.
 */
template <RequirePairEnergy TPairEnergy, typename TPairingPolicy>
void Nonbonded<TPairEnergy, TPairingPolicy>::sync(Energybase* other_energy, const Change& change) {
    if (auto ptr = dynamic_cast<decltype(this)>(other_energy)) {
        energy_accumulator->sync(*(ptr->energy_accumulator), change);
    } else {
        throw std::runtime_error("sync error");
    }
}

/**
 * @brief Computes non-bonded energy contribution from changed particles. Cache group2group energy once calculated,
 * until a new trial configuration is provided. Not for general use as only partially implemented!
 *
 * Original implementation, only refurbished. Generally suboptimal as only PairingPolicy::group2group method
 * may be called.
 * No internal energy is ever computed. Cannot deal with particle count changes. And other unmentioned constrains.
 *
 * @tparam TPairEnergy  a functor to compute non-bonded energy between two particles
 * @tparam TPairingPolicy  pairing policy to effectively sum up the pairwise additive non-bonded energy
 */
template <RequirePairEnergy TPairEnergy, typename TPairingPolicy>
class NonbondedCached : public Nonbonded<TPairEnergy, TPairingPolicy> {
    using Base = Nonbonded<TPairEnergy, TPairingPolicy>;
    using TAccumulator = InstantEnergyAccumulator<TPairEnergy>;
    Eigen::MatrixXf energy_cache;
    using Base::spc;

    template <typename TGroup> double g2g(const TGroup& g1, const TGroup& g2) {
        int i = &g1 - spc.groups.data();
        int j = &g2 - spc.groups.data();
        if (j < i) {
            std::swap(i, j);
        }
        if (Energybase::state == Energybase::MonteCarloState::TRIAL) { // if this is from the trial system
            TAccumulator energy_accumulator(Base::pair_energy);
            Base::pairing.group2group(energy_accumulator, g1, g2);
            energy_cache(i, j) = static_cast<double>(energy_accumulator); // update the cache
        }
        return energy_cache(i, j); // return (cached) value
    }

    template <typename TGroup>
    double g2g(const TGroup& g1, const TGroup& g2, [[maybe_unused]] const std::vector<std::size_t>& index) {
        // index not implemented
        return g2g(g1, g2);
    }

  public:
    NonbondedCached(const json& j, Space& spc, BasePointerVector<Energybase>& pot)
        : Base(j, spc, pot) {
        Base::name += "EM";
        init();
    }

    /**
     * @brief Cache pair interactions in matrix.
     */
    void init() override {
        const auto groups_size = spc.groups.size();
        energy_cache.resize(groups_size, groups_size);
        energy_cache.setZero();
        TAccumulator u(Base::pair_energy);
        for (auto i = 0; i < groups_size - 1; ++i) {
            for (auto j = i + 1; j < groups_size; ++j) {
                u = 0.0;
                Base::pairing.group2group(u, spc.groups.at(i), spc.groups.at(j));
                energy_cache(i, j) = static_cast<double>(u);
            }
        }
    }

    double energy(const Change& change) override {
        // Only g2g may be called there to compute (and cache) energy!
        double energy_sum = 0.0;
        if (change) {
            if (change.everything || change.volume_change) {
                for (auto i = spc.groups.begin(); i < spc.groups.end(); ++i) {
                    for (auto j = std::next(i); j < Base::spc.groups.end(); ++j) {
                        energy_sum += g2g(*i, *j);
                    }
                }
            } else {
                if (change.groups.size() == 1) { // if exactly ONE molecule is changed
                    auto& d = change.groups[0];
                    auto& g1 = spc.groups.at(d.group_index);
                    for (auto g2_it = spc.groups.begin(); g2_it < spc.groups.end(); ++g2_it) {
                        if (&g1 != &(*g2_it)) {
                            energy_sum += g2g(g1, *g2_it, d.relative_atom_indices);
                        }
                    }
                } else {                                     // many molecules are changed
                    auto moved = change.touchedGroupIndex(); // index of moved groups
                    // moved<->moved
                    if (change.moved_to_moved_interactions) {
                        for (auto i = moved.begin(); i < moved.end(); ++i) {
                            for (auto j = std::next(i); j < moved.end(); ++j) {
                                energy_sum += g2g(spc.groups[*i], spc.groups[*j]);
                            }
                        }
                    }
                    // moved<->static
#if true
                    // classic version
                    const auto fixed = indexComplement(spc.groups.size(), moved) | ranges::to_vector; // static groups
                    for (auto i : moved) {
                        for (auto j : fixed) {
                            energy_sum += g2g(spc.groups[i], spc.groups[j]);
                        }
                    }
#else
                    // OMP-ready version
                    auto fixed =
                        indexComplement(spc.groups.size(), moved) | ranges::to<std::vector>; // index of static groups
                    const size_t moved_size = moved.size();
                    const size_t fixed_size = fixed.size();
                    for (auto i = 0; i < moved_size; ++i) {
                        for (auto j = 0; j < fixed_size; ++j) {
                            energy_sum += g2g(spc.groups[moved[i]], spc.groups[fixed[j]]);
                        }
                    }
#endif
                }
            }
            // more todo!
        }
        return energy_sum;
    }

    /**
     * @brief Copy energy matrix from other
     * @param base_ptr
     * @param change
     */
    void sync(Energybase* base_ptr, const Change& change) override {
        auto other = dynamic_cast<decltype(this)>(base_ptr);
        assert(other);
        if (change.everything || change.volume_change) {
            energy_cache.triangularView<Eigen::StrictlyUpper>() =
                (other->energy_cache).template triangularView<Eigen::StrictlyUpper>();
        } else {
            for (const auto& d : change.groups) {
                for (int i = 0; i < d.group_index; i++) {
                    energy_cache(i, d.group_index) = other->energy_cache(i, d.group_index);
                }
                for (size_t i = d.group_index + 1; i < spc.groups.size(); i++) {
                    energy_cache(d.group_index, i) = other->energy_cache(d.group_index, i);
                }
            }
        }
    }
};

} // namespace Energy
} // namespace Faunus