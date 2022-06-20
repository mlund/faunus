#pragma once
#include "energy.h"
#include "nonbonded/accumulate.h"

namespace Faunus {

namespace Potential {
class PairPotentialBase;
}

namespace Energy {
/**
 * @brief Particle pairing to calculate pairẃise interaction using particles' groups internally. Depending on
 * the accumulator provided, raw particle pairs, energy sum, etc. can be obtained.
 *
 * Accumulator is used as the first argument in all methods. Accumulator shall overload '+=' operator to accept a pair
 * of particle references as used in particle2particle method.
 *
 * @remark Method arguments are generally not checked for correctness because of performance reasons.
 *
 * @tparam TCutoff  a cutoff scheme between groups
 * @see InstantEnergyAccumulator, GroupCutoff
 */
template <typename TCutoff> class GroupPairingPolicy {
  protected:
    const Space& spc; //!< a space to operate on
    TCutoff cut;      //!< a cutoff functor that determines if energy between two groups can be ignored

  public:
    /**
     * @param spc
     */
    explicit GroupPairingPolicy(Space& spc)
        : spc(spc)
        , cut(spc.geometry) {}

    void from_json(const json& j) { Energy::from_json(j, cut); }

    void to_json(json& j) const { Energy::to_json(j, cut); }

    /**
     * @brief Add two interacting particles to the accumulator.
     *
     * Due to compiler optimization, the '+=' operator and this function itself may be inlined to significantly
     * improve performance.
     *
     * @tparam TAccumulator  an accumulator with '+=' operator overloaded to add a pair of particles as references
     *                       {T&, T&}
     * @tparam T  an interacting particle
     * @param pair_accumulator  accumulator of interacting pairs of particles
     * @param a  first particle
     * @param b  second particle
     */
    template <RequireEnergyAccumulator TAccumulator, typename T>
    inline void particle2particle(TAccumulator& pair_accumulator, const T& a, const T& b) const {
        pair_accumulator += {std::cref(a), std::cref(b)};
    }

    /**
     * @brief All pairings within a group.
     *
     * All pair interaction within the group are accumulated. The pair exclusions defined in the molecule
     * topology are honoured.
     *
     * @tparam TAccumulator  an accumulator with '+=' operator overloaded to add a pair of particles as references
     *                       {T&, T&}
     * @tparam TGroup
     * @param group
     * @param pair_accumulator  accumulator of interacting pairs of particles
     */
    template <RequireEnergyAccumulator TAccumulator, typename TGroup>
    void groupInternal(TAccumulator& pair_accumulator, const TGroup& group) {
        const auto& moldata = group.traits();
        if (!moldata.rigid) {
            const int group_size = group.size();
            for (int i = 0; i < group_size - 1; ++i) {
                for (int j = i + 1; j < group_size; ++j) {
                    // This compound condition is faster than an outer atomic condition;
                    // tested on bulk example in GCC 9.2.
                    if (group.isAtomic() || !moldata.isPairExcluded(i, j)) {
                        particle2particle(pair_accumulator, group[i], group[j]);
                    }
                }
            }
        }
    }

    /**
     * @brief Pairings of a single particle within the group.
     *
     * The pair exclusions defined in the molecule topology are honoured.
     *
     * @tparam TAccumulator  an accumulator with '+=' operator overloaded to add a pair of particles as references
     *                       {T&, T&}
     * @tparam TGroup
     * @param pair_accumulator  accumulator of interacting pairs of particles
     * @param group
     * @param index  internal index of the selected particle within the group
     */
    template <RequireEnergyAccumulator TAccumulator, typename TGroup>
    void groupInternal(TAccumulator& pair_accumulator, const TGroup& group, const std::size_t index) {
        const auto& moldata = group.traits();
        if (!moldata.rigid) {
            if (group.isAtomic()) {
                // speed optimization: non-bonded interaction exclusions do not need to be checked for atomic groups
                for (int i = 0; i < index; ++i) {
                    particle2particle(pair_accumulator, group[index], group[i]);
                }
                for (int i = index + 1; i < group.size(); ++i) {
                    particle2particle(pair_accumulator, group[index], group[i]);
                }
            } else {
                // molecular group
                for (int i = 0; i < index; ++i) {
                    if (!moldata.isPairExcluded(index, i)) {
                        particle2particle(pair_accumulator, group[index], group[i]);
                    }
                }
                for (int i = index + 1; i < group.size(); ++i) {
                    if (!moldata.isPairExcluded(index, i)) {
                        particle2particle(pair_accumulator, group[index], group[i]);
                    }
                }
            }
        }
    }

    /**
     * @brief Pairing in the group involving only the particles present in the index.
     *
     * Only such non-bonded pair interactions within the group are considered if at least one particle is present
     * in the index. The pair exclusions defined in the molecule topology are honoured.
     *
     * @tparam TAccumulator  an accumulator with '+=' operator overloaded to add a pair of particles as references
     *                       {T&, T&}
     * @tparam TGroup
     * @tparam TIndex
     * @param pair_accumulator  accumulator of interacting pairs of particles
     * @param group
     * @param index  internal indices of particles within the group
     */
    template <RequireEnergyAccumulator TAccumulator, typename TGroup, typename TIndex>
    void groupInternal(TAccumulator& pair_accumulator, const TGroup& group, const TIndex& index) {
        auto& moldata = group.traits();
        if (!moldata.rigid) {
            if (index.size() == 1) {
                groupInternal(pair_accumulator, group, index[0]);
            } else {
                // TODO investigate overhead of `index_complement` filtering;
                // TODO perhaps allow different strategies based on the index-size/group-size ratio
                auto index_complement = indexComplement(group.size(), index);
                // moved <-> static
                for (int i : index) {
                    for (int j : index_complement) {
                        if (!moldata.isPairExcluded(i, j)) {
                            particle2particle(pair_accumulator, group[i], group[j]);
                        }
                    }
                }
                // moved <-> moved
                for (auto i_it = index.begin(); i_it < index.end(); ++i_it) {
                    for (auto j_it = std::next(i_it); j_it < index.end(); ++j_it) {
                        if (!moldata.isPairExcluded(*i_it, *j_it)) {
                            particle2particle(pair_accumulator, group[*i_it], group[*j_it]);
                        }
                    }
                }
            }
        }
    }

    /**
     * @brief Complete cartesian pairing of particles in two groups.
     *
     * group1 × group2
     *
     * If the distance between the groups is greater or equal to the group cutoff distance, no calculation is performed.
     * The group intersection must be an empty set, i.e., no particle is included in both groups. This is not verified
     * for performance reason.
     *
     * @tparam TAccumulator  an accumulator with '+=' operator overloaded to add a pair of particles as references
     *                       {T&, T&}
     * @tparam TGroup
     * @param pair_accumulator  accumulator of interacting pairs of particles
     * @param group1
     * @param group2
     */
    template <RequireEnergyAccumulator TAccumulator, typename TGroup>
    void group2group(TAccumulator& pair_accumulator, const TGroup& group1, const TGroup& group2) {
        if (!cut(group1, group2)) {
            for (auto& particle1 : group1) {
                for (auto& particle2 : group2) {
                    particle2particle(pair_accumulator, particle1, particle2);
                }
            }
        }
    }

    /**
     * @brief Cross pairing of particles in two groups. Only a cartesian subset of the complete cartesian product is
     * considered as the particles in the first group must be also present in the index. The aim is to capture only
     * interactions that involve changing (indexed) particles.
     *
     * ⊕group1 × group2, where ⊕ denotes a filter by an index
     *
     * If the distance between the groups is greater or equal to the group cutoff distance, no calculation is performed.
     * The group intersection must be an empty set, i.e., no particle is included in both groups. This is not verified
     * for performance reason.

     * @tparam TAccumulator  an accumulator with '+=' operator overloaded to add a pair of particles as references
     *                       {T&, T&}
     * @tparam TGroup
     * @param pair_accumulator  accumulator of interacting pairs of particles
     * @param group1
     * @param group2
     * @param index1  list of particle indices in group1 relative to the group beginning
     */
    template <RequireEnergyAccumulator TAccumulator, typename TGroup>
    void group2group(TAccumulator& pair_accumulator, const TGroup& group1, const TGroup& group2,
                     const std::vector<std::size_t>& index1) {
        if (!cut(group1, group2)) {
            for (auto particle1_ndx : index1) {
                for (auto& particle2 : group2) {
                    particle2particle(pair_accumulator, *(group1.begin() + particle1_ndx), particle2);
                }
            }
        }
    }

    /**
     * @brief Cross pairing of particles in two groups. Only a non-cartesian subset of the complete cartesian product
     * is considered as at least one particles in the pair must be also present in the respective index. The aim is
     * to capture only interactions that involve changing (indexed) particles, i.e., to avoid pairs containing only
     * non-indexed particles.
     *
     * (⊕group1 × ∁⊕group2) + (∁⊕group1 × ⊕group2) + (⊕group1 × ⊕group2) =
     * = group1 × group2 − (∁⊕group2 × ∁⊕group2), where ⊕ denotes a filter by an index and ∁ a complement
     *
     * If the distance between the groups is greater or equal to the group cutoff distance, no calculation is performed.
     * The group intersection must be an empty set, i.e., no particle is included in both groups. This is not verified
     * for performance reason.
     *
     * @tparam TAccumulator  an accumulator with '+=' operator overloaded to add a pair of particles as references
     *                       {T&, T&}
     * @tparam TGroup
     * @param pair_accumulator  accumulator of interacting pairs of particles
     * @param group1
     * @param group2
     * @param index1  list of particle indices in group1 relative to the group beginning
     * @param index2  list of particle indices in group2 relative to the group beginning
     */
    template <RequireEnergyAccumulator TAccumulator, typename TGroup>
    void group2group(TAccumulator& pair_accumulator, const TGroup& group1, const TGroup& group2,
                     const std::vector<std::size_t>& index1, const std::vector<std::size_t>& index2) {
        if (!cut(group1, group2)) {
            if (!index2.empty()) {
                // (∁⊕group1 × ⊕group2) + (⊕group1 × ⊕group2) = group1 × ⊕group2
                group2group(pair_accumulator, group2, group1, index2);
                // + (⊕group1 × ∁⊕group2)
                auto index2_complement = indexComplement(group2.size(), index2);
                for (auto particle1_ndx : index1) {
                    for (auto particle2_ndx : index2_complement) {
                        particle2particle(pair_accumulator, group2[particle2_ndx], group1[particle1_ndx]);
                    }
                }
            } else if (!index1.empty()) {
                // (⊕group1 × ∁⊕group2) + (⊕group1 × ⊕group2) = ⊕group1 × group2
                group2group(pair_accumulator, group1, group2, index1);
                // + (∁⊕group1 × ⊕group2) = Ø as ⊕group2 is empty
            } else {
                // both indices empty hence nothing to do
            }
        }
    }

    /**
     * @brief Complete cartesian pairing between particles in a group and a union of groups.
     *
     * group × (∪ groups)
     *
     * If the distance between the groups is greater or equal to the group cutoff distance, the particle pairing between
     * them is skipped. The internal energy of the group is not computed even if the group is also present in the union
     * of groups.
     *
     * @tparam TAccumulator  an accumulator with '+=' operator overloaded to add a pair of particles as references
     *                       {T&, T&}
     * @tparam TGroup
     * @tparam TGroups
     * @param pair_accumulator  accumulator of interacting pairs of particles
     * @param group
     * @param groups
     */
    template <RequireEnergyAccumulator TAccumulator, typename TGroup, typename TGroups>
    void group2groups(TAccumulator& pair_accumulator, const TGroup& group, const TGroups& groups) {
        for (auto& other_group : groups) {
            if (&other_group != &group) {
                group2group(pair_accumulator, group, other_group);
            }
        }
    }

    /**
     * @brief Cross pairing of particles in a group and a union of groups. Only a cartesian subset of the complete
     * cartesian product is considered as the particles of the first group must be also present in the index. The aim
     * is to capture only interactions that involve changing (indexed) particles.
     *
     * ⊕group × (∪ groups), where ⊕ denotes a filter by an index
     *
     * If the distance between the groups is greater or equal to the group cutoff distance, the particle pairing
     * between them is skipped. The internal energy of the group is not computed even if the group is also present
     * in the union of groups.
     *
     * @tparam TAccumulator  an accumulator with '+=' operator overloaded to add a pair of particles as references
     *                       {T&, T&}
     * @tparam TGroup
     * @tparam TGroups
     * @param pair_accumulator  accumulator of interacting pairs of particles
     * @param group
     * @param group_index  groups as indices in Space::groups
     * @param index  list of particle indices in the group relative to the group beginning
     */
    template <RequireEnergyAccumulator TAccumulator, typename TGroup, typename TGroups>
    void group2groups(TAccumulator& pair_accumulator, const TGroup& group, const TGroups& group_index,
                      const std::vector<std::size_t>& index) {
        for (auto other_group_ndx : group_index) {
            const auto& other_group = spc.groups[other_group_ndx];
            if (&other_group != &group) {
                group2group(pair_accumulator, group, other_group, index);
            }
        }
    }

    /**
     * @brief Complete cartesian pairing between particles in a group and particles in other groups in space.
     *
     * group × (space ∖ group)
     *
     * If the distance between the groups is greater or equal to the group cutoff distance, the particle pairing
     * between them is skipped.
     *
     * @tparam TAccumulator  an accumulator with '+=' operator overloaded to add a pair of particles as references
     *                       {T&, T&}
     * @tparam TGroup
     * @param pair_accumulator  accumulator of interacting pairs of particles
     * @param group
     */
    template <RequireEnergyAccumulator TAccumulator, typename Tgroup>
    void group2all(TAccumulator& pair_accumulator, const Tgroup& group) {
        for (auto& other_group : spc.groups) {
            if (&other_group != &group) {
                group2group(pair_accumulator, group, other_group);
            }
        }
    }

    /**
     * @brief Complete cartesian pairing between a single particle in a group and particles in other groups in space.
     *
     * ⊕group × (space ∖ group), where ⊕ denotes a filter by an index (here a single particle)
     *
     * If the distance between the groups is greater or equal to the group cutoff distance, the particle pairing
     * between them is skipped. This method is performance-optimized version of the multiple indices method.
     *
     * @tparam TAccumulator  an accumulator with '+=' operator overloaded to add a pair of particles as references
     *                       {T&, T&}
     * @tparam TGroup
     * @param pair_accumulator  accumulator of interacting pairs of particles
     * @param group
     * @param index  a particle index relative to the group beginning
     */
    template <RequireEnergyAccumulator TAccumulator, typename TGroup>
    void group2all(TAccumulator& pair_accumulator, const TGroup& group, const int index) {
        const auto& particle = group[index];
        for (auto& other_group : spc.groups) {
            if (&other_group != &group) {                      // avoid self-interaction
                if (!cut(other_group, group)) {                // check g2g cut-off
                    for (auto& other_particle : other_group) { // loop over particles in other group
                        particle2particle(pair_accumulator, particle, other_particle);
                    }
                }
            }
        }
    }

    /**
     * @brief Complete cartesian pairing between selected particles in a group and particles in other groups in space.
     *
     * ⊕group × (space ∖ group), where ⊕ denotes a filter by an index
     *
     * If the distance between the groups is greater or equal to the group cutoff distance, the particle pairing
     * between them is skipped.
     *
     * @tparam TAccumulator  an accumulator with '+=' operator overloaded to add a pair of particles as references
     *                       {T&, T&}
     * @param pair_accumulator  accumulator of interacting pairs of particles
     * @param group
     * @param index  list of particle indices in the group relative to the group beginning
     */
    template <RequireEnergyAccumulator TAccumulator, typename Tgroup>
    void group2all(TAccumulator& pair_accumulator, const Tgroup& group, const std::vector<std::size_t>& index) {
        if (index.size() == 1) {
            group2all(pair_accumulator, group, index[0]);
        } else {
            for (auto& other_group : spc.groups) {
                if (&other_group != &group) {
                    group2group(pair_accumulator, group, other_group, index);
                }
            }
        }
    }

    /**
     * @brief Cross pairing of particles among a union of groups. No internal pairs within any group are considered.
     *
     * If the distance between any two groups is greater or equal to the group cutoff distance, the particle pairing
     * between them is skipped.
     *
     * @tparam TAccumulator  an accumulator with '+=' operator overloaded to add a pair of particles as references
     *                       {T&, T&}
     * @tparam T
     * @param pair_accumulator  accumulator of interacting pairs of particles
     * @param group_index  list of groups
     */
    template <RequireEnergyAccumulator TAccumulator, typename T>
    void groups2self(TAccumulator& pair_accumulator, const T& group_index) {
        for (auto group1_ndx_it = group_index.begin(); group1_ndx_it < group_index.end(); ++group1_ndx_it) {
            // no such move exists that the internal energy has to be recalculated
            // groupInternal(pair_accumulator, spc.groups[*group1_ndx_it]);
            for (auto group2_ndx_it = std::next(group1_ndx_it); group2_ndx_it < group_index.end(); group2_ndx_it++) {
                group2group(pair_accumulator, spc.groups[*group1_ndx_it], spc.groups[*group2_ndx_it]);
            }
        }
    }

    /**
     * @brief Cross pairing of particles between a union of groups and its complement in space.
     *
     * If the distance between any two groups is greater or equal to the group cutoff distance, the particle pairing
     * between them is skipped.
     *
     * @tparam TAccumulator  an accumulator with '+=' operator overloaded to add a pair of particles as references
     *                       {T&, T&}
     * @tparam T
     * @param pair_accumulator  accumulator of interacting pairs of particles
     * @param group_index  list of groups
     */
    template <RequireEnergyAccumulator TAccumulator, typename T>
    void groups2all(TAccumulator& pair_accumulator, const T& group_index) {
        groups2self(pair_accumulator, group_index);
        auto index_complement = indexComplement(spc.groups.size(), group_index);
        for (auto group1_ndx : group_index) {
            for (auto group2_ndx : index_complement) {
                group2group(pair_accumulator, spc.groups[group1_ndx], spc.groups[group2_ndx]);
            }
        }
    }

    /**
     * @brief Cross pairing between all particles in the space.
     *
     * If the distance between particles' groups is greater or equal to the group cutoff distance, no calculation is
     * performed.
     *
     * @tparam TAccumulator  an accumulator with '+=' operator overloaded to add a pair of particles as references
     *                       {T&, T&}
     * @param pair_accumulator  accumulator of interacting pairs of particles
     */
    template <RequireEnergyAccumulator TAccumulator> void all(TAccumulator& pair_accumulator) {
        for (auto group_it = spc.groups.begin(); group_it < spc.groups.end(); ++group_it) {
            groupInternal(pair_accumulator, *group_it);
            for (auto other_group_it = std::next(group_it); other_group_it < spc.groups.end(); other_group_it++) {
                group2group(pair_accumulator, *group_it, *other_group_it);
            }
        }
    }

    /**
     * @brief Cross pairing between all particles in the space.
     *
     * If the distance between particles' groups is greater or equal to the group cutoff distance, no calculation is
     * performed.
     *
     * @tparam TAccumulator  an accumulator with '+=' operator overloaded to add a pair of particles as references
     *                       {T&, T&}
     * @tparam TCondition  a function returning bool and having a group as an argument
     * @param pair_accumulator  accumulator of interacting pairs of particles
     * @param condition  a group filter if internal energy of the group shall be added
     */
    template <RequireEnergyAccumulator TAccumulator, typename TCondition>
    void all(TAccumulator& pair_accumulator, TCondition condition) {
        for (auto group_it = spc.groups.begin(); group_it < spc.groups.end(); ++group_it) {
            if (condition(*group_it)) {
                groupInternal(pair_accumulator, *group_it);
            }
            for (auto other_group_it = std::next(group_it); other_group_it < spc.groups.end(); other_group_it++) {
                group2group(pair_accumulator, *group_it, *other_group_it);
            }
        }
    }
};

} // namespace Energy
} // namespace Faunus