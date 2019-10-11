#pragma once

#include "move.h"

namespace Faunus {
namespace Move {
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
class ChainRotationMove : public ChainRotationMovebase {
    using Tbase = ChainRotationMovebase;

  protected:
    Space &spc;
    typename Space::Tgvec::iterator molecule_iter;
    //! Indices of atoms in the spc.p vector that mark the origin and the direction of the axis of rotation.
    std::array<size_t, 2> axis_ndx;
    //! Indices of atoms in the spc.p vector that shall be rotated.
    std::vector<size_t> segment_ndx;

  public:
    explicit ChainRotationMove(Space &spc);

  protected:
    void _from_json(const json &j) override;

  private:
    /**
     * @brief Rotates the chain segment around the axes by the given angle.
     * @param angle
     */
    void rotate_segment(double angle) override;

    /**
     * Stores changes of atoms after the move attempt.
     * @param change
     */
    void store_change(Change &change) override;

    /** In periodic systems (cuboid, slit, etc.) a chain rotational move can cause the molecule to be larger
     *  than half the box length which we catch here.
     *  @throws std::runtime_error
     */
    bool box_big_enough();
};

/**
 * @brief Performs a crankshaft move of a random segment in a polymer chain.
 *
 * Two random atoms are select from a random polymer chain to be the joints of a crankshaft.
 * The chain segment between the joints is then rotated by a random angle around the axis
 * determined by the joints. The extend of the angle is limited by dprot.
 */
class CrankshaftMove : public ChainRotationMove {
    using Tbase = ChainRotationMove;

  public:
    explicit CrankshaftMove(Space &spc) : ChainRotationMove(spc) { this->name = "crankshaft"; }

  protected:
    void _from_json(const json &j) override;

  private:
    size_t joint_max; //!< maximum number of bonds between the joints of a crankshaft
    /** Randomly selects two atoms as joints in a random chain. The joints then determine the axis of rotation
     *  of the chain segment between the joints.
     *  The member vectors containing atoms' indices of the axis and the segment are populated accordingly.
     *  Returns the segment size as atom count.
     *  A non-branched chain is assumed having atom indices in a dense sequence.
     */
    size_t select_segment() override;
};

/**
 * @brief Performs a pivot move of a random tail part of a polymer chain.
 *
 * A random harmonic bond is selected from a random polymer chain. The bond determines the axes the rotation.
 * A part of the chain either before or after the bond (the selection has an equal probability )
 * then constitutes a segment which is rotated by a random angle. The extend of the angle is limited by dprot.
 */
class PivotMove : public ChainRotationMove {
    using Tbase = ChainRotationMove;

  private:
    std::vector<std::shared_ptr<Potential::BondData>> bonds;

  public:
    explicit PivotMove(Space &spc);

  protected:
    void _from_json(const json &j) override;

    /** Selects a random harmonic bond of a random polymer chain which atoms then create an axis of rotation.
     *  Atoms between the randomly selected chain's end and the bond atom compose a segment to be rotated.
     *  The member vectors containing atoms' indices of the axis and the segment are populated accordingly.
     *  Returns the segment size as atom count.
     *  A non-branched chain is assumed having atom indices in a dense sequence.
     */
    size_t select_segment() override;
};

} // namespace Move
} // namespace Faunus
