#pragma once

/**
 * @brief Cell list class templates.
 *
 * Minimal interface to be used in other header file to avoid number of class templates. Full class
 * templates definitions are provided in "celllistimpl.h" file.
 *
 * @author Richard Chudoba
 * @date 2021-02-01
 */

namespace Faunus {

namespace CellList {

namespace Grid {
/**
 * Type manipulation helpers.
 */
template <typename T> using GridTypeOf = typename T::GridType;
template <typename T> using IndexOf = typename T::CellIndex;
template <typename T> using CoordOf = typename T::CellCoord;
template <typename T> using SpaceAxisOf = typename T::SpaceAxis;
template <typename T> using PointOf = typename T::Point;

/**
 * @brief A cell grid type declarations.
 *
 * @tparam VDimension  number of grid dimensions, e.g., 3 for the real 3D world
 * @tparam TCellIndex  type of a single coordinate in cell space; also an absolute cell index
 * @tparam TSpaceAxis  type of a single coordinate in physical space
 */
template <unsigned int VDimension, std::integral TCellIndex = int,
          std::floating_point TSpaceAxis = double>
struct GridType
{
    constexpr static unsigned int Dimension =
        VDimension;               //!< number of dimmensions, e.g., 3 for the real 3D world
    using SpaceAxis = TSpaceAxis; //!< a single coordinate in physical space, e.g., double
    using CellIndex =
        TCellIndex; //!< a single coordinate in cell space, e.g., int; also an absolute cell index
    using Point = Eigen::Array<SpaceAxis, Dimension, 1>;     //!< a point in the VDimensional space
    using CellCoord = Eigen::Array<CellIndex, Dimension, 1>; //!< VDimensional cell coordinates
};

//! @brief Commonly used GridType in Faunus.
using Grid3DType = GridType<3, int, double>;

/**
 * @brief A simple data holder for common offsets in 3D grids, e.g., nearest neighbors or forward
 * neighbors.
 */
class GridOffsets3D
{
    using CellCoord = CoordOf<Grid3DType>; //!< 3-dimension cell coordinates
    using CellIndex = IndexOf<Grid3DType>; //!< cell index
  public:
    const std::vector<CellCoord> self{{0, 0, 0}}; //!< the cell itself
    std::vector<CellCoord> neighbors; //!< all neighbors (excluding self); (1 + 2×distance)³ - 1
    std::vector<CellCoord> forward_neighbors; //!< only the preceding half of all neighbors
    CellIndex distance; //!< maximal distance of neighbours (1 for the nearest neighbors)

    /**
     * @param distance  maximal distance of neighbours (1 for the nearest neighbors)
     */
    explicit GridOffsets3D(CellIndex distance);

  private:
    void initNeighbors();
    void initForwardNeighbors();
};
} // namespace Grid

namespace Container {
/**
 * Type manipulation helpers.
 */
template <typename T> using ContainerTypeOf = typename T::ContainerType;
template <typename T> using IndexOf = typename ContainerTypeOf<T>::Index;
template <typename T> using MemberOf = typename ContainerTypeOf<T>::Member;
template <typename T> using MembersOf = typename ContainerTypeOf<T>::Members;

/**
 * @brief An interface of a container.
 *
 * The container is responsible for storing (and retriving) members under a given cell index.
 * Various storage policies may be used in different implementations, e.g., a vector or a map.
 *
 * @tparam TMember  a member type
 * @tparam TIndex  an index type, see GridType::CellIndex
 */
template <typename TMember, typename TIndex> struct AbstractContainer
{
    using ContainerType = AbstractContainer<TMember, TIndex>;
    using Index = TIndex;                //!< the cell index (key) type
    using Member = TMember;              //!< the cell member (primitive value) type
    using Members = std::vector<Member>; //!< the cell members type

    //! @brief Gets cell members at the given cell index.
    virtual const Members& get(TIndex) const = 0;
    //! @brief Gets cell members at the given cell index.
    virtual Members& get(TIndex) = 0;
    //! @brief Gets cell members of an empty cell. Allows non-existing cells.
    virtual const Members& getEmpty() const = 0;
    //! @brief Gets indicies of non-empty cells.
    //! All non-empty cells are included. However, not all cells included must be non-empty.
    virtual std::vector<TIndex> indices() const = 0;
    virtual ~AbstractContainer() = default;
};
} // namespace Container

/**
 * Type manipulation helpers.
 */
template <typename T> using ContainerOf = typename T::Container;
template <typename T> using GridOf = typename T::Grid;
using Container::ContainerTypeOf;
using Grid::GridTypeOf;
template <typename T>
using IndexOf = typename Grid::IndexOf<GridTypeOf<typename T::AbstractCellList>>;
template <typename T>
using CoordOf = typename Grid::CoordOf<GridTypeOf<typename T::AbstractCellList>>;
template <typename T>
using PointOf = typename Grid::PointOf<GridTypeOf<typename T::AbstractCellList>>;
template <typename T>
using SpaceAxisOf = typename Grid::SpaceAxisOf<GridTypeOf<typename T::AbstractCellList>>;
template <typename T> using MemberOf = typename Container::MemberOf<typename T::AbstractCellList>;
template <typename T> using MembersOf = typename Container::MembersOf<typename T::AbstractCellList>;

/**
 * @brief An interface of a immutable cell list. No methods for manipulating (adding, removing)
 * members are provided.
 * @tparam TContainerType  a container type to obtain Member and Members subtypes
 * @tparam TGridType  a grid type to obtain index and coordinates subtypes
 */
template <class TContainerType, class TGridType> struct AbstractImmutableCellList
{
    using AbstractCellList = AbstractImmutableCellList<TContainerType, TGridType>; //! self
    using ContainerType = TContainerType;                         //!< abstract container types
    using GridType = TGridType;                                   //!< abstract grid type
    using Members = typename Container::MembersOf<ContainerType>; //!< members type
    using Member = typename Container::MemberOf<ContainerType>;   //!< member type
    using CellIndex = typename Grid::IndexOf<GridType>;           //!< grid (cell) index type
    using CellCoord = typename Grid::CoordOf<GridType>;           //!< grid (cell) coordinates type

    //! @brief Gets cell members at given cell coordinates.
    virtual const Members& getMembers(const CellCoord&) = 0;
    //! @brief Gets cell members at given cell coordinates.
    virtual const Members& getNeighborMembers(const CellCoord&, const CellCoord&) = 0;
    //! @brief Get cells count. Valid cell indices are in the range [0, gridSize).
    virtual CellIndex gridSize() const = 0; // todo check if needed
    //! @brief Gets indicies of non-empty cells. See AbstractContainer interface.
    //! @todo refactor as an iterator
    virtual std::vector<CellCoord> getCells() const = 0;
    virtual ~AbstractImmutableCellList() = default;
};

/**
 * @brief An interface to a cell list that can return cell's members in sorted order.
 *
 * Interface is used to efficiently compute differences between cell lists.
 *
 * @tparam TContainerType  a container type to obtain Member and Members subtypes
 * @tparam TGridType  a grid type to obtain index and coordinates subtypes
 */
template <class TContainerType, class TGridType>
struct AbstractSortableCellList
    : virtual public AbstractImmutableCellList<TContainerType, TGridType>
{
    using typename AbstractImmutableCellList<TContainerType, TGridType>::AbstractCellList;
    using typename AbstractCellList::CellCoord;
    using typename AbstractCellList::Members;

    //! @brief Gets sorted cell members at given cell coordinates
    virtual const Members& getSortedMembers(const CellCoord&) = 0;
    virtual const Members& getNeighborSortedMembers(const CellCoord&, const CellCoord&) = 0;
    virtual ~AbstractSortableCellList() = default;
};
} // namespace CellList
} // namespace Faunus
