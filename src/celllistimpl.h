#pragma once
/**
 * @brief Cell list class templates implementation.
 *
 * A minimal interface to be used in other header files is provided in "celllist.h" file.
 *
 * @author Richard Chudoba
 * @date 2020-02-01
 */

#include <vector>
#include <set>
#include <map>
#include <algorithm>
#include <cassert>
#include <Eigen/Core>
#include <range/v3/range/conversion.hpp>
#include <range/v3/view/join.hpp>
#include <range/v3/view/transform.hpp>
#include "celllist.h"

namespace Faunus {

/**
 * @namespace CellList
 * # Cell List
 *
 * A cell list allow to organize members, e.g., particles, into a grid cells hence allowing fast determination of their
 * approximate location in space. Cell lists can be used to quickly compute distances or interaction only between
 * close members. The granuality of distance is determined by the cell size and ussualy arise from some sort of cut off
 * distance.
 *
 * Each cell list is composed of a grid and a container by double inheritence. The grid (Grid) is responsible for
 * the cell indexing and the container (Container) holds he members assigned to respective cells. The class template
 * CellListBase is a bare-bone implementation with the the essential functionality only. It should be extended using
 * templated mixin to insert features as needed.
 *
 * ## Cell Grid
 *
 * The GridType defines basic data types, GridBase calculates the cell index from the cell coordinates and vice versa,
 * generally assuming fixed boundaries. Respective outer classes (FixedBoundaryGrid, PeriodicBoundaryGrid) wrapps
 * boundary conditions.
 *
 * The convenient template aliases GridBoundary, Grid3DBoundary, Grid3DFixed and Grid3DPeriodic are provided. A simple
 * structure GridOffsets3D contains commonly used neighbours' offsets.
 *
 * ## Containers
 *
 * Containers holds the members and allow a simple manipulation (insertion and deletion). A dense (using a vector) and
 * sparse (using a map) implementations are provided.
 */
namespace CellList {

/**
 * @namespace Grid
 */
namespace Grid {

/**
 * @brief Interface of a grid
 * @tparam TGridType  a cell grid types declaration
 */
template <typename TGridType> struct AbstractGrid : private TGridType {
    using GridType = TGridType;
    using typename GridType::CellCoord;
    using typename GridType::CellIndex;
    using typename GridType::Point;
    using typename GridType::SpaceAxis;

    virtual CellCoord coordinates(const CellIndex index) const = 0;
    virtual CellIndex index(const CellCoord& coordinates) const = 0;
    virtual CellCoord coordinatesAt(const Point& position) const = 0;
    virtual CellIndex indexAt(const Point& position) const = 0;
    virtual bool isCell(const CellCoord& coordinates) const = 0;
    virtual bool isCellAt(const Point& point) const = 0;
    virtual CellIndex size() const = 0;
    virtual ~AbstractGrid() = default;
};

template <typename TGridType> struct AbstractNeighborGrid : public virtual AbstractGrid<TGridType> {
    using GridType = TGridType;
    using typename TGridType::CellCoord;

    virtual bool isNeighborCell(const CellCoord& coordinates, const CellCoord& offset) const = 0;
    virtual ~AbstractNeighborGrid() = default;
};

/**
 * @brief A cell grid core class assuming fixed boundary conditions.
 *
 * The grid converts spatial coordinates to cell coordinates and cell coordinates to the scalar cell index. The grid
 * is designed for arbitrary number of dimensions.
 *
 * @tparam TGridType  a cell grid type declaration
 */
template <typename TGridType> class GridBase : public virtual AbstractGrid<TGridType> {
  public:
    using GridType = typename AbstractGrid<TGridType>::GridType;
    using typename GridType::CellCoord;
    using typename GridType::CellIndex;
    using typename GridType::Point;
    using typename GridType::SpaceAxis;

    /**
     * @brief Compute the cell coordinates from the cell index.
     *
     * ⌊ (index mod cell_index_to_coordinates) / cell_coordinates_to_index ⌋
     *
     * @param index
     * @return  cell coordinates
     */
    CellCoord coordinates(const CellIndex index) const override {
        assert(index < size());
        return cell_index_to_coordinates.unaryExpr([index](const auto modulo) { return index % modulo; }) /
               cell_coordinates_to_index;
    }

    /**
     * @brief Compute the cell index from the cell coordinates.
     * @param coordinates  cell coordinates
     * @return  cell index
     */
    CellIndex index(const CellCoord& coordinates) const override {
        assert(isCell(coordinates));
        return (coordinates * cell_coordinates_to_index).sum();
    }

    /**
     * @brief Compute cell coordinates from a spatial position.
     *
     * ⌊ position / cell_dimension ⌋
     *
     * @param position
     * @return  cell coordinates
     */
    CellCoord coordinatesAt(const Point& position) const override {
        assert(isCellAt(position));
        return (position / cell).floor().template cast<CellIndex>();
    }

    /**
     * @brief Compute the cell index from the spatial position.
     * @param position   spatial position
     * @return  cell index
     */
    CellIndex indexAt(const Point& position) const override { return index(coordinatesAt(position)); }

    /**
     * @brief Verifies if coordinates are valid in the grid.
     *
     * Coordinates must be within [0, cell_list_end).
     *
     * @param coordinates
     * @return  true if cell coordinates are valid, false otherwise
     */
    bool isCell(const CellCoord& coordinates) const override {
        return (coordinates < this->getCellListEnd()).all() && (coordinates >= CellCoord::Zero()).all();
    }

    /**
     * @brief Verifies if spatial coordinates are in the box.
     *
     * Coordinates must be within [0, box).
     *
     * @param point
     * @return  true if spatial coordinates are valid, false otherwise
     */
    bool isCellAt(const Point& point) const override {
        return (point < this->getBox()).all() && (point >= Point::Zero()).all();
    }

    /**
     * @return CellIndex  total number of cells
     */
    CellIndex size() const override { return cell_list_end.prod(); };

    /**
     * @param box  a box to be divided into cells, e.g., a simulation box
     * @param minimal_cell_dimension  minimal cell dimension, e.g., a cut off distance
     */
    GridBase(const Point& box, SpaceAxis minimal_cell_dimension) : minimal_cell_dimension(minimal_cell_dimension) {
        updateCellGrid(box);
    }

  protected:
    const CellCoord& getCellListEnd() const { return cell_list_end; }
    const Point& getBox() const { return box; }
    const Point& getCell() const { return cell; }

  private:
    CellCoord cell_list_end; //!< last adressable coordinates in each direction, i.e., (last_x+1, last_y+1, …)
    SpaceAxis minimal_cell_dimension; //!< minimal cell dimension, e.g., a cut off distance
    Point box;                        //!< dimensions of covered cubic box
    Point cell;                       //!< current dimensions of the cell

    //! transformation from CDIM-dimension coordinates space to 1D index space; e.g., [1, i, i×j, …]
    //! index = (coordinates × cell_coordinates_to_index).sum()
    CellCoord cell_coordinates_to_index;
    //! transformation from 1D index space to CDIM-dimension coordinates space, e.g. [i, i×j, i×j×k, …]
    //! coordinates = ⌊ (index mod cell_index_to_coordinates) / cell_coordinates_to_index ⌋
    CellCoord cell_index_to_coordinates;

    /**
     * @brief Setup immutable parameters based on the box dimensions. To be called in constructors only.
     * @param new_box
     */
    void updateCellGrid(const Point& new_box) {
        box = new_box;
        cell_list_end = (box / box.min(minimal_cell_dimension)).floor().template cast<CellIndex>();
        cell = box / cell_list_end.template cast<SpaceAxis>();

        // (x, y, z, …) -> (1, x, x*y, …)
        cell_coordinates_to_index[0] = 1;
        cell_index_to_coordinates[0] = cell_list_end[0];
        for (int i = 1; i < TGridType::Dimension; ++i) {
            cell_coordinates_to_index[i] = cell_list_end[i - 1] * cell_coordinates_to_index[i - 1];
            cell_index_to_coordinates[i] = cell_list_end[i] * cell_index_to_coordinates[i - 1];
        }
    }
};

/**
 * @brief Fixed boundary conditions.
 *
 * @tparam TGridType  a cell grid type declaration
 */
template <typename TGridType>
class FixedBoundaryGrid : public GridBase<TGridType>, virtual public AbstractNeighborGrid<TGridType> {
    using Base = GridBase<TGridType>;

  public:
    using typename AbstractNeighborGrid<TGridType>::GridType;
    using typename GridType::CellCoord;

    /**
     * @brief Verifies if offseted coordinates are valid.
     *
     * Coordinates must be within [0, cell_list_end).
     *
     * @param coordinates
     * @param offset
     * @return  true if cell coordinates are valid, false otherwise
     */
    bool isNeighborCell(const CellCoord& coordinates, const CellCoord& offset) const {
        return Base::isCell(coordinates + offset);
    }

    using Base::Base; // use parent constructors
};

/**
 * @brief Periodic boundary conditions.
 *
 * Both spatial coordinates and cell coordinates are first converted to the original image and then evaluated.
 * To avoid repeated cell inclusion by completing a full circle in periodic boundary conditions, offsets are limited
 * in both directions assuming symetry.
 *
 * @todo Currently only full or none PBC is supperted. To support geometries like a slit, more fine distingtion is
 * needed.
 *
 * @tparam TGridType  a cell grid type declaration
 */
template <typename TGridType>
class PeriodicBoundaryGrid : public GridBase<TGridType>, virtual public AbstractNeighborGrid<TGridType> {
    using Base = GridBase<TGridType>;

  public:
    using GridType = typename AbstractNeighborGrid<TGridType>::GridType;
    using typename GridType::CellCoord;
    using typename GridType::CellIndex;
    using typename GridType::Point;
    using typename GridType::SpaceAxis;

    /**
     * @brief Compute the cell coordinates from the cell index.
     *
     * ⌊ (index mod cell_index_to_coordinates) / cell_coordinates_to_index ⌋
     *
     * @param index
     * @return  cell coordinates
     */
    CellCoord coordinates(const CellIndex index) const override { return Base::coordinates(index); }

    /**
     * @brief Compute the cell index from the cell coordinates.
     *
     * coordinates mod cell_list_end
     *
     * @param coordinates  cell coordinates
     * @return  cell index
     */
    CellIndex index(const CellCoord& coordinates) const override {
        auto pbc_coordinates = coordinates;
        auto& boundary_coords = this->getCellListEnd();
        for (auto i = 0; i < pbc_coordinates.size(); ++i) {
            if (coordinates[i] < 0) {
                pbc_coordinates[i] += boundary_coords[i] * (std::abs(coordinates[i]) / boundary_coords[i] + 1);
            } else if (coordinates[i] >= boundary_coords[i]) {
                pbc_coordinates[i] -= boundary_coords[i] * (coordinates[i] / boundary_coords[i]);
            }
        }
        return Base::index(pbc_coordinates);
    }

    /**
     * @brief Compute cell coordinates from a spatial position.
     *
     * ⌊ position / cell_dimension ⌋
     *
     * @param position
     * @return  cell coordinates
     */
    CellCoord coordinatesAt(const Point& position) const override {
        assert(isCellAt(position));
        return (position / this->getCell()).floor().template cast<CellIndex>();
    }

    /**
     * @brief Compute the cell index from the spatial position.
     * @param position   spatial position
     * @return  cell index
     */
    CellIndex indexAt(const Point& position) const override { return this->index(this->coordinatesAt(position)); }

    /**
     * @brief Verifies if coordinates are valid in the grid.
     *
     * Any coordinates are valid under full PBC.
     *
     * @param coordinates
     * @return  always true
     */
    bool isCell(const CellCoord&) const override { return true; }

    /**
     * @brief Verifies if spatial coordinates are in the box.
     *
     * Any coordinates are valid under full PBC.
     *
     * @param point
     * @return  always true
     */
    bool isCellAt(const Point&) const override { return true; }

    /**
     * @brief Verifies if offseted coordinates are valid taking into account periodic boundary conditions.
     *
     * Symmetric offsets around the reference cell are assumed. To avoid repeated inclusion of a cell, the offset
     * is limited to be at most the half of the cell grid, with an attention paid to the odd and even number of cells.
     *
     * @param coordinates
     * @param offset
     * @return  true if cell coordinates are valid, false otherwise
     */
    bool isNeighborCell([[maybe_unused]] const CellCoord& coordinates, const CellCoord& offset) const override {
        return (offset >= neighbours_cell_offset_min).all() && (offset <= neighbours_cell_offset_max).all();
    }

    PeriodicBoundaryGrid(const Point& box, SpaceAxis cell_dimension) : Base(box, cell_dimension) { updateCellGrid(); }

  private:
    CellCoord neighbours_cell_offset_min; //!< utmost negative offset of neighbours to prevent circular match in PBC
    CellCoord neighbours_cell_offset_max; //!< utmost positive offset of neighbours to prevent circular match in PBC

    /**
     * @brief Setup immutable parameters based on the box dimensions. To be called in constructors only.
     * @param new_box
     */
    void updateCellGrid() {
        const auto cell_list_middle = (this->getCellListEnd() - 1).template cast<double>() * 0.5;
        neighbours_cell_offset_min = -cell_list_middle.floor().template cast<CellIndex>();
        neighbours_cell_offset_max = cell_list_middle.ceil().template cast<CellIndex>();
    }
};

/**
 * @brief Provides a cell grid with or without PBC.
 * @tparam PBC  if periodic boundary condition shall be employed
 */
template <bool PBC, class TGridType>
using GridBoundary =
    typename std::conditional<PBC, PeriodicBoundaryGrid<TGridType>, FixedBoundaryGrid<TGridType>>::type;
/**
 * @brief Provides a 3D cell grid.
 * @tparam PBC  if periodic boundary condition shall be employed
 */
template <bool PBC> using Grid3DBoundary = GridBoundary<PBC, Grid3DType>;
using Grid3DFixed = Grid3DBoundary<false>;
using Grid3DPeriodic = Grid3DBoundary<true>;
} // namespace Grid

/**
 * @namespace Container
 *
 * Containers store members in individual cells with a common minimalistic interface. DenseContainer (based on
 * a vector) and Sparse container (based on a map) are currently provided.
 */
namespace Container {
/**
 * @brief Container based on a vector suitable for system with a uniform and high member density.
 *
 * As all cells (even the empty ones) reside in memory, the vector has size of gridBase::size(). Very fast access at
 * the expense of a bigger memory footprint, expecialy in sparse systems.
 *
 * @tparam TMember  member (value)
 * @tparam TIndex  index (key); corresponds do the CellIndex of GridType
 */
template <typename TMember, typename TIndex> class DenseContainer : virtual public AbstractContainer<TMember, TIndex> {
  public:
    using ContainerType = AbstractContainer<TMember, TIndex>;
    using typename ContainerType::Index;
    using typename ContainerType::Members;

    const Members& get(Index index) const override {
        assert(index >= 0 && index < indexEnd());
        return container[index];
    }

    Members& get(Index index) override {
        assert(index >= 0 && index < indexEnd());
        return container[index];
    }

    const Members& getEmpty() const override { return empty_set; }

    /**
     * @brief Return all cell indices that may contain members.
     *
     * It is not guaranteed that all returned cells are non-empty though. However all non-empty cells are returned.
     *
     * @return  vector of cell indices
     */
    std::vector<Index> indices() const override {
        std::vector<Index> indices;
        indices.reserve(std::distance(container.begin(), container.end()));
        auto iterator = container.begin();
        while ((iterator = std::find_if(iterator, container.end(),
                                        [](const auto& members) { return !members.empty(); })) != container.end()) {
            indices.push_back(std::distance(iterator, container.begin()));
            std::advance(iterator, 1);
        }
        return indices;
    }

    DenseContainer(std::size_t size) { container.resize(size); }

  protected:
    const Index indexEnd() const { return container.size(); }

  private:
    const Members empty_set;        //!< an empty set, e.g., for out of boundary conditions
    std::vector<Members> container; //!< container itself
};

/**
 * @brief Container based on a map suitable for system with a non-uniform or low member density.
 *
 * Only cells with members reside in memory. Smaller memory footprint at expense of slower access.
 *
 * @tparam TMember  member (value)
 * @tparam TIndex  index (key); corresponds do CellIndex of GridBase
 */
template <typename TMember, typename TIndex> class SparseContainer : virtual public AbstractContainer<TMember, TIndex> {
  public:
    using ContainerType = AbstractContainer<TMember, TIndex>;
    using typename ContainerType::Index;
    using typename ContainerType::Members;

    const Members& get(Index index) const override {
        assert(index >= 0 && index < indexEnd());
        try {
            return container.at(index);
        } catch (std::out_of_range& e) { return empty_set; }
    }

    Members& get(Index index) override {
        assert(index >= 0 && index < indexEnd());
        try {
            return container.at(index);
        } catch (std::out_of_range& e) {
            [[maybe_unused]] auto [iterator, flag] = container.emplace(index, Members{});
            return iterator->second;
        }
    }

    const Members& getEmpty() const override { return empty_set; }

    /**
     * @brief Return all cell indicis that may contain members.
     *
     * It is not guaranteed that all selected cells are non-empty though.
     *
     * @return  vector of cell indices
     */
    std::vector<Index> indices() const override {
        std::vector<Index> indices;
        indices.reserve(std::distance(container.begin(), container.end()));
        std::transform(container.begin(), container.end(), back_inserter(indices),
                       [](const auto& pair) { return pair.first; });
        return indices;
    }

    SparseContainer(std::size_t size) : index_end(size) {}

  protected:
    const Index indexEnd() const { return index_end; }

  private:
    Index index_end;                    //!< the lowest index not allowed to appear, i.e., the grid size
    const Members empty_set;            //!< an empty set, e.g., for out of boundary conditions
    std::map<Index, Members> container; //!< container itself
};

/**
 * @brief A mixin allowing to insert and erase members into and from the container, respectively.s
 * @tparam TMember
 * @tparam TIndex
 * @tparam TContainer
 */
template <class TContainer> class Container : public TContainer {
  public:
    using ContainerType = ContainerTypeOf<TContainer>;
    using Member = MemberOf<ContainerType>;
    using Index = IndexOf<ContainerType>;

    void insert(const Member& member, const Index index_new) {
        assert(index_new < this->indexEnd());
        this->get(index_new).push_back(member);
    }

    void erase(const Member& member, const Index index_old) {
        assert(index_old < this->indexEnd());
        auto& members = this->get(index_old);
        // std::erase_if available in C++20
        auto remove_it = std::remove(members.begin(), members.end(), member);
        members.erase(remove_it, members.end());
        // end of std::erase_if
    }

    void move(const Member& member, const Index index_old, const Index index_new) {
        erase(member, index_old);
        insert(member, index_new);
    }

    using TContainer::TContainer; // use parent constructor
};

} // namespace Container

/**
 * @brief A basic class for a cell list allowing to keep members organized into a grid.
 *
 * It provides the minimal functionality that should be extended by mixin and decorators as needed.
 *
 * @tparam TMember
 * @tparam TGrid
 */
template <class TContainer, class TGrid>
class CellListBase : protected TGrid,
                     protected TContainer,
                     virtual public AbstractImmutableCellList<ContainerTypeOf<TContainer>, GridTypeOf<TGrid>> {
  public:
    using AbstractCellList = AbstractImmutableCellList<ContainerTypeOf<TContainer>, GridTypeOf<TGrid>>;
    using Grid = TGrid;
    using Container = TContainer;
    using typename AbstractCellList::CellCoord; //!< grid (cell) coordinates type
    using typename AbstractCellList::CellIndex; //!< grid (cell) index type
    using typename AbstractCellList::GridType;  //!< grid type
    using typename AbstractCellList::Member;    //!< member type
    using typename AbstractCellList::Members;   //!< members type

    /**
     * @brief Get cell members at given cell coordinates.
     * @param cell_coordinates
     * @return  vector of cell members
     */
    const Members& getMembers(const CellCoord& cell_coordinates) override {
        return this->get(this->index(cell_coordinates));
    }

    /**
     * @brief Get cell members at given offseted cell coordinates, or empty if offseted dcoordinates do not exist.
     * @param cell_coordinates
     * @return  vector of cell members or empty vector if cell does not exist
     */
    const Members& getNeighborMembers(const CellCoord& cell_coordinates, const CellCoord& cell_offset) override {
        if (this->isNeighborCell(cell_coordinates, cell_offset)) {
            return getMembers(cell_coordinates + cell_offset);
        } else {
            return getEmpty();
        }
    }

    /**
     * @return  range of containg coordinates of all possibly non-empty cells
     */
    std::vector<CellCoord> getCells() const override {
        const auto indices = this->indices();
        return ranges::cpp20::views::all(indices) |
               ranges::cpp20::views::transform([this](auto index) { return this->coordinates(index); }) |
               ranges::to_vector;
    }

    /**
     * @return  grid size as the total number of cells
     */
    CellIndex gridSize() const override { return Grid::size(); }

    /**
     * @brief Static cast itself to a grid.
     *
     * Useful to avoid using a copy constructor.
     *
     * @return  self statically cast as a grid`
     */
    const TGrid& getGrid() { return *this; }

    /**
     *  @brief Construct an empty cell list from a grid.
     *  @param cell_grid
     */
    CellListBase(const TGrid& cell_grid) : Grid(cell_grid), Container(cell_grid.size()) {}

    /**
     *  @brief Construct an empty cell list knowing a box dimension and minimal cell dimension.
     *  @param box  spatial dimension of the box
     *  @param cell_dimension  minimal length of the cell's edge
     */
    CellListBase(const typename TGrid::Point& box, const typename TGrid::SpaceAxis cell_dimension)
        : Grid(box, cell_dimension), Container(Grid::size()) {}

  protected:
    //    const auto& getMembers(const CellIndex cell_index) const { return this->get(cell_index); }
    //    auto& getMembers(const CellIndex cell_index) { return this->get(cell_index); }
    //    const auto& getEmptyCell() const { return this->getEmpty(); }

    using Container::erase;
    using Container::get;
    using Container::getEmpty;
    using Container::insert;
    using Container::move;
};

/**
 * @brief A template type alias for CellList.
 */
template <typename TMember, class TGrid, template <typename, typename> class TCellList = CellListBase,
          template <typename, typename> class TContainer = Container::DenseContainer>
using CellListType = TCellList<Container::Container<TContainer<TMember, Grid::IndexOf<TGrid>>>, TGrid>;

/**
 * @brief A mixing returning members in a cell in sorted order, e.g., for lookup or filtering.
 *
 * Sorting is performed on as-needed basis, i.e., only when requested and if modified since the last request.
 *
 * @tparam TBase
 */
template <class TContainer, class TGrid>
class SortableCellList : public CellListBase<TContainer, TGrid>,
                         virtual public AbstractSortableCellList<ContainerTypeOf<TContainer>, GridTypeOf<TGrid>> {
    using TBase = CellListBase<TContainer, TGrid>;

  public:
    using AbstractCellList =
        typename AbstractSortableCellList<ContainerTypeOf<TContainer>, GridTypeOf<TGrid>>::AbstractCellList;
    using Grid = TGrid;
    using Container = TContainer;
    using typename AbstractCellList::CellCoord; //!< grid (cell) coordinates type
    using typename AbstractCellList::CellIndex; //!< grid (cell) index type
    using typename AbstractCellList::Member;    //!< member type
    using typename AbstractCellList::Members;   //!< members type

    /**
     * @brief Get members in a sorted order.
     * @param cell_coordinates
     * @return  sorted members
     */
    const Members& getSortedMembers(const CellCoord& cell_coordinates) override {
        return getSorted(this->index(cell_coordinates));
    }

    /**
     * @brief Get members of a neighbor cell in a sorted order.
     * @param cell_coordinates
     * @param cell_offset
     * @return  sorted members
     */
    const Members& getNeighborSortedMembers(const CellCoord& cell_coordinates, const CellCoord& cell_offset) override {
        if (this->isNeighborCell(cell_coordinates, cell_offset)) {
            return getSortedMembers(cell_coordinates + cell_offset);
        } else {
            return this->getEmpty();
        }
    }

    /**
     *  @brief Construct an empty cell list from a grid.
     *  @param cell_grid
     */
    SortableCellList(const GridOf<TBase>& cell_grid) : TBase(cell_grid) { is_sorted.resize(this->size(), false); }

    /**
     *  @brief Construct an empty cell list knowing a box dimension and minimal cell dimension.
     *  @param box  spatial dimension of the box
     *  @param cell_dimension  minimal length of the cell's edge
     */
    SortableCellList(const typename GridOf<TBase>::Point& box, const typename GridOf<TBase>::SpaceAxis cell_dimension)
        : TBase(box, cell_dimension) {
        is_sorted.resize(this->size(), false);
    }

  protected:
    void insert(const Member& member, const CellIndex& new_cell_index) {
        TBase::insert(member, new_cell_index);
        if (is_sorted[new_cell_index]) {
            is_sorted[new_cell_index] = false;
        }
    }

    void remove(const Member& member, const CellIndex& old_cell_index) {
        TBase::remove(member, old_cell_index);
        if (is_sorted[old_cell_index]) {
            is_sorted[old_cell_index] = false;
        }
    }

    /**
     * @brief Get members in a sorted order.
     * @param cell_index
     * @return  sorted members
     */
    const Members& getSorted(const CellIndex cell_index) {
        auto& members = this->get(cell_index);
        if (!is_sorted[cell_index]) {
            std::sort(members.begin(), members.end());
            is_sorted[cell_index] = true;
        }
        assert(std::is_sorted(members.begin(), members.end()));
        return members;
    }

  private:
    // TODO perhaps a map for sparse containers
    std::vector<bool> is_sorted;
};

/**
 * @brief Reverse mapping (from particle to its cell) is stored.  Particles can be thus removed or update without
 * providing their previous cell.
 *
 * @tparam TBase
 */
template <class TBase> class CellListReverseMap : public TBase {
  public:
    using typename TBase::AbstractCellList;     //!< a structure with types definitions
    using typename AbstractCellList::CellCoord; //!< grid (cell) coordinates type
    using typename AbstractCellList::CellIndex; //!< grid (cell) index type
    using typename AbstractCellList::Member;    //!< member type
    using typename TBase::Container;
    using typename TBase::Grid;

    void addMember(const Member& member, const CellCoord& new_cell_coordinates) {
        this->insert(member, this->index(new_cell_coordinates));
    }

    void removeMember(const Member& member) {
        const auto old_cell_index = member2cell.at(member);
        this->erase(member, old_cell_index);
        member2cell.erase(member);
    }

    void updateMember(const Member& member, const CellCoord& new_cell_coordinates) {
        this->update(member, this->index(new_cell_coordinates));
    }

    /**
     * @brief returns true if member is present in the cell list false if not
     * @param member
     */
    bool containsMember(const Member& member) { return member2cell.count(member) != 0; }

    /**
     * @brief Imports members from other list without computing cell coordinates from member positions.
     * @tparam T
     * @param source
     * @param members
     */
    template <typename T> void importMembers(CellListReverseMap& source, const T& members) {
        for (const auto& member : members) {
            insert(member, source.member2cell[member]);
        }
    }

    using TBase::TBase;

  protected:
    void insert(const Member& member, const CellIndex& new_cell_index) {
        assert(member2cell.count(member) == 0); // FIXME C++20 contains
        TBase::insert(member, new_cell_index);
        member2cell.insert({member, new_cell_index});
    }

    void update(const Member& member, const CellIndex& new_cell_index) {
        const auto old_cell_index = member2cell.at(member);
        if (new_cell_index != old_cell_index) {
            TBase::move(member, old_cell_index, new_cell_index);
            member2cell.at(member) = new_cell_index;
        }
    }

  private:
    // TODO as a vector because size of Particle vector is constant
    std::map<Member, CellIndex> member2cell; //!< mapping from a member to its cell
};

/**
 * @brief A mixin using spatial coordinates instead of cell coordinates when inserting or removing a member.
 *
 * @tparam TBase
 */
template <class TBase> class CellListSpatial : public CellListReverseMap<TBase> {
  public:
    using typename CellListReverseMap<TBase>::AbstractCellList; //!< a structure with types definitions
    using typename CellListReverseMap<TBase>::Grid;
    using typename CellListReverseMap<TBase>::Container;
    using typename AbstractCellList::Member;          //!< member type
    using typename AbstractCellList::GridType::Point; //!< point

    void insertMember(const Member& member, const Point& new_position) {
        const auto new_cell_index = this->indexAt(new_position);
        this->insert(member, new_cell_index);
    }

    void updateMemberAt(const Member& member, const Point& new_position) {
        const auto new_cell_index = this->indexAt(new_position);
        this->update(member, new_cell_index);
    }

    using CellListReverseMap<TBase>::CellListReverseMap;
};

/**
 * @brief Wrapper of two cell list to efficiently provide difference of members, i.e., members presented in the minuend
 * and not in the subtrahend.
 *
 * @todo Allow cache invalidation when underlying cell lists change
 *
 * @tparam TCellListMinuend
 * @tparam TCellListSubtrahend
 */
template <typename TCellListMinuend, typename TCellListSubtrahend>
class CellListDifference : virtual public AbstractSortableCellList<ContainerTypeOf<ContainerOf<TCellListMinuend>>,
                                                                   GridTypeOf<GridOf<TCellListMinuend>>> {
    using CellListMinuend = TCellListMinuend;
    using CellListSubtrahend = TCellListSubtrahend;

  public:
    using AbstractCellList = typename CellListMinuend::AbstractCellList;
    using Grid = typename CellListMinuend::Grid;
    using Container = typename CellListMinuend::Container;
    using typename AbstractCellList::CellCoord; //!< grid (cell) coordinates type
    using typename AbstractCellList::CellIndex; //!< grid (cell) index type
    using typename AbstractCellList::Member;    //!< member type
    using typename AbstractCellList::Members;   //!< members type

    const Members& getSortedMembers(const CellCoord& cell_coordinates) override {
        const auto subtrahend_cell = subtrahend->getSortedMembers(cell_coordinates);
        if (std::size(subtrahend_cell) == 0) {
            // nothing to subtract
            return minuend->getSortedMembers(cell_coordinates);
        } else {
            const auto cell_index = minuend->getGrid().index(cell_coordinates);
            auto difference_cell_it = difference_cache.find(cell_index);
            if (difference_cell_it == difference_cache.end()) {
                // not cached yet
                auto minuend_cell = minuend->getSortedMembers(cell_coordinates);
                std::tie(difference_cell_it, std::ignore) = difference_cache.emplace(cell_index, Members{});
                std::set_difference(minuend_cell.begin(), minuend_cell.end(), subtrahend_cell.begin(),
                                    subtrahend_cell.end(), std::back_inserter(difference_cell_it->second));
            }
            return difference_cell_it->second;
        }
    }

    const Members& getMembers(const CellCoord& cell_coordinates) override {
        if (subtrahend->getMembers(cell_coordinates).empty()) {
            // avoid difference and soring when not necessary
            return minuend->getMembers(cell_coordinates);
        } else {
            return getSortedMembers(cell_coordinates);
        }
    }

    const Members& getNeighborSortedMembers(const CellCoord& cell_coordinates, const CellCoord& cell_offset) override {
        if (minuend->getGrid().isNeighborCell(cell_coordinates, cell_offset)) {
            // avoid difference and soring when not necessary
            return getSortedMembers(cell_coordinates + cell_offset);
        } else {
            return minuend->getNeighborMembers(cell_coordinates, cell_offset);
        }
    }

    const Members& getNeighborMembers(const CellCoord& cell_coordinates, const CellCoord& cell_offset) override {
        if (minuend->getGrid().isNeighborCell(cell_coordinates, cell_offset)) {
            return getMembers(cell_coordinates + cell_offset);
        } else {
            return minuend->getNeighborMembers(cell_coordinates, cell_offset);
        }
    }

    CellIndex gridSize() const override { return minuend->gridSize(); }

    auto& getGrid() { return minuend->getGrid(); }

    std::vector<CellCoord> getCells() const override { return minuend->getCells(); }

    CellListDifference(std::shared_ptr<TCellListMinuend> minuend, std::shared_ptr<TCellListSubtrahend> subtrahend)
        : minuend(minuend), subtrahend(subtrahend) {}

  private:
    std::shared_ptr<CellListMinuend> minuend;
    std::shared_ptr<CellListSubtrahend> subtrahend;
    std::map<CellIndex, Members> difference_cache;
};

} // namespace CellList
} // namespace Faunus
