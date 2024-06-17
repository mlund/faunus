#include <doctest/doctest.h>
#include "celllistimpl.h"
namespace Faunus {
namespace CellList {
namespace Grid {

GridOffsets3D::GridOffsets3D(CellIndex distance) : distance(distance) {
    initNeighbors();
    initForwardNeighbors();
}

void GridOffsets3D::initNeighbors() {
    const auto number_of_neighbors = std::pow(1 + 2 * distance, 3) - 1;
    neighbors.clear();
    neighbors.reserve(number_of_neighbors);
    for (auto i = -distance; i <= distance; ++i) {
        for (auto j = -distance; j <= distance; ++j) {
            for (auto k = -distance; k <= distance; ++k) {
                if (i == 0 && j == 0 && k == 0) {
                    continue;
                }
                neighbors.emplace_back(i, j, k);
            }
        }
    }
    assert(neighbors.size() == number_of_neighbors);
}

void GridOffsets3D::initForwardNeighbors() {
    const auto forward_neighbors_count = neighbors.size() / 2;
    forward_neighbors.clear();
    forward_neighbors.reserve(forward_neighbors_count);
    std::copy_if(neighbors.begin(), neighbors.end(), std::back_inserter(forward_neighbors), [](auto& offset) {
        return (offset > CellCoord{0, 0, 0}).all();
    });
    assert(forward_neighbors.size() == forward_neighbors_count);
}

/**
 * Tests
 */
TEST_CASE("Grid3DFixed") {
    using TestGridBase = Grid3DFixed;
    using CellCoord = TestGridBase::CellCoord;
    TestGridBase grid({20., 30., 40.}, 1.5);
    CellCoord end{(200 / 15), (300 / 15), (400 / 15)};
    CellCoord first{0, 0, 0};
    CellCoord last = end - CellCoord{1, 1, 1};
    CellCoord coord1{5, 6, 7};

    REQUIRE_EQ(grid.size(), end.prod());

    SUBCASE("Cell coordinates are limited to [origin, last]") {
        CHECK(grid.isCell(first));
        CHECK_FALSE(grid.isCell(first - CellCoord{1, 0, 0}));
        CHECK_FALSE(grid.isCell(first - CellCoord{0, 1, 0}));
        CHECK_FALSE(grid.isCell(first - CellCoord{0, 0, 1}));

        CHECK(grid.isCell(last));
        CHECK_FALSE(grid.isCell(last + CellCoord{1, 0, 0}));
        CHECK_FALSE(grid.isCell(last + CellCoord{0, 1, 0}));
        CHECK_FALSE(grid.isCell(last + CellCoord{0, 0, 1}));

        CHECK_FALSE(grid.isCell(end));

        CHECK(grid.isCell(coord1));
        CHECK(grid.isCell(2 * coord1));
        CHECK_FALSE(grid.isCell(4 * coord1));

        CHECK_EQ(grid.index(first), 0);
        CHECK_EQ(grid.index(last), grid.size() - 1);

        CHECK((grid.coordinates(grid.index(coord1)) == coord1).all());
        CHECK_NE(grid.index(coord1.reverse()), grid.index(coord1));

        SUBCASE("Offsets cannot wrap around") {
            CHECK(grid.isNeighborCell(first, {0, 0, 0}));
            CHECK(grid.isNeighborCell(first, {0, 0, 1}));
            CHECK_FALSE(grid.isNeighborCell(first, {0, 0, -1}));
            CHECK(grid.isNeighborCell(last, {0, 0, 0}));
            CHECK(grid.isNeighborCell(last, {0, 0, -1}));
            CHECK_FALSE(grid.isNeighborCell(last, {0, 0, 1}));
        }

        SUBCASE("Spatial coordinates has to be inside the box") {
            CHECK(grid.isCellAt({0., 0., 0.}));
            CHECK(grid.isCellAt({1., 1., 1.}));
            CHECK(grid.isCellAt({19., 29., 39.}));
            CHECK_FALSE(grid.isCellAt({20., 30., 40.}));

            CHECK((grid.coordinatesAt({1., 1., 1.}) == first).all());
            CHECK((grid.coordinatesAt({19., 29., 39.}) == last).all());

            CHECK((grid.coordinatesAt({3., 2.5, 3.}) == CellCoord{1, 1, 1}).all());
            CHECK((grid.coordinatesAt({3., 3., 3.}) == CellCoord{1, 2, 1}).all());
        }
    }
}

TEST_CASE("Grid3DPeriodic") {
    using TestGridBase = Grid3DPeriodic;
    using CellCoord = TestGridBase::CellCoord;
    TestGridBase grid({20., 30., 40.}, 1.5);
    CellCoord end{(200 / 15), (300 / 15), (400 / 15)};
    CellCoord first{0, 0, 0};
    CellCoord last = end - CellCoord{1, 1, 1};
    CellCoord coord1{5, 6, 7};

    REQUIRE_EQ(grid.size(), end.prod());

    SUBCASE("All cell coordinates are allowed") {
        CHECK(grid.isCell(first));
        CHECK(grid.isCell(first - CellCoord{1, 0, 0}));
        CHECK(grid.isCell(first - CellCoord{0, 1, 0}));
        CHECK(grid.isCell(first - CellCoord{0, 0, 1}));

        CHECK(grid.isCell(last));
        CHECK(grid.isCell(last + CellCoord{1, 0, 0}));
        CHECK(grid.isCell(last + CellCoord{0, 1, 0}));
        CHECK(grid.isCell(last + CellCoord{0, 0, 1}));

        CHECK(grid.isCell(end));

        CHECK(grid.isCell(coord1));
        CHECK(grid.isCell(2 * coord1));
        CHECK(grid.isCell(4 * coord1));

        CHECK_EQ(grid.index(first), 0);
        CHECK_EQ(grid.index(last), grid.size() - 1);

        CHECK((grid.coordinates(grid.index(coord1)) == coord1).all());
        CHECK_NE(grid.index(coord1.reverse()), grid.index(coord1));

        SUBCASE("Offsets wrap around, but not double around") {
            CHECK(grid.isNeighborCell(first, {0, 0, 0}));
            CHECK(grid.isNeighborCell(first, {0, 0, 1}));
            CHECK(grid.isNeighborCell(first, {0, 0, -1}));
            CHECK_FALSE(grid.isNeighborCell(first, {0, 0, -200 / 15}));
            CHECK(grid.isNeighborCell(last, {0, 0, 0}));
            CHECK(grid.isNeighborCell(last, {0, 0, -1}));
            CHECK(grid.isNeighborCell(last, {0, 0, 1}));
            CHECK(grid.isNeighborCell(first, {0, 0, 200 / 15}));
            CHECK_FALSE(grid.isNeighborCell(first, {0, 0, 1 + 200 / 15}));
        }

        SUBCASE("Spatial coordinates wrap around the box") {
            CHECK(grid.isCellAt({0., 0., 0.}));
            CHECK(grid.isCellAt({1., 1., 1.}));
            CHECK(grid.isCellAt({19., 29., 39.}));
            CHECK(grid.isCellAt({20., 30., 40.}));

            CHECK((grid.coordinatesAt({1., 1., 1.}) == first).all());
            CHECK((grid.coordinatesAt({19., 29., 39.}) == last).all());

            CHECK((grid.coordinatesAt({3., 2.5, 3.}) == CellCoord{1, 1, 1}).all());
            CHECK((grid.coordinatesAt({3., 3., 3.}) == CellCoord{1, 2, 1}).all());

            CHECK_EQ(grid.indexAt({0., 0., 0.}), grid.indexAt({20., 30., 40.}));
            CHECK_EQ(grid.indexAt({21., 10., 10.}), grid.indexAt({1., 10., 10.}));
        }
    }
}
} // namespace Grid

namespace Container {
TEST_CASE_TEMPLATE("Cell List Container", TContainer, Container<DenseContainer<int, int>>,
                   Container<SparseContainer<int, int>>) {
    const int max = 1000;
    const int cell_ndx = 20;
    TContainer container(max);

    SUBCASE("Empty cell") {
        CHECK_EQ(container.indices().size(), 0);
        CHECK_EQ(container.get(cell_ndx).size(), 0);

        SUBCASE("Add two elements to a single cell") {
            container.insert(1, cell_ndx);
            container.insert(2, cell_ndx);
            CHECK_EQ(container.indices().size(), 1);
            CHECK_EQ(container.get(cell_ndx).size(), 2);

            SUBCASE("Remove elements from the cell") {
                container.erase(1, cell_ndx);
                CHECK_EQ(container.get(cell_ndx).size(), 1);
            }

            SUBCASE("Remove non-existing element") {
                container.erase(-1, cell_ndx);
                CHECK_EQ(container.get(cell_ndx).size(), 2);
            }
        }
    }
}
} // namespace Container

TEST_CASE("CellListReverseMap") {
    CellListReverseMap<CellListType<int, Grid::Grid3DFixed>> cell_list({10., 10., 10.}, 1.);
    using CellCoord = GridOf<decltype(cell_list)>::CellCoord;
    const CellCoord cell{2, 3, 2};

    REQUIRE_EQ(cell_list.getMembers(cell).size(), 0);
    SUBCASE("when insert") {
        cell_list.addMember(10, cell);
        cell_list.addMember(11, cell);
        CHECK_EQ(cell_list.getMembers(cell).size(), 2);
        SUBCASE("when removed") {
            cell_list.removeMember(10);
            CHECK_EQ(cell_list.getMembers(cell).size(), 1);
        }
    }
}

TEST_CASE("CellListSpatial") {
    CellListSpatial<CellListType<int, Grid::Grid3DFixed>> cell_list({10., 10., 10.}, 1.);
    using CellCoord = GridOf<decltype(cell_list)>::CellCoord;
    using Point = GridOf<decltype(cell_list)>::Point;
    const Point pos{2.5, 3.5, 2.5};
    const CellCoord cell{2, 3, 2};

    REQUIRE_EQ(cell_list.getMembers(cell).size(), 0);
    SUBCASE("when insert") {
        cell_list.insertMember(10, pos);
        cell_list.insertMember(11, pos);
        CHECK_EQ(cell_list.getMembers(cell).size(), 2);
        SUBCASE("when removed") {
            cell_list.removeMember(10);
            CHECK_EQ(cell_list.getMembers(cell).size(), 1);
        }
    }
}

TEST_CASE("CellListDifference") {
    using CellList = CellListReverseMap<CellListType<int, Grid::Grid3DFixed, SortableCellList>>;
    CellList x({10., 10., 10.}, 1.);
    auto minuend = std::make_shared<CellList>(PointOf<CellList>{10., 10., 10.}, 1.);
    auto subtraend = std::make_shared<CellList>(PointOf<CellList>{10., 10., 10.}, 1.);
    using CellCoord = GridOf<CellList>::CellCoord;
    const CellCoord cell{2, 3, 4};
    const CellCoord other_cell{2, 3, 5};

    for (auto i : std::vector<int>{10, 8, 6, 12})
        minuend->addMember(i, cell);
    minuend->addMember(33, other_cell);

    REQUIRE_EQ(minuend->getMembers(cell).size(), 4);
    REQUIRE_EQ(minuend->getMembers(other_cell).size(), 1);

    subtraend->addMember(8, cell);
    subtraend->addMember(30, other_cell);
    REQUIRE_EQ(subtraend->getMembers(cell).size(), 1);
    REQUIRE_EQ(subtraend->getMembers(other_cell).size(), 1);

    CellListDifference diff(minuend, subtraend);

    auto cell_diff = diff.getMembers(cell);
    CHECK_EQ(cell_diff.size(), 3);
    WARN(std::is_sorted(cell_diff.begin(), cell_diff.end()));
    CHECK_EQ(diff.getMembers(other_cell).size(), 1);
}

} // namespace CellList
} // namespace Faunus
