 //-+--------------------------------------------------------------------
// Igatools a general purpose Isogeometric analysis library.
// Copyright (C) 2012-2014  by the igatools authors (see authors.txt).
//
// This file is part of the igatools library.
//
// The igatools library is free software: you can use it, redistribute
// it and/or modify it under the terms of the GNU General Public
// License as published by the Free Software Foundation, either
// version 3 of the License, or (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.
//-+--------------------------------------------------------------------

#ifndef GRID_TOOLS_H_
#define GRID_TOOLS_H_

#include <igatools/geometry/cartesian_grid.h>

IGA_NAMESPACE_OPEN

namespace grid_tools
{
/**
 * Given one grid <tt>grid_coarse</tt> and a refinement <tt>grid_fine</tt>,
 * this function builds and returns the one-to-one mapping between the elements on the
 * fine and the elements on the coarse.
 * @code{.cpp}
   fine_to_coarse_elements = build_map_elements_between_cartesian_grids(grid_fine,grid_coarse);
   // fine_to_coarse_elements[i] is the flat id of the element on the coarse grid that fully contains the
   // element on the fine grid with flat_id equal to i.
   @endcode
 *
 * @warning The grid must be defined on the same domain and each element on the fine grid must be
 * FULLY contained in one element of the coarse grid, otherwise an exception will be raised
 * (in Debug mode).
 *
 * @relates CartesianGrid
 */
template <int dim>
std::vector<Index>
build_map_elements_between_cartesian_grids(
    const CartesianGrid<dim> &grid_fine,
    const CartesianGrid<dim> &grid_coarse);


/**
 * Given two grids <tt>grid_1</tt> and <tt>grid_2</tt> defined over the same domain,
 * this function returns the grid that contains both. Moreover, this function gives back
 * the one-to-one mapping between the flat id of the elements in the CartesianGrid union with the
 * elements in the two starting grids.
 *
 * @param[in] grid_1 First grid.
 * @param[in] grid_2 Second grid.
 * @param[out] map_elem_grid_union_to_elem_grid_1 One-to-one mapping between the elements in the
 * CartesianGrid union and the elements in the first grid.
 * @param[out] map_elem_grid_union_to_elem_grid_2 One-to-one mapping between the elements in the
 * CartesianGrid union and the elements in the second grid.
 * @return The CartesianGrid that contains the all the knots from <tt>grid_1</tt> and <tt>grid_2</tt>.
 *
 * @warning In Debug mode an assertion will be raised if
 * the two grids <tt>grid_1</tt> and <tt>grid_2</tt> are not defined on the same domain.
 *
 * @relates CartesianGrid
 */
template <int dim>
std::shared_ptr<CartesianGrid<dim>>
build_cartesian_grid_union(
    const CartesianGrid<dim> &grid_1,
    const CartesianGrid<dim> &grid_2,
    std::vector<Index> &map_elem_grid_union_to_elem_grid_1,
    std::vector<Index> &map_elem_grid_union_to_elem_grid_2);
}

IGA_NAMESPACE_CLOSE

#endif /* GRID_TOOLS_H_ */
