//-+--------------------------------------------------------------------
// Igatools a general purpose Isogeometric analysis library.
// Copyright (C) 2012-2015  by the igatools authors (see authors.txt).
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

#ifndef __GRID_TOOLS_H_
#define __GRID_TOOLS_H_

#include <igatools/geometry/cartesian_grid.h>
IGA_NAMESPACE_OPEN

namespace grid_tools
{

/**
 * Returns true if all the values in <tt>knots_coarse</tt> are also present
 * in <tt>knot_fine</tt>.
 */
bool test_if_knots_fine_contains_knots_coarse(
  const SafeSTLVector<Real> &knots_fine,
  const SafeSTLVector<Real> &knots_coarse);

/**
 * Given one grid <tt>grid_coarse</tt> and a refinement <tt>grid_fine</tt>,
 * this function builds and returns the one-to-one mapping between the 1D intervals on the
 * fine grid and the 1D intervals on the coarse grid, for each coordinate direction.
 * @code{.cpp}
   fine_to_coarse_intervals = build_map_intervals_between_cartesian_grids(grid_fine,grid_coarse);
   // fine_to_coarse_intervals[dir][i] is the id of the interval along the dir direction
   // of the coarse grid that fully contains the
   // interval along the dir direction of the fine grid with flat_id equal to i.
   @endcode
 *
 * @warning The grid must be defined on the same domain and each element on the fine grid must be
 * FULLY contained in one element of the coarse grid, otherwise an exception will be raised
 * (in Debug mode).
 *
 * @relates CartesianGrid
 */
template <int dim>
SafeSTLArray<SafeSTLVector<Index>,dim>
build_map_intervals_id_between_cartesian_grids(const CartesianGrid<dim> &grid_fine,
                                               const CartesianGrid<dim> &grid_coarse);

#if 0
/**
 * This function returns the maximum numbers of intervals in the <tt>fine_grid</tt>
 * that are fully contained in an interval of the <tt>coarse_grid</tt>.
 *
 * @warning The grid must be defined on the same domain and each element on the fine grid must be
 * FULLY contained in one element of the coarse grid, otherwise an exception will be raised
 * (in Debug mode).
 *
 * @relates CartesianGrid
 */
template <int dim>
SafeSTLArray<Index,dim>
get_max_num_fine_intervals_in_coarse_interval(const CartesianGrid<dim> &grid_fine,
                                              const CartesianGrid<dim> &grid_coarse);


#endif

/**
 * Type alias for the container that associates element indices between two CartesianGrid.
 *
 * @relates CartesianGrid
 */
template<int dim>
using InterGridMap = std::map<typename CartesianGrid<dim>::IndexType,typename CartesianGrid<dim>::IndexType>;


/**
 * Given one CartesianGrid <tt>grid_coarse</tt> and a refinement <tt>grid_fine</tt>,
 * this function builds and returns the one-to-one mapping between the elements on the
 * fine and the elements on the coarse.
 * @code{.cpp}
   fine_to_coarse_elements_id = build_map_elements_id_between_cartesian_grids(grid_fine,grid_coarse);
   // fine_to_coarse_elements_id[i] is the flat id of the element on the coarse grid that fully contains the
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
InterGridMap<dim>
build_map_elements_id_between_cartesian_grids(const CartesianGrid<dim> &grid_fine,
                                              const CartesianGrid<dim> &grid_coarse);

#if 0

template<int dim>
using InterGridMap = std::map<typename CartesianGrid<dim>::ElementConstIterator,
      typename CartesianGrid<dim>::ElementConstIterator>;

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
InterGridMap<dim>
build_map_elements_between_cartesian_grids(const CartesianGrid<dim> &grid_fine,
                                           const CartesianGrid<dim> &grid_coarse);
#endif




/**
 * Given two CartesianGrid <tt>grid_1</tt> and <tt>grid_2</tt> defined over the same domain,
 * this function returns the grid that contains both. Moreover, this function gives back
 * the one-to-one mapping between the elements in the CartesianGrid union with the
 * elements in the two starting grids.
 *
 * @param[in] grid_1 First grid.
 * @param[in] grid_2 Second grid.
 * @param[out] map_elem_id_grid_union_to_elem_id_grid_1 One-to-one mapping between the elements id in the
 * CartesianGrid union and the elements id in the first grid.
 * @param[out] map_elem_id_grid_union_to_elem_id_grid_2 One-to-one mapping between the elements id in the
 * CartesianGrid union and the elements id in the second grid.
 * @return The CartesianGrid that contains the all the knots from <tt>grid_1</tt> and <tt>grid_2</tt>.
 *
 * @warning In Debug mode an assertion will be raised if
 * the two grids <tt>grid_1</tt> and <tt>grid_2</tt> are not defined on the same domain.
 *
 * @relates CartesianGrid
 */
template <int dim>
std::shared_ptr<CartesianGrid<dim> >
build_cartesian_grid_union(
  const CartesianGrid<dim> &grid_1,
  const CartesianGrid<dim> &grid_2,
  InterGridMap<dim> &map_elem_grid_union_to_elem_grid_1,
  InterGridMap<dim> &map_elem_grid_union_to_elem_grid_2);

}

IGA_NAMESPACE_CLOSE

#endif /* GRID_TOOLS_H_ */
