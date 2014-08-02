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

#include <igatools/geometry/grid_tools.h>
#include <set>

using std::vector;
using std::array;
using std::set;

IGA_NAMESPACE_OPEN

namespace grid_tools
{
template <int dim>
vector<Index>
build_map_elements_between_cartesian_grids(
    const CartesianGrid<dim> &grid_fine,
    const CartesianGrid<dim> &grid_coarse)
{
    //---------------------------------------------------------
    // checks that the grid are on the same domain
    Assert(grid_fine.get_bounding_box() == grid_coarse.get_bounding_box(),
           ExcMessage("Grids on different domains."));
    //---------------------------------------------------------

    //---------------------------------------------------------
    array<vector<int>,dim> map_interv_fid_fine_coarse;
    for (int i = 0 ; i < dim ; ++i)
    {
        const auto &coords_coarse = grid_coarse.get_knot_coordinates(i);
        const auto &coords_fine = grid_fine.get_knot_coordinates(i);

#ifndef NDEBUG
        const int n_intervals_coarse = coords_coarse.size() - 1;
#endif
        const int n_intervals_fine = coords_fine.size() - 1;

        for (int fid_fine = 0 ; fid_fine < n_intervals_fine ; ++fid_fine)
        {
            int fid_coarse = 0;
            while (!(coords_fine[fid_fine] >= coords_coarse[fid_coarse] &&
                     coords_fine[fid_fine+1] <= coords_coarse[fid_coarse+1]))
            {
                ++fid_coarse;
            }
            Assert(fid_coarse < n_intervals_coarse,
                   ExcMessage("Impossible to find an interval "
                              "on the coarse grid that fully contains the interval " +
                              std::to_string(fid_fine) + " along the direction "
                              + std::to_string(i) +
                              "of the fine grid."));

            map_interv_fid_fine_coarse[i].push_back(fid_coarse);
        }
    }

    const int n_elems_fine = grid_fine.get_num_active_elems();
    vector<int> map_elem_fine_to_elem_coarse(n_elems_fine);
    for (int elem_fine_fid = 0 ; elem_fine_fid < n_elems_fine ; ++elem_fine_fid)
    {
        TensorIndex<dim> elem_fine_tid =
            grid_fine.flat_to_tensor_element_index(elem_fine_fid);

        TensorIndex<dim> elem_coarse_tid;
        for (int i = 0 ; i < dim ; ++i)
            elem_coarse_tid[i] = map_interv_fid_fine_coarse[i][elem_fine_tid[i]];

        const int elem_coarse_fid =
            grid_coarse.tensor_to_flat_element_index(elem_coarse_tid);

        map_elem_fine_to_elem_coarse[elem_fine_fid] = elem_coarse_fid;
    }

    return map_elem_fine_to_elem_coarse;
}



template <int dim>
std::shared_ptr<CartesianGrid<dim>>
                                 build_cartesian_grid_union(
                                     const CartesianGrid<dim> &grid_1,
                                     const CartesianGrid<dim> &grid_2,
                                     vector<Index> &map_elem_grid_union_to_elem_grid_1,
                                     vector<Index> &map_elem_grid_union_to_elem_grid_2)
{
    //---------------------------------------------------------
    // checks that the grid are on the same domain
    Assert(grid_1.get_bounding_box() == grid_2.get_bounding_box(),
           ExcMessage("Grids on different domains."));
    //---------------------------------------------------------


    //---------------------------------------------------------
    // getting the coordinates from the two grids and building the grid union
    array<vector<Real>,dim> knots_union;
    for (int i = 0 ; i < dim ; ++i)
    {
        const auto &coords_grid_1 = grid_1.get_knot_coordinates(i);
        const auto &coords_grid_2 = grid_2.get_knot_coordinates(i);

        // here we remove the duplicates (if any)
        set<Real> coords_unique(coords_grid_1.begin(),coords_grid_1.end());
        std::copy(coords_grid_2.begin(), coords_grid_2.end(),
                  std::inserter(coords_unique, coords_unique.end()));

        std::copy(coords_unique.begin(),coords_unique.end(),
                  std::back_inserter(knots_union[i]));
    }
    auto grid_union = CartesianGrid<dim>::create(knots_union);
    //---------------------------------------------------------


    //---------------------------------------------------------
    map_elem_grid_union_to_elem_grid_1 =
        build_map_elements_between_cartesian_grids(*grid_union,grid_1);
    map_elem_grid_union_to_elem_grid_2 =
        build_map_elements_between_cartesian_grids(*grid_union,grid_2);
    //---------------------------------------------------------

    return grid_union;
}

}

IGA_NAMESPACE_CLOSE

#include <igatools/geometry/grid_tools.inst>
