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

#include <igatools/geometry/grid_tools.h>
#include <set>



using std::set;

IGA_NAMESPACE_OPEN


namespace grid_tools
{
bool test_if_knots_fine_contains_knots_coarse(
  const SafeSTLVector<Real> &knots_fine,
  const SafeSTLVector<Real> &knots_coarse)
{
  Assert(std::is_sorted(knots_fine.begin(),knots_fine.end()),ExcMessage("Vector not sorted."));

  return std::all_of(
           knots_coarse.begin(),
           knots_coarse.end(),
           [&](const Real &knot_coarse)
  {
    return std::binary_search(knots_fine.begin(),knots_fine.end(),knot_coarse);
  }
         );
}

template <int dim>
SafeSTLArray<SafeSTLVector<Index>,dim>
build_map_intervals_id_between_grids(const Grid<dim> &grid_fine,
                                     const Grid<dim> &grid_coarse)
{
  //---------------------------------------------------------
  // checks that the grid are on the same domain
  Assert(grid_fine.get_bounding_box() == grid_coarse.get_bounding_box(),
         ExcMessage("Grids on different domains."));
  //---------------------------------------------------------

  //---------------------------------------------------------
  SafeSTLArray<SafeSTLVector<Index>,dim> map_interv_id_fine_coarse;
  for (int i = 0 ; i < dim ; ++i)
  {
    const auto &coords_coarse = grid_coarse.get_knot_coordinates(i);
    const auto &coords_fine = grid_fine.get_knot_coordinates(i);

    Assert(test_if_knots_fine_contains_knots_coarse(coords_fine,coords_coarse),
           ExcMessage("The knots of the fine grid does not contains the knots of the coarse grid"
                      " along the direction " + std::to_string(i)));

    const int n_intervals_fine = coords_fine.size() - 1;

    int fid_coarse = 0;
    for (int fid_fine = 0 ; fid_fine < n_intervals_fine ; ++fid_fine)
    {
      while (!(coords_fine[fid_fine] >= coords_coarse[fid_coarse] &&
               coords_fine[fid_fine+1] <= coords_coarse[fid_coarse+1]))
      {
        ++fid_coarse;
      }

      map_interv_id_fine_coarse[i].push_back(fid_coarse);
    }
  }
  return map_interv_id_fine_coarse;
}

#if 0
template <int dim>
SafeSTLArray<Index,dim>
get_max_num_fine_intervals_in_coarse_interval(const Grid<dim> &grid_fine,
                                              const Grid<dim> &grid_coarse)
{
  const auto map_interv_id_fine_coarse =
    build_map_intervals_id_between_grids(grid_fine,grid_coarse);

  //---------------------------------------------------------
  SafeSTLArray<Index,dim> max_num_fine_intervals_in_coarse_interval(0);
  for (int i = 0 ; i < dim ; ++i)
  {
    const int n_intervals_coarse = grid_coarse.get_knot_coordinates(i).size() - 1;

    for (int fid_coarse = 0 ; fid_coarse < n_intervals_coarse ; ++fid_coarse)
    {
      const int n_fine_intervals_this_coarse_interval =
        std::count(map_interv_id_fine_coarse[i].begin(),
                   map_interv_id_fine_coarse[i].end(),
                   fid_coarse);
      max_num_fine_intervals_in_coarse_interval[i] =
        std::max(max_num_fine_intervals_in_coarse_interval[i],n_fine_intervals_this_coarse_interval);
    }
  } // end loop i

  return max_num_fine_intervals_in_coarse_interval;
}
#endif


template <int dim>
InterGridMap<dim>
build_map_elements_id_between_grids(const Grid<dim> &grid_fine,
                                    const Grid<dim> &grid_coarse)
{
  const auto map_interv_id_fine_coarse =
    build_map_intervals_id_between_grids(grid_fine,grid_coarse);

  InterGridMap<dim> res;
  for (const auto &elem_fine : grid_fine)
  {
    const auto elem_fine_tid = elem_fine.get_index();

    TensorIndex<dim> elem_coarse_tid;
    for (int i = 0 ; i < dim ; ++i)
      elem_coarse_tid[i] = map_interv_id_fine_coarse[i][elem_fine_tid[i]];

    res.emplace(elem_fine_tid, elem_coarse_tid);
  }

  return res;
}

#if 0
template <int dim>
InterGridMap<dim>
build_map_elements_between_grids(const Grid<dim> &grid_fine,
                                 const Grid<dim> &grid_coarse)
{
  const auto map_interv_id_fine_coarse =
    build_map_intervals_id_between_grids(grid_fine,grid_coarse);

  InterGridMap<dim> res;
  const int n_elems_fine = grid_fine.get_num_all_elems();

  auto f_elem = grid_fine.begin();
  auto c_elem = grid_coarse.begin();
  for (int elem_fine_fid = 0 ; elem_fine_fid < n_elems_fine ; ++elem_fine_fid)
  {
    f_elem.move_to(elem_fine_fid);

    TensorIndex<dim> elem_coarse_tid;
    for (int i = 0 ; i < dim ; ++i)
      elem_coarse_tid[i] = map_interv_id_fine_coarse[i][f_elem.get_tensor_index()[i]];

    c_elem.move_to(elem_coarse_tid);

    res.emplace(f_elem, c_elem);
  }

  return res;
}
#endif


template <int dim>
std::shared_ptr<Grid<dim> >
build_grid_union(
  const Grid<dim> &grid_1,
  const Grid<dim> &grid_2,
  InterGridMap<dim> &map_elem_id_grid_union_to_elem_id_grid_1,
  InterGridMap<dim> &map_elem_id_grid_union_to_elem_id_grid_2)
{
  //---------------------------------------------------------
  // checks that the grid are on the same domain
  Assert(grid_1.get_bounding_box() == grid_2.get_bounding_box(),
         ExcMessage("Grids on different domains."));
  //---------------------------------------------------------


  //---------------------------------------------------------
  // getting the coordinates from the two grids and building the grid union
  SafeSTLArray<SafeSTLVector<Real>,dim> knots_union;
  for (int i = 0 ; i < dim ; ++i)
  {
    const auto &coords_grid_1 = grid_1.get_knot_coordinates(i);
    const auto &coords_grid_2 = grid_2.get_knot_coordinates(i);

    auto &coords_grid_union = knots_union[i];
    coords_grid_union.resize(coords_grid_1.size() + coords_grid_2.size());

    auto it = std::set_union(coords_grid_1.begin(),coords_grid_1.end(),
                             coords_grid_2.begin(),coords_grid_2.end(),
                             coords_grid_union.begin());

    coords_grid_union.resize(it-coords_grid_union.begin());
  }
  auto grid_union = Grid<dim>::create(knots_union);
  //---------------------------------------------------------


  //---------------------------------------------------------
  map_elem_id_grid_union_to_elem_id_grid_1 =
    build_map_elements_id_between_grids(*grid_union,grid_1);
  map_elem_id_grid_union_to_elem_id_grid_2 =
    build_map_elements_id_between_grids(*grid_union,grid_2);
  //---------------------------------------------------------

  return grid_union;
}

}

IGA_NAMESPACE_CLOSE

#include <igatools/geometry/grid_tools.inst>

