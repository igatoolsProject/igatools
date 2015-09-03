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

/**
 *  @file
 *  @brief  Grid::get_sub_grid
 *  @author pauletti
 *  @date  2015-08-19
 */

#include "../tests.h"
#include <igatools/geometry/grid.h>
#include <igatools/geometry/grid_element.h>

template<int dim, int sdim = dim-1>
void get_subgrid(const TensorSize<dim> &n_knots)
{
  OUTSTART

  using Grid =  Grid<dim>;
  auto grid = Grid::const_create(n_knots);
  out.begin_item("Grid:");
  grid->print_info(out);
  out.end_item();

  for (auto &i : UnitElement<dim>::template elems_ids<sdim>())
  {
    typename Grid<dim>::template SubGridMap<sdim> map;
    out.begin_item("Sub element: " + to_string(i));
    auto sub_grid = grid->template get_sub_grid<sdim>(i, map);
    sub_grid->print_info(out);
    map.print_info(out);
    out.end_item();
  }

  OUTEND
}



int main()
{
  get_subgrid<1>(TensorSize<1>(sequence<1>(2)));
  get_subgrid<2>(TensorSize<2>(sequence<2>(2)));
  get_subgrid<3>(TensorSize<3>(sequence<3>(2)));

  return  0;
}
