//-+--------------------------------------------------------------------
// Igatools a general purpose Isogeometric analysis library.
// Copyright (C) 2012-2016  by the igatools authors (see authors.txt).
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
 *  @brief  elem.is_boundary()
 *  @author pauletti
 *  @date   2015-08-19
 */

#include "../tests.h"

#include <igatools/base/quadrature_lib.h>
#include <igatools/geometry/grid.h>
#include <igatools/geometry/grid_handler.h>
#include <igatools/geometry/grid_element.h>

template <int dim, int sdim>
void boundaryId(const TensorSize<dim> &n_knots)
{
  OUTSTART

  using Grid = Grid<dim>;
  auto grid = Grid::const_create(n_knots);

  for (const auto &elem : *grid)
  {
    elem.print_info(out);
    out.begin_item("Boundary subelements:");
    for (auto &s_id : UnitElement<dim>::template elems_ids<sdim>())
    {
      if (elem.template is_boundary<sdim>(s_id))
        out << s_id << ",";
    }
    out.end_item();
  }

  OUTEND
}


int main()
{
  boundaryId<0,0>(TensorSize<0>(3));
  boundaryId<1,1>(TensorSize<1>(3));
  boundaryId<2,2>(TensorSize<2>(3));
  boundaryId<3,3>(TensorSize<3>(3));

  boundaryId<1,0>(TensorSize<1>(3));
  boundaryId<2,1>(TensorSize<2>(3));
  boundaryId<3,2>(TensorSize<3>(3));

  boundaryId<1,0>(TensorSize<1>(sequence<1>(2)));
  boundaryId<2,1>(TensorSize<2>(sequence<2>(2)));
  boundaryId<3,2>(TensorSize<3>(sequence<3>(2)));

  return  0;
}
