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
 *  @brief  get_measure
 *  @author pauletti
 *  @date   2015-08-19
 */

#include "../tests.h"

#include <igatools/geometry/grid.h>
#include <igatools/geometry/grid_element.h>

#include "../tests_utils.h"

template <int dim, int sdim = dim>
void iterate()
{
  OUTSTART

  auto grid = non_uniform_grid<dim>();

  for (auto &elem : *grid)
  {
    for (auto &s_id : UnitElement<dim>::template elems_ids<sdim>())
      out << elem.template get_measure<sdim>(s_id) << endl;
    out << endl;
  }

  OUTEND
}



int main()
{
  iterate<0>();
  iterate<1>();
  iterate<2>();
  iterate<3>();

  iterate<1,0>();
  iterate<2,1>();
  iterate<3,2>();

  return  0;
}
