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
 *  @brief Grid get_boundary_normals()
 *  @author pauletti
 *  @date  2015-08-19
 */

#include "../tests.h"
#include <igatools/geometry/grid.h>

template <int dim, int sdim>
void boundary_normals()
{
  OUTSTART

  auto grid = Grid<dim>::create();
  const int n_elems = UnitElement<dim>::template num_elem<sdim>();
  for (int j = 0; j < n_elems; ++j)
  {
    out << "Sub element index: " << j << endl;
    auto normals = grid->template get_boundary_normals<sdim>(j);
    out << "Outer normals: ";
    for (int i=0; i<dim-sdim; ++i)
      out << normals[i] << " " << endl;
  }

  OUTEND
}



int main()
{
  boundary_normals<1,0>();
  boundary_normals<2,1>();
  boundary_normals<3,2>();

  return 0;
}
