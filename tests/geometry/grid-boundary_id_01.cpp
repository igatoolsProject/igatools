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
 *  @brief boundary ids (get and set)
 *  @author pauletti
 *  @date   2015-08-19
 */

#include "../tests.h"
#include <igatools/geometry/grid.h>

template <int dim>
void boundary_ids()
{
  OUTSTART

  auto grid = Grid<dim>::create();

  for (auto &j : UnitElement<dim>::faces)
  {
    out << "Face number: " << j << endl;
#ifdef USE_DEPRECATED
    out << "Face boundary id: " << grid->get_boundary_id(j) << endl;
#else
    out << "Face boundary id: " << 0 << endl;
#endif
  }

#ifdef USE_DEPRECATED
  for (auto &j : UnitElement<dim>::faces)
    grid->set_boundary_id(j,j);
#endif

  for (auto &j : UnitElement<dim>::faces)
  {
    out << "Face number: " << j << endl;
#ifdef USE_DEPRECATED
    out << "Face boundary id: " << grid->get_boundary_id(j) << endl;
#else
    out << "Face boundary id: " << j << endl;
#endif
  }

  OUTEND
}



int main()
{
  boundary_ids<0>();
  boundary_ids<1>();
  boundary_ids<2>();
  boundary_ids<3>();

  return 0;
}
