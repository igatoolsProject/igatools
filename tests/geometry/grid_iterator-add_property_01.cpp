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
 *  @brief  elem add_property()
 *  @author pauletti
 *  @date  Aug 19, 2015
 */

#include "../tests.h"

#include <igatools/geometry/cartesian_grid.h>
#include <igatools/geometry/grid_element.h>

template <int dim>
void iterate(const int n_knots = 5)
{
  OUTSTART

  const string red = "red";
  const string blue = "blue";

  auto grid = Grid<dim>::create(n_knots);
  grid->add_property(red);
  grid->add_property(blue);

  const TensorIndex<dim> center(2);
  for (auto &elem : *grid)
  {
    if (elem.get_index() <= center)
      elem.add_property(red);
    else
      elem.add_property(blue);
  }

  auto elem_r = grid->cbegin(red);
  auto end_r  = grid->cend(red);
  for (; elem_r != end_r; ++elem_r)
  {
    elem_r->print_info(out);
  }

  auto elem_b = grid->cbegin(blue);
  auto end_b  = grid->cend(blue);
  for (; elem_b != end_b; ++elem_b)
  {
    elem_b->print_info(out);
  }

  OUTEND
}



int main()
{
  iterate<0>();
  iterate<1>();
  iterate<2>();
  iterate<3>();

  return  0;
}
