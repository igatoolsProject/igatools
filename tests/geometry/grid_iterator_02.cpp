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
 *  @brief  One cache handler with two different elements
 *  @author pauletti
 *  @date   2015-08-19
 */

#include "../tests.h"

#include <igatools/base/quadrature_lib.h>
#include <igatools/geometry/grid.h>
#include <igatools/geometry/grid_handler.h>
#include <igatools/geometry/grid_element.h>


template <int dim>
void handler_two_elems(const int n_knots = 3)
{
  OUTSTART

  using Grid = Grid<dim>;
  using Flags = typename Grid::ElementAccessor::Flags;
  auto grid = Grid::const_create(n_knots);

  auto flag = Flags::point;

  auto cache_handler = grid->create_cache_handler();
  cache_handler->template set_flags<dim>(flag);

  auto quad1 = QGauss<dim>::create(2);
  auto quad2 = QGauss<dim>::create(1);

  auto elem1 = grid->cbegin();
  auto elem2 = grid->cbegin();

  cache_handler->template init_cache<dim>(elem1, quad1);
  cache_handler->template init_cache<dim>(elem2, quad2);

  cache_handler->template fill_cache<dim>(elem1, 0);
  elem1->template get_points<dim>(0).print_info(out);
  out << endl;

  ++elem1;
  cache_handler->template fill_cache<dim>(elem1, 0);
  elem1->template get_points<dim>(0).print_info(out);
  out << endl;

  cache_handler->template fill_cache<dim>(elem2, 0);
  elem2->template get_points<dim>(0).print_info(out);
  out << endl;

  OUTEND
}


int main()
{

  handler_two_elems<1>();
  handler_two_elems<2>();
  handler_two_elems<3>();

  return  0;
}
