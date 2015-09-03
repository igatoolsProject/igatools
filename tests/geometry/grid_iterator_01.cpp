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
 *  @brief  Simultaneous use of element with two sub elements
 *  @author pauletti
 *  @date   2015-08-19
 */

#include "../tests.h"

#include <igatools/base/quadrature_lib.h>
#include <igatools/geometry/grid.h>
#include <igatools/geometry/grid_cache_handler.h>
#include <igatools/geometry/grid_element.h>


template <int dim, int sdim = dim-1>
void iterate(const int n_knots = 5)
{
  OUTSTART

  using Grid = Grid<dim>;
  using Flags = typename Grid::ElementAccessor::Flags;
  auto grid = Grid::const_create(n_knots);

  auto flag = Flags::weight;
  auto s_flag = Flags::point;
  auto cache_handler = grid->create_cache_handler();
  cache_handler->template set_flags<dim>(flag);
  cache_handler->template set_flags<sdim>(s_flag);

  auto quad   = QGauss<dim>::create(2);
  auto s_quad = QGauss<sdim>::create(1);

  auto elem = grid->cbegin();
  cache_handler->template init_cache<dim>(elem, quad);
  cache_handler->template init_cache<sdim>(elem, s_quad);

  for (; elem != grid->cend(); ++elem)
  {
    cache_handler->template fill_cache<dim>(elem, 0);
    elem->template get_weights<dim>(0).print_info(out);
    out << endl;

    for (auto &s_id : UnitElement<dim>::template elems_ids<sdim>())
    {
      cache_handler->template fill_cache<sdim>(elem, s_id);
      elem->template get_points<sdim>(s_id).print_info(out);
      out << endl;
    }
    out << endl;
  }

  OUTEND
}



int main()
{
  iterate<1>();
  iterate<2>();
  iterate<3>();

  return  0;
}
