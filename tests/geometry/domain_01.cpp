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
 *  @brief Developing domain class
 *  @author pauletti
 *  @date 2015
 */
#include <igatools/base/quadrature_lib.h>
#include <igatools/geometry/domain.h>
#include <igatools/geometry/domain_element.h>
#include <igatools/geometry/domain_handler.h>
#include <igatools/geometry/grid_function_lib.h>
#include "../tests.h"


template<int dim, int codim>
void domain()
{
  OUTSTART

  using Grid = Grid<dim>;
  //using Domain = Domain<dim, codim>;

  using Func = grid_functions::BallGridFunction<dim>;
  auto grid = Grid::const_create();
  auto func = Func::const_create(grid);
  auto domain = Domain<dim,codim>::const_create(func);

  using Flags = typename domain_element::Flags;

  auto handler = domain->create_cache_handler();

  auto flag = Flags::w_measure | Flags::point;
  handler->template set_flags<dim>(flag);

  auto quad   = QGauss<dim>::create(2);

  auto elem = domain->cbegin();
  auto end = domain->cend();
  handler->init_cache(elem, quad);

  for (; elem != end; ++elem)
  {
    handler->template fill_cache<dim>(elem,0);
    elem->template get_points<dim>(0).print_info(out);
//    elem->get_element_values_D1().print_info(out);
    elem->template get_w_measures<dim>(0).print_info(out);
    out << endl;
  }

  OUTEND
}


int main()
{
  //domain<0,0>();
  domain<1,0>();
  domain<2,0>();
  domain<3,0>();
  //domain<2,1>();

  return 0;
}

