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

/*
 *  Test cylindrical annulus mapping values.
 *  author: antolin
 *  date: 2014-03-18
 *
 */

#include "../tests.h"

#include <igatools/base/quadrature_lib.h>
#include <igatools/geometry/grid_function_lib.h>
#include <igatools/geometry/grid_function_element.h>
#include <igatools/geometry/domain.h>
#include <igatools/geometry/domain_element.h>
#include <igatools/geometry/domain_handler.h>


std::shared_ptr<const GridFunction<3,3>>
                                      create_cylinder_function()
{
  auto grid = Grid<3>::create();

  using Cyl = grid_functions::CylindricalAnnulusGridFunction;
  auto cylinder = Cyl::create(grid,1, 2, 0, 2.0, 0.0, numbers::PI/2.0);

  return cylinder;
}

template <int dim>
void domain_values(const std::shared_ptr<const GridFunction<dim,dim>> &func)
{
  auto domain = Domain<dim,0>::const_create(func);


  auto domain_handler = domain->create_cache_handler();

  using Flags = domain_element::Flags;
  auto flag = Flags::point |
              Flags::jacobian |
              Flags::hessian |
              Flags::measure |
              Flags::w_measure;
  domain_handler->set_element_flags(flag);

  auto elem = domain->begin();
  auto end = domain->end();

  auto quad = QGauss<dim>::create(1);
  domain_handler->init_element_cache(elem,quad);
  int elem_id = 0;
  for (; elem != end; ++elem, ++elem_id)
  {
    out.begin_item("Element " + std::to_string(elem_id));

    out << "Element ID: " << elem->get_index() <<std::endl;

    domain_handler->fill_element_cache(elem);

    out.begin_item("Ref. Points:");
    elem->get_grid_function_element().get_grid_element().get_element_points().print_info(out);
    out.end_item();

    out.begin_item("Points:");
    elem->get_element_points().print_info(out);
    out.end_item();

    out.begin_item("Jacobians:");
    elem->get_element_jacobians().print_info(out);
    out.end_item();

    out.begin_item("Hessians:");
    elem->get_element_hessians().print_info(out);
    out.end_item();

    out.begin_item("Measure:");
    elem->get_element_measures().print_info(out);
    out.end_item();

    out.begin_item("Weight * Measure:");
    elem->get_element_w_measures().print_info(out);
    out.end_item();

    out.end_item();
  }
}


int main()
{
  out.depth_console(10);

  auto cyl = create_cylinder_function();
  domain_values<3>(cyl);

  return 0;
}


