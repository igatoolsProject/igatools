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
 *  @brief Sphere Domain
 *  @author pauletti
 *  @date 2015-09-02
 */

#include <igatools/base/quadrature_lib.h>
#include <igatools/functions/grid_function_lib.h>
#include <igatools/geometry/domain.h>
#include <igatools/geometry/domain_element.h>

#include "../tests.h"

template <int dim>
void sphere_domain()
{
  OUTSTART

  const int codim = 1;
  using Grid = Grid<dim>;
  using Domain = Domain<dim, codim>;
  using GridFunc = grid_functions::SphereGridFunction<dim>;

  BBox<dim> bbox;
  for (int i=0; i<dim; ++i)
    bbox[i]= {M_PI/6, M_PI/2};
  auto grid = Grid::const_create(bbox, 2);
  auto grid_func = GridFunc::const_create(grid);
  auto domain = Domain::const_create(grid_func);

  using Flags = typename Domain::ElementAccessor::Flags;

  auto flag = Flags::measure | Flags::w_measure | Flags::point |
              Flags::ext_normal;

  auto handler = domain->create_cache_handler();
  handler->template set_flags<dim>(flag);
  auto quad = QUniform<dim>::create(3);
  auto elem = domain->cbegin();
  handler->init_cache(elem, quad);
  auto  end = domain->cend();

  for (; elem != end; ++elem)
  {
    handler->template fill_cache<dim>(elem, 0);

    out << "Points:" << endl;
    elem->template get_points<dim>(0).print_info(out);
    out << endl;
    out << "Measure:" << endl;
    elem->template get_measures<dim>(0).print_info(out);
    out << endl;
    out << "weight * measure:" << endl;
    elem->template get_w_measures<dim>(0).print_info(out);
    out << endl;
    out << "Exterior normal:" << endl;
    elem->template get_exterior_normals().print_info(out);
    out << endl;
  }

  OUTEND
}


int main()
{
  sphere_domain<1>();
  sphere_domain<2>();

  return 0;
}
