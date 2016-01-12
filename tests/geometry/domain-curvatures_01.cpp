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
 *  Test for the curvature of SphericalGridFunction
 *
 *  author: pauletti
 *  date: 2014-10-24
 *
 */

#include "../tests.h"

#include <igatools/functions/grid_function_lib.h>
#include <igatools/base/quadrature_lib.h>
#include <igatools/geometry/domain.h>
#include <igatools/geometry/domain_element.h>
#include <igatools/geometry/domain_handler.h>


template <int dim>
void principal_curvatures()
{
  OUTSTART

  out.begin_item("principal_curvatures<" + std::to_string(dim) +">");

  using Sphere = grid_functions::SphereGridFunction<dim>;

  using Flags = domain_element::Flags;
  auto flag = Flags::ext_normal |  Flags::curvature;


  BBox<dim> box;
  for (int i=0; i<dim-1; ++i)
    box[i] = {0.+M_PI/8, M_PI-M_PI/8};
  if (dim>=1)
    box[dim-1] = {0., M_PI};
  auto grid = Grid<dim>::const_create(box, 2);

  auto F = Sphere::const_create(grid);

  auto domain = Domain<dim,1>::const_create(F);

  auto domain_handler = domain->create_cache_handler();

  domain_handler->set_element_flags(flag);

  auto elem = domain->begin();
  auto end = domain->end();

  auto quad = QUniform<dim>::create(3);
  domain_handler->init_cache(elem,quad);
  int elem_id = 0;
  for (; elem != end; ++elem, ++elem_id)
  {
    out.begin_item("Element " + std::to_string(elem_id));

    out << "Element ID: " << elem->get_index() << std::endl;

    domain_handler->fill_element_cache(elem);

    out.begin_item("Normals:");
    elem->get_exterior_normals().print_info(out);
    out.end_item();

    out.begin_item("Curvature:");
    elem->get_principal_curvatures().print_info(out);
    out.end_item();

    out.end_item();
  }

  out.end_item();
  OUTEND
}


int main()
{
  out.depth_console(10);

  principal_curvatures<1>();
  principal_curvatures<2>();

  return 0;
}
