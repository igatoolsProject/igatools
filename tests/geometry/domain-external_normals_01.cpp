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

/*
 *  Test for the SphericalFunction class as a mapping
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

#include <igatools/linear_algebra/dense_matrix.h>
#include "Teuchos_LAPACK.hpp"

template <int dim>
void ext_normal_values()
{
  OUTSTART

  out.begin_item("ext_normal_values<" + std::to_string(dim) + ">");

  using Sph = grid_functions::SphereGridFunction<dim>;



  BBox<dim> box;
  for (int i=0; i<dim-1; ++i)
    box[i] = {0.+M_PI/8, M_PI-M_PI/8};
  if (dim>=1)
    box[dim-1] = {0., M_PI};
  auto grid = Grid<dim>::const_create(box, 2);

  auto sphere = Sph::const_create(grid);


  auto domain = Domain<dim,1>::const_create(sphere);

  auto domain_handler = domain->create_cache_handler();

  using Flags = domain_element::Flags;
  auto flag = Flags::point | Flags::jacobian | Flags::ext_normal;
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

    out.begin_item("Ref. Points:");
    elem->get_grid_function_element().get_grid_element().get_element_points().print_info(out);
    out.end_item();

    out.begin_item("Points:");
    elem->get_element_points().print_info(out);
    out.end_item();

    out.begin_item("Jacobians:");
    elem->get_element_jacobians().print_info(out);
    out.end_item();

    out.end_item();
  }

  out.end_item();

  OUTEND
}


int main()
{
  ext_normal_values<1>();
  ext_normal_values<2>();

  return 0;
}
