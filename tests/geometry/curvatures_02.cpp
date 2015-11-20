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
 *  Test for the SphereGridFunction class as a mapping
 *
 *  author: pauletti
 *  date: 2014-10-24
 *
 */

#include "../tests.h"

#include <igatools/geometry/grid_function_lib.h>
#include <igatools/base/quadrature_lib.h>
#include <igatools/geometry/domain.h>
#include <igatools/geometry/domain_element.h>
#include <igatools/geometry/domain_handler.h>


template <int dim>
void normal_derivatives()
{
  OUTSTART

  out.begin_item("normal_derivatives<" + std::to_string(dim) + ">");

  BBox<dim> box;
  for (int i=0; i<dim-1; ++i)
    box[i] = {0.+M_PI/8, M_PI-M_PI/8};
  if (dim>=1)
    box[dim-1] = {0., M_PI};
  auto grid = Grid<dim>::create(box, 2);

  using Sphere = grid_functions::SphereGridFunction<dim>;

  auto sphere_func = Sphere::create(grid);



  auto sphere_domain = Domain<dim,1>::create(sphere_func);

  auto domain_handler = sphere_domain->create_cache_handler();

  using Flags = domain_element::Flags;
  auto flag = Flags::ext_normal |  Flags::ext_normal_D1;
//  auto flag = Flags::ext_normal_D1;
  domain_handler->set_element_flags(flag);

  auto elem = sphere_domain->begin();
  auto end = sphere_domain->end();

  auto quad = QUniform<dim>::create(3);

  domain_handler->init_cache(elem,quad);
  for (; elem != end; ++elem)
  {
    domain_handler->fill_element_cache(elem);


    out.begin_item("Normals:");
    auto normals = elem->get_exterior_normals();
    normals.print_info(out);
    out.end_item();
    //*/

    out.begin_item("Der normal:");
    auto D_normals = elem->get_exterior_normals_D1();
    D_normals.print_info(out);
    out.end_item();

    out.begin_item("Dn^t on n:");
    for (int pt=0; pt<normals.get_num_points(); ++pt)
//            out << action(transpose(D_normals[pt]), normals[pt]) << endl;
      out << action(co_tensor(transpose(D_normals[pt])), normals[pt][0]) << endl;
    out.end_item();

    //*/
  }

  out.end_item();
  OUTEND
}


int main()
{
  out.depth_console(10);

  normal_derivatives<1>();
  normal_derivatives<2>();

  return 0;
}
