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
 *  Test of for the faces normals using a cylindrical annulus mapping.
 *  author: antolin
 *  date: 2014-03-18
 *
 */

#include "../tests.h"

#include <igatools/base/quadrature_lib.h>
#include <igatools/functions/grid_function_lib.h>
#include <igatools/geometry/domain.h>
#include <igatools/geometry/domain_element.h>
#include <igatools/geometry/domain_handler.h>


auto create_mapping1(const shared_ptr<const Grid<3>> &grid)
{
  using Cyl = grid_functions::CylindricalAnnulusGridFunction;

  return Cyl::const_create(grid);
}

template <int dim>
auto create_mapping2(const shared_ptr<const Grid<dim>> &grid)
{
  using Identity = grid_functions::IdentityGridFunction<dim>;

  return Identity::const_create(grid);
}

template <int sub_dim>
void boundary_normals()
{
  auto grid = Grid<3>::const_create({{1.,2.},{0.,numbers::PI/2.0},{0.,2.}});
  auto map_func = create_mapping1(grid);
  auto domain = Domain<3,0>::const_create(map_func);

  auto domain_handler = domain->create_cache_handler();

  using Flags = domain_element::Flags;
  auto flag = Flags::w_measure| Flags::point| Flags::boundary_normal;

  auto quad = QGauss<sub_dim>::create(1);

  domain_handler->template set_flags<sub_dim>(flag);

  auto elem = domain->begin();
  auto end = domain->end();

  domain_handler->init_cache(elem,quad);

  int elem_id = 0;
  for (; elem != end; ++elem, ++elem_id)
  {
    out.begin_item("Element " + std::to_string(elem_id));

    out << "Element ID: " << elem->get_index() << std::endl;

    for (auto &s_id : UnitElement<3>::template elems_ids<sub_dim>())
    {
      out.begin_item("Face: " + std::to_string(s_id));

      domain_handler->template fill_cache<sub_dim>(elem, s_id);


      out.begin_item("Normal vector:");
      elem->template get_boundary_normals<sub_dim>(s_id).print_info(out);
      out.end_item();

      out.end_item();
    }

    out.end_item();
  }



}

int main()
{
  boundary_normals<2>();
  return 0;
}

