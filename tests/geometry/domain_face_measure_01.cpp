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
 *  Test for the Domain build with CylindricalAnnulusGridFunction class.
 *
 *  author: antolin
 *  date: 2014-03-19
 *
 */

#include "../tests.h"

#include <igatools/base/quadrature_lib.h>
#include <igatools/functions/grid_function_lib.h>
#include <igatools/geometry/domain.h>
#include <igatools/geometry/domain_element.h>
#include <igatools/geometry/domain_handler.h>

// TODO (pauletti, Nov 20, 2014): this test is more about cylindrical annulus
void test_evaluate()
{
  const int dim = 3;
  const int sdim = 2;

  auto grid = Grid<dim>::const_create({{0.0,numbers::PI/2.0},{1.,2.},{0.0,1.0}});

  using Cyl = grid_functions::CylindricalAnnulusGridFunction;
  auto cyl = Cyl::const_create(grid);

  auto domain = Domain<3,0>::const_create(cyl);

  auto domain_handler = domain->create_cache_handler();

  using Flags = domain_element::Flags;
  auto flag = Flags::w_measure;
  domain_handler->template set_flags<sdim>(flag);


  auto elem = domain->begin();
  auto end  = domain->end();


  SafeSTLArray<Real, UnitElement<dim>::template num_elem<sdim>()> face_area ;
  std::fill(face_area.begin(), face_area.end(), 0.0) ;

  auto quad = QGauss<sdim>::create(3);

  domain_handler->init_cache(elem,quad);

  for (; elem != end; ++elem)
  {
    const auto &grid_elem = elem->get_grid_function_element().get_grid_element();
    if (grid_elem.is_boundary())
    {
      for (auto &s_id : UnitElement<dim>::template elems_ids<sdim>())
      {
        if (grid_elem.is_boundary(s_id))
        {
          domain_handler->template fill_cache<sdim>(elem, s_id);
          auto w_meas = elem->template get_w_measures<sdim>(s_id);
          for (auto &w : w_meas)
            face_area[s_id] += w;
        }
      }
    }
  }

  out << "Dimension " << dim << endl;
  for (Index face_id = 0; face_id < UnitElement<dim>::template num_elem<sdim>(); ++face_id)
  {
    out << "Area of face " << face_id << " : " << face_area[face_id] << endl;
  }


}

int main()
{
  out.depth_console(10);

  test_evaluate();

}
