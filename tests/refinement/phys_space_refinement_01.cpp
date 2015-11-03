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
 *  Test refinement of a basic PhysicalSpaceBasis using the BSpline as reference space
 *  and the IdentityFunction as mapping.
 *
 *  author: pauletti
 *  date: 2013-10-02
 *
 */

#include "../tests.h"

#include <igatools/base/quadrature_lib.h>
#include <igatools/geometry/push_forward.h>
#include <igatools/basis_functions/bspline.h>
#include <igatools/basis_functions/physical_space_basis.h>
#include <igatools/basis_functions/physical_space_element.h>
#include <igatools/geometry/grid_function_lib.h>





template <int dim>
void test_evaluate()
{
  auto grid = Grid<dim>::create();
  grid->refine();
  out << endl;

  const int deg = 2;

  using RefSpace = ReferenceSpace<dim>;
  using RefSpacePtr = std::shared_ptr<RefSpace>;
  RefSpacePtr ref_space = BSpline<dim>::create(SplineSpace<dim>::create(deg,grid));
  auto phys_space =
    PhysicalSpaceBasis<dim,1,1,0>::create(
      ref_space,
      Domain<dim,0>::create(grid_functions::IdentityGridFunction<dim>::create(grid)));


  out << endl;
  out << endl;

  out << "===============================================================" << endl;
  out.begin_item("O R I G I N A L     S P A C E");
  phys_space->print_info(out);
  out.end_item();
  out << "===============================================================" << endl;
  out << endl;

  out << "===============================================================" << endl;
  out.begin_item("R E F I N E D     S P A C E");
  phys_space->refine_h();
  phys_space->print_info(out);
  out.end_item();
  out << "===============================================================" << endl;
  out << endl;
}

int main()
{
  out.depth_console(10);

  test_evaluate<1>();
  out << endl ;

  test_evaluate<2>();
  out << endl ;

  test_evaluate<3>();
  out << endl ;

  return 0;
}
