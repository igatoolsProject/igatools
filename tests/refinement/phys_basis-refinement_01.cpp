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
 *  Test refinement of a basic PhysicalBasis using the BSpline as reference basis
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
#include <igatools/basis_functions/physical_basis.h>
#include <igatools/basis_functions/physical_basis_element.h>
#include <igatools/functions/grid_function_lib.h>





template <int dim>
void test_evaluate()
{
  auto grid = Grid<dim>::create();
  grid->refine();
  out << endl;

  const int deg = 2;

  using RefBasis = ReferenceBasis<dim>;
  using RefBasisPtr = std::shared_ptr<RefBasis>;
  RefBasisPtr ref_basis = BSpline<dim>::create(SplineSpace<dim>::create(deg,grid));
  auto phys_basis =
    PhysicalBasis<dim,1,1,0>::create(
      ref_basis,
      Domain<dim,0>::create(grid_functions::IdentityGridFunction<dim>::create(grid)));


  out << endl;
  out << endl;

  out << "===============================================================" << endl;
  out.begin_item("O R I G I N A L     B A S I S");
  phys_basis->print_info(out);
  out.end_item();
  out << "===============================================================" << endl;
  out << endl;

  out << "===============================================================" << endl;
  out.begin_item("R E F I N E D     B A S I S");
  phys_basis->refine_h();
  phys_basis->print_info(out);
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
