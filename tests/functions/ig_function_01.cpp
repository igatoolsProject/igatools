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
 *  Test for Function class, as a prototype for an spline function
 *  author: pauletti
 *  date: Oct 11, 2014
 */

#include "../tests.h"

#include <igatools/functions/ig_function.h>
#include <igatools/functions/grid_function_lib.h>
#include <igatools/base/quadrature_lib.h>
#include <igatools/functions/function_element.h>
#include <igatools/basis_functions/bspline.h>
#include <igatools/basis_functions/bspline_element.h>


#include "function_test.h"

template<int dim, int range>
void test()
{
  auto grid = Grid<dim>::create(3);
  const int deg = 2;
  auto space = SplineSpace<dim,range>::create(deg, grid);
  auto ref_basis = BSpline<dim,range>::create(space);

  auto grid_func = grid_functions::IdentityGridFunction<dim>::create(grid);
  auto domain = Domain<dim>::create(grid_func);

  auto phys_basis = PhysicalBasis<dim,range>::create(ref_basis,domain);


  Epetra_SerialComm comm;
  auto map = EpetraTools::create_map(*phys_basis, "active", comm);
  auto coeff = EpetraTools::create_vector(*map);
  (*coeff)[0] = 1.;

  auto func = IgFunction<dim,0,range,1>::create(phys_basis, *coeff);

  function_values(*func);
}


int main()
{
  test<2,1>();
//    test<3,3>();

  return 0;
}

