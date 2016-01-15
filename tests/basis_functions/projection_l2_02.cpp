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
 *  Test for the boundary projection function.
 *  Physical bases version
 *  author: pauletti
 *  date: 2013-10-10
 *
 */

#include "../tests.h"
#include "./common_functions.h"

#include <igatools/base/quadrature_lib.h>
//#include <igatools/functions/formula_function.h>

#include <igatools/functions/grid_function_lib.h>

#include <igatools/basis_functions/bspline.h>
#include <igatools/basis_functions/bspline_element.h>

#include <igatools/basis_functions/space_tools.h>

#include <igatools/io/writer.h>
#include <igatools/basis_functions/physical_basis.h>
#include <igatools/basis_functions/physical_basis_element.h>
#include <igatools/basis_functions/physical_basis_handler.h>

using numbers::PI;




template<int dim, int codim, int range>
void do_test(const int p, const int num_knots = 10)
{
  auto grid = Grid<dim>::const_create(num_knots);
  auto space = SplineSpace<dim,range>::const_create(p, grid);
  auto ref_basis = BSpline<dim,range>::const_create(space);

  using F = grid_functions::LinearGridFunction<dim,dim + codim>;
  typename F::Value    b;
  typename F::Gradient A;
  for (int i = 0; i < dim; ++i)
  {
    A[i][i] = 1+i;
  }
  auto map = F::const_create(grid,A, b);
  auto domain = Domain<dim,codim>::const_create(map);

  auto basis = PhysicalBasis<dim,range,1,codim>::const_create(ref_basis,domain);

  const int n_qpoints = 4;
  auto quad = QGauss<dim>::const_create(n_qpoints);

  auto f = TestFunction<dim,range>::const_create(domain);
  auto coeffs_func = space_tools::projection_l2_function<dim,codim,range,1>(*f, *basis, quad);

  auto proj_func = IgFunction<dim,codim,range,1>::const_create(basis,coeffs_func);
  proj_func->print_info(out);

}



int main()
{
  out.depth_console(20);
  // do_test<1,1,1>(3);
  do_test<2,0,1>(3);
  //do_test<3,1,1>(1);

  return 0;
}

