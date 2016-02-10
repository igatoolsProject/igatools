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
 *  Test for the l2_projection function.
 *  Bspline bases, constant function, SafeSTLVector valued case
 *
 *  author: pauletti
 *  date: 2013-10-31
 */

#include "../tests.h"


#include <igatools/base/quadrature_lib.h>
#include <igatools/functions/grid_function_lib.h>

#include <igatools/basis_functions/bspline.h>
#include <igatools/basis_functions/bspline_element.h>

#include <igatools/basis_functions/basis_tools.h>

template<int dim , int range=1>
void test_proj(const int deg, const int n_knots = 4)
{
  OUTSTART


  auto grid = Grid<dim>::const_create(n_knots);
  auto space = SplineSpace<dim,range>::const_create(deg, grid);
  auto basis = BSpline<dim,range>::const_create(space);


  using Func = typename grid_functions::ConstantGridFunction<dim,range>;
  typename Func::Value val;
  for (int i=0; i<range; ++i)
    val[i] = i+3;

  auto f = Func::const_create(grid,val);

  const int n_qp = 4;
  auto quad = QGauss<dim>::create(n_qp);

  auto coeffs_func = basis_tools::projection_l2_grid_function<dim,range>(*f,*basis,quad);

  auto proj_func = IgGridFunction<dim,range>::const_create(basis,coeffs_func);
  proj_func->print_info(out);

  OUTEND

}



int main()
{

//  test_proj<0,1>(3);
  test_proj<1,1>(3);
  test_proj<2,1>(3);
  test_proj<3,1>(1);

  test_proj<2,3>(1);

  return 0;
}

