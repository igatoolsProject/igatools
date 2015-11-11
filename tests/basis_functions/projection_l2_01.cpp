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
 *  Bspline spaces case
 *
 *  author: pauletti
 *  date: 2013-10-10
 */

#include "../tests.h"
#include "./common_functions.h"

#include <igatools/base/quadrature_lib.h>

#include <igatools/basis_functions/bspline.h>
#include <igatools/basis_functions/bspline_element.h>

#include <igatools/basis_functions/space_tools.h>






template<int dim , int range>
void project_l2(const int p, const int num_knots = 10)
{
  OUTSTART

  auto knots = Grid<dim>::const_create(num_knots);
  auto space = SplineSpace<dim,range>::const_create(p, knots) ;
  auto basis = BSpline<dim,range>::const_create(space) ;


  const int n_qpoints = 4;
  auto quad = QGauss<dim>::const_create(n_qpoints);

  auto f = TestGridFunction<dim,range>::const_create(knots);
  auto coeffs_func = space_tools::projection_l2_grid_function<dim,range>(*f,*basis,quad);

  auto proj_func = IgGridFunction<dim,range>::const_create(basis,coeffs_func);
  proj_func->print_info(out);

  OUTEND
}



int main()
{
//  project_l2<0,1>(1);
  project_l2<1,1>(3);
  project_l2<2,1>(3);
  project_l2<3,1>(1);

  return 0;
}

