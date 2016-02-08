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
 *  Test for the integrate_difference function.
 *
 *  author: pauletti
 *  date: 26 Jun 2014
 */

#include "../tests.h"

#include "common_functions.h"
#include <igatools/base/quadrature_lib.h>
#include <igatools/basis_functions/space_tools.h>
#include <igatools/functions/grid_function_lib.h>

using std::to_string;

template<int dim>
void integrate_grid_function(const int deg,  const int n_knots)
{
  out.begin_item("integrate_grid_function<"
                 + to_string(dim) + ">("
                 + to_string(deg) + ","
                 + to_string(n_knots) +")");

  auto grid = Grid<dim>::const_create(n_knots);
  const int n_qpoints = ceil((2*dim + 1)/2.);
  auto quad = QGauss<dim>::const_create(n_qpoints);


  auto f = NormGridFunction<dim>::const_create(grid);
  auto g = grid_functions::ConstantGridFunction<dim,1>::const_create(grid, {0.});

  SafeSTLMap<ElementIndex<dim>,Real> elem_err_l2;
  Real err_l2 = space_tools::l2_norm_difference<dim,1>(*f, *g, quad, elem_err_l2);
  out << "Error L2 = "<< err_l2 << endl;

  SafeSTLMap<ElementIndex<dim>,Real> elem_err_h1;
  Real err_h1 = space_tools::h1_norm_difference<dim,1>(*f, *g, quad, elem_err_h1);
  out << "Error H1 = "<< err_h1 << endl;

  out.end_item();
}



int main()
{
  out.depth_console(20);

  integrate_grid_function<1>(1,3);
  integrate_grid_function<2>(1,3);
  integrate_grid_function<3>(1,3);


  return 0;
}

