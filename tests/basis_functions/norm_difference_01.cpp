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
 *  Test for the norm difference function.
 *
 *  author: pauletti
 *  date: 2015-03-17
 */

#include "../tests.h"

#include "common_functions.h"
#include <igatools/base/quadrature_lib.h>
#include <igatools/basis_functions/space_tools.h>
#include <igatools/geometry/grid_function_lib.h>


template<int dim, int range = 1>
void norm_difference(const int deg, const int n_knots = 10)
{
  out.begin_item("norm_difference<"
                 + to_string(dim) + ">("
                 + to_string(deg) + ","
                 + to_string(n_knots) + ")");
  const Real p=2.;


  auto grid = Grid<dim>::const_create(n_knots);

  const int n_qpoints = ceil((2*dim + 1)/2.);
  auto quad = QGauss<dim>::const_create(n_qpoints);

  auto f = ProductGridFunction<dim,range>::const_create(grid);
  auto g = grid_functions::ConstantGridFunction<dim,range>::const_create(grid, {0.});

  SafeSTLMap<ElementIndex<dim>,Real> elem_err;
  auto err = space_tools::l2_norm_difference<dim,range>(*f,*g,quad,elem_err);
  out.begin_item("Error L2:");
  out << "Expected: " << std::pow(p+1, -dim/p) << endl;
  out << "Computed: " << err << endl;
  out.end_item();

  SafeSTLMap<ElementIndex<dim>,Real> elem_err_h1;
  auto err_h1 = space_tools::h1_norm_difference<dim,range>(*f,*g,quad,elem_err_h1);
  out.begin_item("Error H1:");
  out <<  err_h1 << endl;
  out.end_item();

  SafeSTLMap<ElementIndex<dim>,Real> elem_err_inf;
  auto err_inf = space_tools::inf_norm_difference<dim,range>(*f,*g,quad,elem_err_inf);
  out.begin_item("Error inf:");
  out <<  err_inf << endl;
  out.end_item();

  out.end_item();

}

template<int dim>
void norm_difference_p(const int deg, const int n_knots, const Real p)
{
  out.begin_item("norm_difference_p<"
                 + to_string(dim) + ">("
                 + to_string(deg) + ","
                 + to_string(n_knots) + ","
                 + to_string(p) + ")");
  auto grid = Grid<dim>::const_create(n_knots);

  const int n_qpoints = ceil((2*dim + 1)/2.);
  auto quad = QGauss<dim>::const_create(n_qpoints);

  auto f = ProductGridFunction<dim,1>::const_create(grid);
  auto g = grid_functions::ConstantGridFunction<dim,1>::const_create(grid, {0.});

  SafeSTLMap<ElementIndex<dim>,Real> elem_err;
  space_tools::norm_difference_grid_functions<0,dim>(*f, *g, quad, p, elem_err);

  Real err = 0.;
  for (const auto &loc_err : elem_err)
    err += loc_err.second;

  err = std::pow(err,1./p);

  out.begin_item("Error L2:");
  out << "Expected: " << std::pow(p+1, -dim/p) << endl;
  out << "Computed: " << err << endl;
  out.end_item();

  out.end_item();
}



int main()
{
  out.depth_console(20);

  norm_difference<1>(3);
  norm_difference<2>(3);
  norm_difference<3>(1);

  norm_difference_p<1>(3, 10, 1);
  norm_difference_p<2>(3, 10, 1);
  norm_difference_p<3>(3, 10, 1);
  return 0;
}

