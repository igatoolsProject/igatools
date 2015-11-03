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
#include <igatools/basis_functions/bspline.h>
#include <igatools/basis_functions/space_tools.h>


template<int dim, int range = 1, int rank = 1>
void norm_difference(const int deg, const int n_knots = 10)
{
  const Real p=2.;
  using Basis = BSpline<dim, range, rank>;


  auto grid = Grid<dim>::const_create(n_knots);
  auto space = Basis::const_create(deg, grid);

  const int n_qpoints = ceil((2*dim + 1)/2.);
  QGauss<dim> quad(n_qpoints);

  auto f = std::shared_ptr<ProductFunction<dim> >(new ProductFunction<dim>(grid, IdentityFunction<dim>::const_create(grid)));
  typename functions::ConstantFunction<dim,0,1>::Value val {0.};
  auto g = functions::ConstantFunction<dim,0,1>::const_create(grid,
                                                              IdentityFunction<dim>::const_create(grid), val);

  SafeSTLVector<Real> elem_err(grid->get_num_all_elems());
  auto err = space_tools::l2_norm_difference(*f, *g, quad, elem_err);
  out << std::pow(p+1, -dim/p) << "\t" << err << endl;

  SafeSTLVector<Real> elem_err_h1(grid->get_num_all_elems());
  auto err_h1 = space_tools::h1_norm_difference(*f, *g, quad, elem_err_h1);
  out <<  err_h1 << endl;

  SafeSTLVector<Real> elem_err_inf(grid->get_num_all_elems());
  auto err_inf = space_tools::inf_norm_difference(*f, *g, quad, elem_err_inf);
  out <<  err_inf << endl;

}

template<int dim, int range = 1, int rank = 1>
void norm_difference_p(const int deg, const int n_knots, const Real p)
{
  using Basis = BSpline<dim, range, rank>;


  auto grid = Grid<dim>::const_create(n_knots);
  auto space = Basis::const_create(deg, grid);

  const int n_qpoints = ceil((2*dim + 1)/2.);
  QGauss<dim> quad(n_qpoints);

  auto f = std::shared_ptr<ProductFunction<dim> >(new ProductFunction<dim>(grid, IdentityFunction<dim>::const_create(grid)));
  typename functions::ConstantFunction<dim,0,1>::Value val {0.};
  auto g = functions::ConstantFunction<dim,0,1>::const_create(grid,
                                                              IdentityFunction<dim>::const_create(grid), val);

  SafeSTLVector<Real> elem_err(grid->get_num_all_elems());
  space_tools::norm_difference<0,dim>(*f, *g, quad, p, elem_err);

  Real err = 0;
  for (const Real &loc_err : elem_err)
    err += loc_err;

  err = std::pow(err,1./p);

  out << std::pow(p+1, -dim/p) << "\t" << err << endl;

}



int main()
{
  out.depth_console(20);

  norm_difference<1,1,1>(3);
  norm_difference<2,1,1>(3);
  norm_difference<3,1,1>(1);

  norm_difference_p<1,1,1>(3, 10, 1);
  norm_difference_p<2,1,1>(3, 10, 1);
  norm_difference_p<3,1,1>(3, 10, 1);
  return 0;
}

