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
 *  Test for LinearFunction defined over an Ig Domain
 *  author: antolin
 *  date: Jan 07, 2016
 */

#include "../tests.h"

#include <igatools/base/quadrature_lib.h>
#include <igatools/functions/ig_grid_function.h>
#include <igatools/functions/function_lib.h>
#include <igatools/functions/function_element.h>
#include <igatools/basis_functions/bspline.h>

#include "function_test.h"


template<int dim, int codim>
shared_ptr<const Domain<dim, codim>> create_ig_domain()
{
  static const int range = dim + codim;
  static const int rank = 1;
  auto grid = Grid<dim>::const_create(2);
  const int deg = 1;
  auto space = SplineSpace<dim, range,rank>::const_create(deg, grid);
  auto basis = BSpline<dim,range,rank>::const_create(space);

  using Function = IgGridFunction<dim, range>;

  Epetra_SerialComm comm;
  auto map = EpetraTools::create_map(*basis, "active", comm);
  auto coeff = EpetraTools::create_vector(*map);
  (*coeff)[0] = 1.;

  auto func = Function::const_create(basis, *coeff);

  return Domain<dim, codim>::const_create(func);
}



template<int dim, int codim, int range>
void test_linear_function()
{
  using std::to_string;
  out.begin_item("test_linear_function<dim=" + to_string(dim) +
                 ",codim=" + to_string(codim)+ ",range=" + to_string(range));

  using LinearFunc = functions::LinearFunction<dim, codim, range>;

  typename LinearFunc::Value    b;
  typename LinearFunc::Gradient A;
  for (int i=0; i<range; i++)
  {
    for (int j=0; j<dim; j++)
      if (j == i)
        A[j][j] = 2.;
    b[i] = i;
  }

  const auto domain = create_ig_domain<dim, codim>();

  auto F = LinearFunc::const_create(domain, A, b);
  function_values(*F);

  out.end_item();
}




int main()
{
  test_linear_function<1, 0, 1>();
//  test_linear_function<2, 0, 2>();
//  test_linear_function<3, 0, 3>();

  return 0;
}

