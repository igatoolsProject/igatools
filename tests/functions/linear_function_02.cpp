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
 *  Test for LinearFunction
 *  author: pauletti
 *  date: Oct 12, 2014
 */

#include "../tests.h"

#include <igatools/functions/identity_function.h>
#include <igatools/base/quadrature_lib.h>
#include <igatools/geometry/grid_function_lib.h>
#include <igatools/functions/function_lib.h>
#include <igatools/functions/function_element.h>

#include "function_test.h"





template<int dim, int codim, int range>
void test_linear_function()
{
  using std::to_string;
  out.begin_item("test_linear_function<dim=" + to_string(dim) +
                 ",codim=" + to_string(codim)+ ",range=" + to_string(range));

  using Function = functions::LinearFunction<dim, codim, range>;

  typename Function::Value    b;
  typename Function::Gradient A;
  for (int i=0; i<range; i++)
  {
    for (int j=0; j<dim; j++)
      if (j == i)
        A[j][j] = 2.;
    b[i] = i;
  }


  auto grid = Grid<dim>::create(3);
  auto grid_func = grid_functions::IdentityGridFunction<dim>::create(grid);
  auto domain = Domain<dim>::create(grid_func);
  auto F = Function::create(domain, A, b);
  function_values(*F);

  out.end_item();
}




int main()
{
  test_linear_function<1, 0, 1>();
  test_linear_function<2, 0, 2>();
  test_linear_function<3, 0, 3>();

  return 0;
}

