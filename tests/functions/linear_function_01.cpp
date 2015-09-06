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

/**
 *  @file
 *  @brief ConstantFunction
 *  @author pauletti
 *  @date 2013-11-01
 */
#include <igatools/geometry/domain_lib.h>
#include <igatools/functions/function_element.h>
#include <igatools/functions/function_lib.h>

#include <igatools/base/quadrature_lib.h>

#include "../tests.h"
#include "function_test.h"

using namespace functions;

template<int dim, int codim, int range>
void linear_func()
{
  using Grid = Grid<dim>;
  using Domain = domains::BallDomain<dim>;
  using Function = functions::LinearFunction<dim, codim, range>;

  auto grid = Grid::const_create(3);
  auto domain = Domain::const_create(grid);

  typename Function::template Derivative<1> A;
  typename Function::Value b;


  for (int j=0; j<range; ++j)
  {
    b[j] = j;
    for (int i=0; i<dim; ++i)
      A[i][j] = i+j;
  }
  auto func = Function::const_create(domain, A, b);
  function_values<dim, codim, range>(func);
}

int main()
{
  linear_func<1,0,1>();
  linear_func<2,0,2>();

  return 0;
}

