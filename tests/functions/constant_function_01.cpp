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
#include <igatools/geometry/grid_function_lib.h>
#include <igatools/geometry/domain.h>

#include <igatools/base/quadrature_lib.h>

#include "../tests.h"
#include "function_test.h"

using namespace functions;

template<int dim, int codim, int range, int rank>
void constant_func()
{
  using Function = functions::ConstantFunction<dim, codim, range, rank>;

  auto grid = Grid<dim>::const_create(3);
  auto ball_func = grid_functions::BallGridFunction<dim>::const_create(grid);
  auto ball_domain = Domain<dim,0>::const_create(ball_func);

  typename Function::Value b;
  for (int i=0; i<range; ++i)
    b[i] = i;

  auto func = Function::const_create(ball_domain, b);
  function_values(*func);
}

int main()
{
  constant_func<1,0,1,1>();
  constant_func<2,0,2,1>();

  return 0;
}

