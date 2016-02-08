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
 *  Test for a domain built using a LinearGridFunction
 *
 *  author: pauletti
 *  date: Oct 11, 2014
 */

#if 0
#include "../tests.h"

#include <igatools/geometry/domain_element.h>
#include <igatools/base/quadrature_lib.h>
#include <igatools/functions/grid_function_lib.h>
#endif

#include "domain_values.h"

template<int dim, int codim>
void test()
{
  using std::to_string;
  out.begin_item("test<dim=" + to_string(dim) + ",codim=" + to_string(codim));

  const int space_dim = dim + codim;
  using F = grid_functions::LinearGridFunction<dim,space_dim>;

  typename F::Value    b;
  typename F::Gradient A;
  for (int i=0; i<space_dim; i++)
  {
    for (int j=0; j<dim; j++)
      if (j == i)
        A[j][j] = 2.;
    b[i] = i;
  }

  auto quad = QGauss<dim>::const_create(2);
  auto grid = Grid<dim>::create(3);
  auto grid_func = F::create(grid,A,b);
  auto domain = Domain<dim,codim>::create(grid_func);


  domain_values<dim,codim>(*domain,quad);

  out.end_item();
}


int main()
{
  test<2,0>();
//    test<3,3>();

  return 0;
}

