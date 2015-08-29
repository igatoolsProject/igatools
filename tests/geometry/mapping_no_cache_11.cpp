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
 *  Test for linear mapping class, without the use of the cache.
 *  It is based on mapping_11.cpp
 *  author: martinelli
 *  date: Feb 03, 2015
 */

#include "../tests.h"

#include <igatools/geometry/mapping.h>
#include <igatools/geometry/mapping_element.h>
#include <igatools/functions/function_lib.h>
#include <igatools/functions/identity_function.h>
#include <igatools/base/quadrature_lib.h>
#include <igatools/functions/function_element.h>

template<int dim, int codim>
void test()
{
  OUTSTART

  const int space_dim = dim+codim;
  using Function = functions::LinearFunction<dim, 0, space_dim>;
  using Mapping  = Mapping<dim, codim>;

  typename Function::Value    b;
  typename Function::Gradient A;
  for (int i=0; i<space_dim; i++)
  {
    for (int j=0; j<dim; j++)
      if (j == i)
        A[j][j] = 2.;
    b[i] = i;
  }
  auto grid = CartesianGrid<dim>::create(3);
  auto F = Function::create(grid, IdentityFunction<dim>::create(grid), A, b);

  auto quad = QGauss<dim>(2);

  Mapping map(F);

  auto elem = map.begin();
  auto end  = map.end();

  for (; elem != end; ++elem)
  {
    out << "Values:" << endl;
    elem->template evaluate_at_points<_Value>(quad).print_info(out);
    out << endl;
    out << "Gradients:" << endl;
    elem->template evaluate_at_points<_Gradient>(quad).print_info(out);
    out << endl;
    out << "Hessians:" << endl;
    elem->template evaluate_at_points<_Hessian>(quad).print_info(out);
    out << endl;
  }

  OUTEND
}


int main()
{
  test<2,0>();
  test<3,0>();
  test<2,1>();

  return 0;
}

