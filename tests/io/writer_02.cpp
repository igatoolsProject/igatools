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
 * Testing the writer, this test the mapping
 * author: pauletti
 * date: Jun 21, 2014
 *
 */

#include "../tests.h"
#include "igatools/io/writer.h"
#include "igatools/functions/identity_function.h"
#include "igatools/functions/function_lib.h"
//#include "igatools/geometry/mapping_element.h"

template<int dim>
void
test()
{
  const int n_knots = 4;
  auto grid = CartesianGrid<dim>::create(n_knots);
  auto identity_function = IdentityFunction<dim>::create(grid);
  auto ball_function = functions::BallFunction<dim>::create(grid,identity_function);
  Writer<dim> writer(ball_function, 4);

  string filename = "map" + to_string(dim);
  writer.save(filename);
  writer.save(filename,"appended");
  writer.print_info(out);
}


int main()
{
  test<1>();
  test<2>();
  test<3>();

  return 0;
}

