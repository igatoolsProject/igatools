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
 * Testing the writer, this test the BallGridFunction
 * author: pauletti
 * date: Jun 21, 2014
 *
 */

#include "../tests.h"
#include "igatools/io/writer.h"
#include "igatools/geometry/grid_function_lib.h"

template<int dim>
void
test()
{
  const int n_knots = 4;
  auto grid = Grid<dim>::create(n_knots);
  auto ball_function = grid_functions::BallGridFunction<dim>::create(grid);
  Writer<dim> writer(ball_function, 4);

  string filename = "ball_dim" + to_string(dim);
  writer.save(filename, true);
  writer.print_info(out);
}


int main()
{
  test<1>();
  test<2>();
  test<3>();

  return 0;
}

