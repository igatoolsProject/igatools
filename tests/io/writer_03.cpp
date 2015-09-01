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
 * Testing the writer, this the cell data
 * author: pauletti
 * date: Jun 21, 2014
 *
 */

#include "../tests.h"
#include "igatools/io/writer.h"

template<int dim>
void
test()
{
  const int n_knots = 4;
  auto grid = CartesianGrid<dim>::create(n_knots);
  Writer<dim> writer(grid);

  SafeSTLVector<Real> cell_data(grid->get_num_all_elems());
  int n=1;
  for (auto elem : *grid)
  {
    cell_data[elem.get_flat_index()] = n;
    n *= -1;
  }
  writer.add_element_data(cell_data, "chess board");

  string filename = "grid" + to_string(dim);
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

