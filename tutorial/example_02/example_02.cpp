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
 * Example for looping on the elements of a grid-like container
 * and accessing information that do not require the use of cache
 */
// [old_includes]
#include <igatools/geometry/cartesian_grid.h>
#include <igatools/basis_functions/bspline_space.h>
// [old_includes]
// [acc_includes]
#include <igatools/geometry/grid_element.h <igatools/basis_functions/bspline_element.h>
// [acc_includes]
#include <igatools/base/logstream.h>

using namespace iga;
using namespace std;

LogStream out;

// [templated_function]
template <int dim>
void loop_on_grid()
{
// [templated_function]
  // [create_grid]
  out << "Traversing the elements of a " << dim;
  out << "-dimensional grid." << endl;
  const int n_knots = 3;
  auto grid = CartesianGrid<dim>::create(n_knots);
  // [create_grid]
  // [iter_grid]
  for (auto elem : *grid)
  {
    out << "The tensor index of element: " << elem.get_flat_index();
    out << " is: "<< elem.get_tensor_index() << endl;
  }
  // [iter_grid]
  out << endl;
}


template <int dim>
void loop_on_space()
{
  out << "Traversing the elements of a " << dim;
  out << "-dimensional B-spline space." << endl;
  const int n_knots = 3;
  auto grid = CartesianGrid<dim>::create(n_knots);
  const int degree = 2;
  auto space = BSplineSpace<dim>::create(degree, grid);

  for (auto elem : *space)
  {
    out << "Element: " << elem.get_flat_index();
    out << " has global basis: ";
    elem.get_local_to_global(DofProperties::active).print_info(out);
    out << endl;
  }
  out << endl;
}


int main()
{
  loop_on_grid<1>();
  loop_on_grid<2>();
  loop_on_grid<3>();

  out << endl;

  loop_on_space<1>();
  loop_on_space<2>();
  loop_on_space<3>();

  return 0;
}



