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

// [includes]
#include <igatools/basis_functions/bspline.h>
#include <igatools/base/logstream.h>
// [includes]

// [using]
using namespace iga;
using namespace std;
// [using]

int main()
{
  // [dim]
  const int dim = 2;
  // [dim]

  // [grid]
  const int n_knots = 3;

  // the type of the variable grid is shared_ptr<const Grid<dim>>
  auto grid = Grid<dim>::const_create(n_knots);
  // [grid]

  // [logstream]
  LogStream out;
  // [logstream]

  // [grid_print]
  out << "Grid: " << endl;
  out << "  number of elements:  " << endl;
  //
  // total number of elements in the grid
  out << "    " << grid->get_num_all_elems();
  //
  // number of intervals in each coordinate direction
  out << " = " << grid->get_num_intervals() << endl;
  //
  out << "  knots:  " << endl;
  //
  // the type for knots is SafeSTLArray<shared_ptr<SafeSTLVector<Real>>,dim>
  auto knots = grid->get_knots();
  //
  // SafeSTLArray can be used ad std::array --> the range-based for loops works fine
  for (const auto &knot_vect : knots)
  {
    out << "    ";
    knot_vect->print_info(out);
    out << endl;
  }
  out << endl;
  // [grid_print]

  // [space]
  const int degree = 2;
  auto space = SplineSpace<dim>::const_create(degree, grid);
  // [space]

  // [basis]
  auto basis = BSpline<dim>::const_create(space);
  out.begin_item("BSpline:");
  basis->print_info(out);
  out.end_item();
  // [basis]

  return 0;
}



