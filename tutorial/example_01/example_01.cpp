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

// [logstream]
LogStream out;
// [logstream]

int main()
{
  // [dim]
  const int dim = 2;
  // [dim]

  // [grid]
  const int n_knots = 3;
  shared_ptr<const Grid<dim>> grid = Grid<dim>::const_create(n_knots);
  // [grid]

  // [grid_print]
  out << "Grid: " << endl;
  out << "  number of elements:  " << endl;
  out << "    " << grid->get_num_all_elems();
  out << " = " << grid->get_num_intervals() << endl;
  out << "  knots:  " << endl;
  auto knots = grid->get_knots();
  for (const auto &knot_vect : knots)
  {
    out << "   ";
    for (const auto &knot : *knot_vect)
      out << " " << knot;
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



