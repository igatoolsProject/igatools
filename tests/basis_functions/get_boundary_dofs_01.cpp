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
 *  Test for space_tools::get_boundary_dofs
 *
 *  author: pauletti
 *  date: 2015-03-27
 *
 */

#include "../tests.h"

#include <igatools/basis_functions/bspline_space.h>
#include <igatools/basis_functions/space_tools.h>

using space_tools::get_boundary_dofs;

template<int dim, int range = 1, int rank = 1>
void get_bdry_dof(const int deg = 1, const int n_knots = 3)
{
  OUTSTART
  using RefSpace = ReferenceSpace<dim, range, rank>;
  using Space = BSplineSpace<dim, range, rank>;
  auto grid = CartesianGrid<dim>::create(n_knots);
  grid->set_boundary_id(0, 1);

  auto space = Space::create(deg, grid);

  std::set<boundary_id>  piece_one  = {1};
  std::set<boundary_id>  piece_zero = {0};
  auto one_dofs = get_boundary_dofs<RefSpace>(space, piece_one);
  auto zero_dofs = get_boundary_dofs<RefSpace>(space, piece_zero);

  // TODO (pauletti, Mar 27, 2015): we should create iga::set with print_info
  for (auto &x : one_dofs)
    out << x << " ";
  out << endl;

  for (auto &x : zero_dofs)
    out << x << " ";
  out << endl;

  OUTEND
}



int main()
{
  get_bdry_dof<1>();
  get_bdry_dof<2>();
  get_bdry_dof<3>();

  get_bdry_dof<1>(2);
  get_bdry_dof<2>(2);

  return 0;
}
