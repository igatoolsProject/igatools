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
 *  Test for the BSpline element iterator local to global
 *  author: pauletti
 *  date: Aug 21, 2014
 *
 */

#include "../tests.h"

#include <igatools/base/quadrature_lib.h>
#include <igatools/basis_functions/bspline_element.h>
#include <igatools/basis_functions/bspline_element_handler.h>



template <int dim, int range=1, int rank=1>
void elem_dofs(const int n_knots = 4, const int deg=1)
{
  OUTSTART

  auto grid  = Grid<dim>::const_create(n_knots);
  auto space = SplineSpace<dim, range, rank>::const_create(deg, grid);
  auto basis = BSpline<dim, range, rank>::const_create(space);

  auto elem = basis->begin();
  auto end = basis->end();

  for (; elem != end; ++elem)
  {
    out << "Element index: " << elem->get_index() << endl;
    out << "Global dofs: ";
    elem->get_local_to_global(DofProperties::active).print_info(out);
    out << endl;
  }

  OUTEND
}


int main()
{
  out.depth_console(10);

  elem_dofs< 1, 1, 1 >();
  elem_dofs< 1, 2, 1 >();
  elem_dofs< 1, 3, 1 >();
  elem_dofs< 2, 1, 1 >();
  elem_dofs< 2, 2, 1 >();
  elem_dofs< 2, 3, 1 >();
  elem_dofs< 3, 1, 1 >();
  elem_dofs< 3, 3, 1 >();

  for (int p=0; p<3; p++)
  {
    elem_dofs< 1, 1, 1 >(4,p);
    elem_dofs< 1, 2, 1 >(4,p);
    elem_dofs< 1, 3, 1 >(4,p);
  }

  return  0;
}
