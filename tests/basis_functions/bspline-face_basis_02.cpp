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
 *  Test for the BSpline class subbasis extraction
 *
 *  author: pauletti
 *  date: 2014-11-18
 */

#include "../tests.h"
#include <igatools/base/quadrature_lib.h>
#include <igatools/basis_functions/bspline.h>


template<int sub_dim, int dim, int range=1, int rank=1>
void sub_basis(TensorSize<dim> n, const int degree = 1)
{
  OUTSTART


  using Space = SplineSpace<dim, range, rank>;
  using Basis = BSpline<dim, range, rank>;

  auto grid = Grid<dim>::const_create(n);
  auto space = Space::const_create(degree, grid);
  auto basis = Basis::const_create(space);

  typename Basis::template InterBasisMap<sub_dim> dof_map;

  for (auto i : UnitElement<dim>::template elems_ids<sub_dim>())
  {
    typename Grid<dim>::template SubGridMap<sub_dim> elem_map;
//    auto sub_grid = space->get_grid()->template get_sub_grid<sub_dim>(i, elem_map);
    out.begin_item(to_string(i) + "-th " + "sub basis:");
    auto sub_basis =
      basis->template get_sub_basis<sub_dim>(i, dof_map, elem_map);
    out.begin_item("Basis:");
    sub_basis->print_info(out);
    out.end_item();

    out.begin_item("Dofs sub element to basis mapping:");
    dof_map.print_info(out);
    out.end_item();
    out.end_item();
  }

  OUTEND
}



int main()
{

  sub_basis<0,1>(TensorSize<1>(sequence<1>(2)));
  sub_basis<1,2>(TensorSize<2>(sequence<2>(2)));
  sub_basis<2,3>(TensorSize<3>(sequence<3>(2)));

  return  0;
}

