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
 *  Test for the NURBS class subbasis extraction
 *
 *  author: martinelli
 *  date: 19 Oct 2015
 */

#include "../tests.h"
//#include <igatools/base/quadrature_lib.h>
#include <igatools/basis_functions/nurbs.h>


template<int sub_dim, int dim, int range=1, int rank=1>
void sub_basis(TensorSize<dim> n, const int degree = 1)
{
  OUTSTART

  auto grid = Grid<dim>::create(n);

  auto deg = TensorIndex<dim>(degree);
  auto bsp_basis = BSpline<dim,range,rank>::const_create(
                     SplineSpace<dim,range,rank>::create(degree,grid));

  auto scalar_bsp_basis = BSpline<dim,1,1>::create(
                            SplineSpace<dim,1,1>::create(degree,grid));
  const auto n_scalar_basis = scalar_bsp_basis->get_num_basis();

  using WeightFunc = IgGridFunction<dim,1>;

  IgCoefficients weights;
  for (int dof = 0 ; dof < n_scalar_basis ; ++dof)
    weights[dof] = dof * 1.0;

  const auto w_func = WeightFunc::const_create(scalar_bsp_basis,weights);

  using Basis = NURBS<dim,range,rank>;
  auto nrb_basis = Basis::const_create(bsp_basis,w_func);


  typename Basis::template InterBasisMap<sub_dim> dof_map;

  out.begin_item("Dim: " + std::to_string(dim) + "     Sub-Dim: " + std::to_string(sub_dim));


  out.begin_item("Original NURBS:");
  nrb_basis->print_info(out);
  out.end_item();

  for (auto i : UnitElement<dim>::template elems_ids<sub_dim>())
  {
    typename Grid<dim>::template SubGridMap<sub_dim> elem_map;
    out.begin_item(to_string(i) + "-th " + "sub basis:");
    auto sub_basis =
      nrb_basis->template get_sub_basis<sub_dim>(i, dof_map, elem_map);
    out.begin_item("Basis:");
    sub_basis->print_info(out);
    out.end_item();
  }
  out.end_item();

  OUTEND
}



int main()
{

  sub_basis<0,1>(TensorSize<1>(sequence<1>(2)));
  sub_basis<1,2>(TensorSize<2>(sequence<2>(2)));
  sub_basis<2,3>(TensorSize<3>(sequence<3>(2)));

  return  0;
}

