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
 *  Test for the NURBS basis iterator
 *
 *  author: pauletti
 *  date: Jun 11, 2014
 *
 *  author: martinelli
 *  date: Dec 02, 2014
 *
 *
 */

#include "../tests.h"

#include <igatools/base/quadrature_lib.h>
#include <igatools/basis_functions/nurbs.h>
#include <igatools/basis_functions/nurbs_element.h>


template< int dim, int range, int rank = 1>
void test()
{
  OUTSTART
  const int r = 2;

  using Basis = NURBS< dim, range, rank >;
  auto  knots = Grid<dim>::const_create(3);

  auto degree = TensorIndex<dim>(r);
  auto bsp_basis = BSpline<dim,range,rank>::const_create(
                     SplineSpace<dim,range,rank>::const_create(degree,knots));

  using ScalarBspBasis = BSpline<dim,1,1>;
  auto scalar_bsp_basis = ScalarBspBasis::const_create(
                            SplineSpace<dim,1,1>::const_create(degree,knots));

  const auto n_scalar_basis = scalar_bsp_basis->get_num_basis();


  IgCoefficients weights;
  for (int dof = 0 ; dof < n_scalar_basis ; ++dof)
    weights[dof] = (dof + 1) * (1.0 / n_scalar_basis) ;

  using WeightFunc = IgGridFunction<dim,1>;
  const auto w_func = WeightFunc::const_create(scalar_bsp_basis,weights);

  auto basis = Basis::const_create(bsp_basis,w_func);

  const int n_points = 3;
  auto quad = QGauss<dim>::create(n_points);

  auto elem     = basis->begin();
  auto end_element = basis->end();

  using Flags = basis_element::Flags;
  const auto flag = Flags::value|
                    Flags::gradient|
                    Flags::hessian;



  auto elem_handler = basis->create_cache_handler();
  elem_handler->template set_flags<dim>(flag);

  using _Value = basis_element::_Value;
  using _Gradient = basis_element::_Gradient;
  using _Hessian = basis_element::_Hessian;

  elem_handler->template init_cache<dim>(*elem,quad);

  int elem_id = 0;
  for (; elem != end_element; ++elem)
  {
    elem_handler->fill_element_cache(elem);
    std::cout << "Element flat id: " << elem_id << endl << endl;
    out << "Element flat id: " << elem_id << endl << endl;

    out.begin_item("Values basis functions:");
    auto values = elem->template get_basis_data<_Value,dim>(0,DofProperties::active);
    values.print_info(out);
    out.end_item();

    out.begin_item("Gradients basis functions:");
    auto gradients = elem->template get_basis_data<_Gradient,dim>(0,DofProperties::active);
    gradients.print_info(out);
    out.end_item();

    out.begin_item("Hessians basis functions:");
    auto hessians = elem->template get_basis_data<_Hessian,dim>(0,DofProperties::active);
    hessians.print_info(out);
    out.end_item();

    ++elem_id;
  }
  OUTEND
}


int main()
{
  test<1, 1>();
  test<1, 2>();
  test<1, 3>();
  test<2, 1>();
  test<2, 2>();
  test<2, 3>();
  test<3, 1>();
  test<3, 3>();

  return 0;
}
