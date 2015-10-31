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

#include "../tests.h"

#include <igatools/base/quadrature_lib.h>
#include <igatools/basis_functions/nurbs_space.h>
#include <igatools/basis_functions/nurbs_element.h>

/*
 *  Test for the NURBS space iterator using the no cache evaluations
 *
 *  author: pauletti
 *  date: Jun 11, 2014
 *
 */

template< int dim, int range, int rank = 1>
void test()
{
  OUTSTART

  const int r = 2;
//    out << "test<" << dim << "," << range << ">" << endl;

  using Space = NURBSSpace< dim, range, rank >;
  auto  knots = Grid<dim>::const_create();

  auto degree = TensorIndex<dim>(r);

  auto bsp_space = BSpline<dim,range,rank>::const_create(
                     SplineSpace<dim,range,rank>::const_create(degree,knots));

  using ScalarSpSpace = BSpline<dim,1,1>;
  auto scalar_bsp_space = ScalarSpSpace::const_create(
                            SplineSpace<dim,1,1>::const_create(degree,knots));


  const auto n_scalar_basis = scalar_bsp_space->get_num_basis();

  using WeightFunc = IgGridFunction<dim,1>;

  IgCoefficients weights;
  for (int dof = 0 ; dof < n_scalar_basis ; ++dof)
    weights[dof] = 1.0;

  const auto w_func = WeightFunc::const_create(scalar_bsp_space,weights);

  auto space = Space::const_create(bsp_space,w_func);

  const int n_points = 3;
  auto quad = QGauss<dim>::create(n_points);

  auto elem     = space->begin();
  auto end = space->end();

  using _Value = space_element::_Value;
  using _Gradient = space_element::_Gradient;
  using _Hessian = space_element::_Hessian;

  int elem_id = 0;
  for (; elem != end; ++elem)
  {
    out << "Element: " << elem_id<< endl;

    out << "Values basis functions:" << endl ;
    const auto values = elem->template evaluate_basis_at_points<_Value>(quad,DofProperties::active);
    values.print_info(out) ;

    out << "Gradients basis functions:" << endl ;
    const auto gradients = elem->template evaluate_basis_at_points<_Gradient>(quad,DofProperties::active);
    gradients.print_info(out) ;

    out << "Hessians basis functions:" << endl ;
    const auto hessians = elem->template evaluate_basis_at_points<_Hessian>(quad,DofProperties::active);
    hessians.print_info(out) ;

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
