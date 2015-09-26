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
 *  Test for the bspline element iterator
 *  Computes values and derivatives of the basis functions without the use of
 *  the cache.
 *  This test computes the same quantities of bspline_element_iterator_03.cpp
 *
 *  author: martinelli
 *  date: May 08, 2013
 *
 */
//TODO (pauletti, Apr 3, 2015): rename to bspline_element_iterato_??
#include "../tests.h"

#include <igatools/base/quadrature_lib.h>
#include <igatools/basis_functions/bspline_space.h>
#include <igatools/basis_functions/bspline_element.h>
//#include <igatools/linear_algebra/distributed_vector.h>
#include <igatools/linear_algebra/epetra_vector.h>
#include <igatools/linear_algebra/dof_tools.h>

using namespace EpetraTools;

template< int dim, int range >
void do_test()
{
  auto grid = Grid<dim>::const_create();

  const int degree = 1;
  const int rank =  1 ;
  using Space = BSplineSpace< dim, range, rank >;
  auto space = Space::const_create(degree, grid);

  const auto  num = space->get_num_basis();
  Epetra_SerialComm comm;
  Epetra_Map map(num, num, 0, comm);

  Vector u(map);
  {
    int id = 0 ;
    u[id++] = 0.0 ;
    u[id++] = 1.0 ;

    u[id++] = 0.0 ;
    u[id++] = 1.0 ;

    u[id++] = 0.0 ;
    u[id++] = 0.0 ;

    u[id++] = 1.0 ;
    u[id++] = 1.0 ;
  }

  QGauss< dim > quad_tensor_prod(2) ;
  const auto eval_points = quad_tensor_prod.get_points();

  auto quad_non_tensor_prod = Quadrature<dim>::create(eval_points);

  auto element1 = space->begin();

  using Elem = typename Space::ElementAccessor;
  using _Value = typename Elem::_Value;
  using _Gradient = typename Elem::_Gradient;
  using _Hessian = typename Elem::_Hessian;
  using _Divergence = typename Elem::_Divergence;

  out.begin_item("Basis values using QGauss<" + std::to_string(dim) + "> with 2 points in each coordinate direction.");
  auto values1 = element1->template evaluate_basis_at_points<_Value>(quad_non_tensor_prod,DofProperties::active);
  values1.print_info(out);
  out.end_item();

  out.begin_item("Basis gradients using QGauss<" + std::to_string(dim) + "> with 2 points in each coordinate direction.");
  auto gradients1 = element1->template evaluate_basis_at_points<_Gradient>(quad_non_tensor_prod,DofProperties::active);
  gradients1.print_info(out);
  out.end_item();

//    QUniform< dim > quad_scheme_2(3) ;
  auto quad_scheme_2 = Quadrature< dim >::create(QUniform< dim >(3).get_points()) ;

  out.begin_item("Basis values using QUniform<" + std::to_string(dim) + "> with 2 points in each coordinate direction.");
  auto values2 = element1->template evaluate_basis_at_points<_Value>(quad_scheme_2,DofProperties::active);
  values2.print_info(out);
  out.end_item();

  out.begin_item("Basis gradients using QUniform<" + std::to_string(dim) + "> with 2 points in each coordinate direction.");
  auto gradients2 = element1->template evaluate_basis_at_points<_Gradient>(quad_scheme_2,DofProperties::active);
  gradients2.print_info(out);
  out.end_item();



}


int main(int argc, char *argv[])
{
  out.depth_console(10);

//   do_test<1,1>();
//    do_test<1,2>();
//    do_test<1,3>();
//
//    do_test<2,1>();
  do_test<2,2>();
//    do_test<2,3>();
//
//    do_test<3,1>();
  do_test<3,3>();

  return 0;
}
