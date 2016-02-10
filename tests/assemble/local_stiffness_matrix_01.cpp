//-+--------------------------------------------------------------------
// Igatools a general purpose Isogeometric analysis library.
// Copyright (C) 2012-2016  by the igatools authors (see authors.txt).
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
 *  Test the assembling of local stiffness matrix using a BSpline basis
 *
 *  author: pauletti
 *  date: Oct 09, 2014
 *
 */

#include "../tests.h"

#include <igatools/basis_functions/physical_basis_handler.h>
#include <igatools/basis_functions/bspline.h>
#include <igatools/basis_functions/bspline_handler.h>
#include <igatools/basis_functions/bspline_element.h>
#include <igatools/base/quadrature_lib.h>


//#define TIME_PROFILING

using std::to_string;

template<int dim,int range>
void loc_stiff_matrix(const int n_knots, const int deg)
{

  OUTSTART

  out.begin_item("local_stiff_matrix<" + to_string(dim) + "," + to_string(range) + ">");

  TensorIndex<dim> deg_dir;
  TensorSize<dim> num_pts;
  for (int i = 0 ; i < dim ; ++i)
  {
    deg_dir[i] = deg+i;
    num_pts[i] = deg_dir[i] + 1;
  }


  auto grid = Grid<dim>::create(n_knots);
  using Basis = BSpline<dim,range>;
  auto basis = Basis::create(SplineSpace<dim,range>::create(deg_dir, grid)) ;

  auto elem_handler = basis->create_cache_handler();


  using Flags = basis_element::Flags;
  auto flag = Flags::gradient | Flags::w_measure;
  elem_handler->template set_flags<dim>(flag);

  auto elem           = basis->begin();
  const auto elem_end = basis->end();

  auto quad = QGauss<dim>::create(num_pts);
  elem_handler->init_element_cache(elem,quad);

  int elem_id = 0;
  for (; elem != elem_end; ++elem, ++elem_id)
  {
    out.begin_item("Element " + std::to_string(elem_id));

    out << "Element ID: " << elem->get_index() << std::endl;

    elem_handler->fill_element_cache(elem);

    auto loc_mat = elem->template integrate_gradu_gradv<dim>(0);
//    auto loc_mat = elem->integrate_gradu_gradv_sum_factorization_impl(Topology<dim>(),0);

#ifndef TIME_PROFILING
    out.begin_item("Stiffnes matrix:");
    loc_mat.print_info(out);
    out.end_item();
#endif

    out.end_item();
  }

  out.end_item();

  OUTEND

}




int main()
{
#ifdef TIME_PROFILING

  const int n_knots = 2;
  const int deg = 20;

//    loc_stiff_matrix<1>(n_knots, deg);
//    loc_stiff_matrix<2>(n_knots, deg);
  loc_stiff_matrix<3>(n_knots, deg);

#else

  const int n_knots = 3;
  const int deg = 3;

  loc_stiff_matrix<1,1>(n_knots, deg);
  loc_stiff_matrix<2,1>(n_knots, deg);
  loc_stiff_matrix<3,1>(n_knots, deg);

  loc_stiff_matrix<2,2>(n_knots, deg);
  loc_stiff_matrix<2,3>(n_knots, deg);

  loc_stiff_matrix<3,3>(n_knots, deg);

#endif

  return  0;
}
