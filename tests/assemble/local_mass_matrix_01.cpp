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
 *  Test the assembling of local mass matrix using a BSpline basis
 *
 *  author: martinelli
 *  date: Nov 06, 2015
 *
 */

#include "../tests.h"

#include <igatools/basis_functions/physical_basis_handler.h>
#include <igatools/basis_functions/bspline.h>
#include <igatools/basis_functions/bspline_handler.h>
#include <igatools/basis_functions/bspline_element.h>
#include <igatools/base/quadrature_lib.h>

//#include <igatools/linear_algebra/dense_matrix.h>


template<int dim,int range>
void loc_mass_matrix(const int n_knots, const int deg)
{

  OUTSTART

  out.begin_item("local_mass_matrix<" + std::to_string(dim) + "," + std::to_string(range) + ">");

  auto grid = Grid<dim>::create(n_knots);
  auto space = SplineSpace<dim,range>::create(deg, grid)
               auto basis = BSpline<dim,range>::create(space) ;

  auto elem_handler = basis->create_cache_handler();


  using Flags = basis_element::Flags;
  auto flag = Flags::value | Flags::w_measure;
  elem_handler->template set_flags<dim>(flag);

  auto elem           = basis->begin();
  const auto elem_end = basis->end();

  auto quad = QGauss<dim>::create(deg+1);
  elem_handler->init_element_cache(elem,quad);

  int elem_id = 0;
  for (; elem != elem_end; ++elem, ++elem_id)
  {
    out.begin_item("Element " + std::to_string(elem_id));

    out << "Element ID: " << elem->get_index() << std::endl;

    elem_handler->fill_element_cache(elem);

//    auto loc_mat = elem->template integrate_u_v<dim>(0);
    auto loc_mat = elem->integrate_element_u_v();
    out.begin_item("Mass matrix:");
    loc_mat.print_info(out);
    out.end_item();

    out.end_item();
  }

  out.end_item();

  OUTEND

}




int main()
{
  const int n_knots = 6;
  const int deg = 1;

  loc_mass_matrix<1,1>(n_knots, deg);
  loc_mass_matrix<2,1>(n_knots, deg);
  loc_mass_matrix<3,1>(n_knots, deg);

  loc_mass_matrix<2,2>(n_knots, deg);
  loc_mass_matrix<3,3>(n_knots, deg);

  return  0;
}
