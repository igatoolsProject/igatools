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
 *  Test the assembling of local stiffness matrix
 *
 *  author: pauletti
 *  date: Oct 09, 2014
 *
 */

#include "../tests.h"

#include <igatools/basis_functions/physical_space.h>
#include <igatools/basis_functions/phys_space_element_handler.h>
#include <igatools/basis_functions/bspline_space.h>
#include <igatools/basis_functions/physical_space_element.h>
//#include <igatools/basis_functions/space_uniform_quad_cache.h>
#include <igatools/base/quadrature_lib.h>
#include <igatools/functions/identity_function.h>

#include <igatools/linear_algebra/dense_matrix.h>


template<int dim>
void loc_stiff_matrix(const int n_knots, const int deg)
{
  OUTSTART

  auto grid = CartesianGrid<dim>::create(n_knots);
  using Space = BSplineSpace<dim>;
  auto ref_space = Space::create(deg, grid) ;

  using PhysSpace = PhysicalSpace<dim,1,1,0,Transformation::h_grad>;
  auto identity_mapping = IdentityFunction<dim>::create(grid);
  auto phys_space = PhysSpace::create(ref_space, identity_mapping) ;

  auto elem_handler = phys_space->create_cache_handler();

  auto quad = QGauss<dim>(deg+1);
  auto flag = ValueFlags::value | ValueFlags::gradient | ValueFlags::w_measure;
  elem_handler->reset(flag,quad);

  const int n_qpoints =  quad.get_num_points();

  auto elem           = phys_space->begin();
  const auto elem_end = phys_space->end();

  const int n_basis = elem->get_num_basis(DofProperties::active);
  DenseMatrix loc_mat(n_basis,n_basis);

  elem_handler->init_element_cache(elem);
  for (; elem != elem_end; ++elem)
  {
    elem_handler->fill_element_cache(elem);

    loc_mat = 0.0;

    const auto w_meas = elem->template get_w_measures<dim>(0);
    const auto grad = elem->template get_basis<_Gradient,dim>(0,DofProperties::active);

    for (int i=0; i<n_basis; ++i)
    {
      const auto grad_i = grad.get_function_view(i);

      for (int j=0; j<n_basis; ++j)
      {
        const auto grad_j = grad.get_function_view(j);
        for (int qp = 0; qp < n_qpoints; ++qp)
          loc_mat(i,j) += scalar_product(grad_i[qp], grad_j[qp]) * w_meas[qp];
      }
    }
    loc_mat.print_info(out);
    out << endl;
  }

  OUTEND

}




int main()
{
  const int n_knots = 6;
  const int deg = 1;

  loc_stiff_matrix<1>(n_knots, deg);
  loc_stiff_matrix<2>(n_knots, deg);
  loc_stiff_matrix<3>(n_knots, deg);

  return  0;
}
