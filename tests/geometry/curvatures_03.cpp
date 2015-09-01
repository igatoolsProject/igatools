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
 *  Measure of a physical domain (mapping)
 *
 *  author: pauletti
 *  date: 2015-08-05
 *
 */

#include "../tests.h"

#include <igatools/functions/function_lib.h>
#include <igatools/base/quadrature_lib.h>
#include <igatools/functions/function_element.h>
#include <igatools/functions/function_lib.h>
#include <igatools/functions/identity_function.h>
#include <igatools/geometry/mapping.h>
#include <igatools/geometry/mapping_element.h>
#include <igatools/basis_functions/space_tools.h>


template <int dim>
void intregrate_on_sphere(const int n_knots)
{
  const int range=1;
  OUTSTART

  using MapFunction       = functions::SphereFunction<dim>;
  using IntegrandFunction = functions::ConstantFunction<dim,1,1,1>;

  BBox<dim> box;
  for (int i=0; i<dim-1; ++i)
    box[i] = {0., M_PI/2.};
  if (dim>=1)
    box[dim-1] = {0., M_PI/2.};
  auto grid = CartesianGrid<dim>::create(box, n_knots);

  auto F = MapFunction::create(grid, IdentityFunction<dim>::create(grid));
  typename IntegrandFunction::Value val {1.};
  auto C = IntegrandFunction::create(grid, F, val);


  auto quad = QGauss<dim>(3);

  SafeSTLVector<typename IntegrandFunction::Value> vec(F->get_grid()->get_num_all_elems());
  auto area = space_tools::integrate<0,dim, 1, range, 1 >(*C, quad, vec);
  vec.print_info(out);
  out << endl;


  //n-sphere volume
  //auto exact = std::pow(M_PI, dim/2.) / tgamma(dim/ 2. + 1);

  //n-sphere area
  auto exact = 2. * std::pow(M_PI, (dim+1)/2.) / tgamma((dim+1)/ 2.);
  auto ans = exact * 0.5 * std::pow(0.25, dim-1.);

  out << exact<< "    "<< ans << "    " << area << endl;
//    SafeSTLVector<typename IntegrandFunction::template Derivative<1>> vec_der(F->get_grid()->get_num_all_elems());
//    space_tools::integrate<1, dim, 1, range, 1 >(*C,  quad, vec_der);
//    vec_der.print_info(out);
//    out << endl;

  OUTEND
}


int main()
{
  out.depth_console(10);

  //intregrate_on_sphere<1>(4);
  intregrate_on_sphere<2>(2);

  return 0;
}
