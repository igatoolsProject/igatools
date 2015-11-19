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
 *  Measure of a physical domain using space_tools::integrate
 *
 *  author: pauletti
 *  date: 2015-08-05
 *
 */

#include "../tests.h"

#include <igatools/geometry/grid_function_lib.h>
#include <igatools/base/quadrature_lib.h>
#include <igatools/functions/function_element.h>
#include <igatools/functions/function_lib.h>
#include <igatools/geometry/domain.h>
#include <igatools/geometry/domain_element.h>
#include <igatools/basis_functions/space_tools.h>


template <int dim>
void intregrate_on_sphere(const int n_knots)
{
  const int range=1;
  OUTSTART

  out.begin_item("integrate_on_sphere<dim="
                 + std::to_string(dim) + ">(n_elems_1D="
                 + std::to_string(n_knots-1) + ")");

  using Sph = grid_functions::SphereGridFunction<dim>;
  using IntegrandFunction = functions::ConstantFunction<dim,1,1,1>;

  BBox<dim> box;
  for (int i=0; i<dim-1; ++i)
    box[i] = {0., M_PI/2.};
  if (dim>=1)
    box[dim-1] = {0., M_PI/2.};
  auto grid = Grid<dim>::create(box, n_knots);

  auto sph_func = Sph::create(grid);
  auto sph_domain = Domain<dim,1>::create(sph_func);
  typename IntegrandFunction::Value val {1.};
  auto C = IntegrandFunction::create(sph_domain, val);


  auto quad = QGauss<dim>::create(2);

  SafeSTLMap<ElementIndex<dim>,
             typename IntegrandFunction::Value> elem_contrib;
  auto area = space_tools::integrate<0,dim,1,range,1>(*C, quad, elem_contrib);
  out.begin_item("Area: contribution from the elements");
  elem_contrib.print_info(out);
  out.end_item();


  //n-sphere volume
  //auto exact = std::pow(M_PI, dim/2.) / tgamma(dim/ 2. + 1);

  //n-sphere area
  auto exact = 2. * std::pow(M_PI, (dim+1)/2.) / tgamma((dim+1)/ 2.);
  auto ans = exact * 0.5 * std::pow(0.25, dim-1.);

  out.begin_item("Area:");
  out << "Exact      = " << ans << endl;
  out << "Quadrature = " << area << endl;
  out.end_item();

  out.end_item();

  OUTEND
}


int main()
{
  out.depth_console(10);


  for (int n_knots = 2 ; n_knots <= 5 ; ++n_knots)
  {
//  intregrate_on_sphere<1>(n_knots);
    intregrate_on_sphere<2>(n_knots);
  }

  return 0;
}
