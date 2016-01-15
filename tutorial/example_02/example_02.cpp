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

#include <igatools/basis_functions/bspline.h>
#include <igatools/base/logstream.h>

// [headers]
#include <igatools/basis_functions/nurbs.h>
#include <igatools/functions/grid_function_lib.h>
#include <igatools/functions/ig_function.h>
#include <igatools/io/writer.h>
// [headers]

using namespace iga;
using namespace std;

LogStream out;

// [hypercube]
template <int dim>
shared_ptr<const Domain<dim>> reference_domain(Size nel)
{
  auto grid = Grid<dim>::const_create(nel+1);
  auto id_funct = grid_functions::IdentityGridFunction<dim>::const_create(grid);
  return Domain<dim>::const_create(id_funct);
}
// [hypercube]

// [annulus_init]
shared_ptr<Domain<2>> quarter_annulus(TensorSize<2> nel)
{

  TensorIndex<2> deg = {1,2};
  IgCoefficients control_points;
  IgCoefficients weights;
  control_points[ 0] = 1.0;
  control_points[ 6] = 0.0;
  control_points[ 1] = 2.0;
  control_points[ 7] = 0.0;
  control_points[ 2] = 1.0;
  control_points[ 8] = 1.0;
  control_points[ 3] = 2.0;
  control_points[ 9] = 2.0;
  control_points[ 4] = 0.0;
  control_points[10] = 1.0;
  control_points[ 5] = 0.0;
  control_points[11] = 2.0;
  weights[0] = 1.0;
  weights[1] = 1.0;
  weights[2] = sqrt(2.0)/2.0;
  weights[3] = sqrt(2.0)/2.0;
  weights[4] = 1.0;
  weights[5] = 1.0;

  auto grid         = Grid<2>::create(2);
  // [annulus_init]

  // [weight_funct]
  auto scal_space   = SplineSpace<2>::create(deg, grid);
  auto scal_bspline = BSpline<2>::create(scal_space);
  auto weight_funct = IgGridFunction<2,1>::create(scal_bspline, weights);
  // [weight_funct]

  // [geom_funct]
  auto vect_space   = SplineSpace<2,2>::create(deg, grid);
  auto vect_bspline = BSpline<2,2>::create(vect_space);
  auto vect_nurbs   = NURBS<2,2>::create(vect_bspline, weight_funct);
  auto geom_funct   = IgGridFunction<2,2>::create(vect_nurbs, control_points);
  // [geom_funct]

  // [refinement]
  grid->refine_directions({true,true},nel);
  return Domain<2>::create(geom_funct);
}
// [refinement]

// [main]
int main()
{
  auto square  = reference_domain<2>(8);
  auto cube    = reference_domain<3>(4);
  auto annulus = quarter_annulus({2,2});
  // [main]

  // [cube_plot]
  int cube_plot_points = 2;
  Writer<3> writer_cube(cube,cube_plot_points);
  writer_cube.save("cube");
  // [cube_plot]

  // [basis]
  auto grid  = annulus->get_grid_function()->get_grid();
  auto space = SplineSpace<2>::const_create(3,grid);
  auto ref_basis = BSpline<2>::const_create(space);
  auto phy_basis = PhysicalBasis<2>::const_create(ref_basis,annulus);
  // [basis]

  // [basis_funct]
  auto dof_distribution = phy_basis->get_dof_distribution();
  auto coefficients = IgCoefficients(dof_distribution->get_global_dofs());
  TensorIndex<2> tensor_id = {1,1};
  auto flat_id = dof_distribution->get_global_dof_id(tensor_id,0);
  coefficients[flat_id] = 1.0;
  auto basis_funct  = IgFunction<2,0,1,1>::const_create(phy_basis,coefficients);
  // [basis_funct]

  // [annulus_plot]
  int annulus_plot_points = 10;
  Writer<2> writer_annulus(annulus,annulus_plot_points);
  const string fieldname = "basis function sample";
  writer_annulus.add_field(*basis_funct,fieldname);
  writer_annulus.save("annulus");

  return 0;
}
// [annulus_plot]
