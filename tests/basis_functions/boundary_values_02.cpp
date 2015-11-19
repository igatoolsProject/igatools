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
 *  Test for the boundary l2 projection function.
 *  On a BSpline (a reference space)
 *
 *  author: pauletti
 *  date: 2014-11-14
 *
 */

#include "../tests.h"
#include <igatools/basis_functions/space_tools.h>
#include <igatools/linear_algebra/dof_tools.h>

#include "common_functions.h"

#include <igatools/base/quadrature_lib.h>
#include <igatools/basis_functions/bspline.h>

#include <igatools/geometry/grid_function_lib.h>





template<int dim , int range ,int rank>
void do_test(const int p, const int num_knots = 10)
{
  out << "Dimension: " << dim << endl;


  auto grid = Grid<dim>::create(num_knots);
  auto space = SplineSpace<dim,range,rank>::create(p,grid);
  auto ref_basis = BSpline<dim,range,rank>::create(space) ;
  auto map = grid_functions::IdentityGridFunction<dim>::create(grid);
  auto domain = Domain<dim,0>::create(map);
  auto basis = PhysicalSpaceBasis<dim,range,rank,0>::create(ref_basis, domain);


  const int sdim = dim-1;
  const int s_id = 0;

  using SubGridElemMap = typename Grid<dim>::template SubGridMap<sdim>;
  SubGridElemMap sub_grid_elem_map;
  const std::shared_ptr<const Grid<sdim>> sub_grid = grid->template get_sub_grid<sdim>(s_id,sub_grid_elem_map);

  auto bndry_domain = domain->get_sub_domain(s_id,sub_grid_elem_map,sub_grid);

  using BndFunc = Function<dim-1,1,range,1>;
  SafeSTLMap<int,std::shared_ptr<const BndFunc>> boundary_functions;
  boundary_functions[s_id] = TestBoundaryFunction<dim-1,range>::const_create(bndry_domain);


  const int n_qpoints = 4;
  auto quad = QGauss<sdim>::const_create(n_qpoints);


  std::map<Index,Real> boundary_values;
  space_tools::project_boundary_values<dim,0,range,rank>(
    boundary_functions, *basis, quad, boundary_values);

  out << "basis index \t value" << endl;
  for (auto entry : boundary_values)
    out << entry.first << "\t" << entry.second << endl;

}



int main()
{

  // do_test<1,1,1>(3);

  //TODO (martinelli, Nov 18,2015): with gcc.5.2.0 it crash calling the two function in this order
  do_test<2,1,1>(3);
  do_test<3,1,1>(2);


//  do_test<2,1,1>(3);

  return 0;
}

