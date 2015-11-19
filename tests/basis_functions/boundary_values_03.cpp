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
 *  Test for the boundary projection function.
 *  This test ....
 *  author: pauletti
 *  date: 2013-03-19
 *
 */
//TODO: this test should be merge into the other ones

#include "../tests.h"
#include <igatools/basis_functions/space_tools.h>
#include <igatools/linear_algebra/dof_tools.h>

#include "common_functions.h"

#include <igatools/base/quadrature_lib.h>
#include <igatools/basis_functions/bspline.h>

#include <igatools/geometry/grid_function_lib.h>

template<int dim>
class XProject : public FormulaFunction<dim>
{
private:
  using base_t = Function<dim>;
  using parent_t = FormulaFunction<dim>;
  using self_t = XProject<dim>;
public:
  using typename parent_t::Point;
  using typename parent_t::Value;
  template <int order>
  using Derivative = typename parent_t::template Derivative<order>;
public:
  XProject(const SharedPtrConstnessHandler<Domain<dim,0>> &domain)
    : FormulaFunction<dim>(domain,"XProject")
  {}

  static std::shared_ptr<self_t>
  const_create(const std::shared_ptr<Domain<dim,0>> &domain)
  {
    return std::make_shared<self_t>(
             SharedPtrConstnessHandler<Domain<dim,0>>(domain));
  }

  void evaluate_0(const ValueVector<Point> &points,
                  ValueVector<Value> &values) const override
  {
    for (int i = 0; i<points.size(); ++i)
    {
      Points<dim> p = points[i];
      values[i][0] = p[0];
    }
  }
  void evaluate_1(const ValueVector<Point> &points,
                  ValueVector<Derivative<1>> &values) const override
  {}

  void evaluate_2(const ValueVector<Point> &points,
                  ValueVector<Derivative<2>> &values) const override
  {}
};



template<int dim , int range ,int rank>
void do_test(const int p, TensorSize<dim> n_knots)
{
  out << "Dimension: " << dim << endl;

  /*
    auto grid = Grid<dim>::create(n_knots);
    auto space = Basis::create(SplineSpace<dim,range,rank>::const(p, grid)) ;
    auto f = XProject<dim>::const_create(grid);
  //*/

  auto grid = Grid<dim>::create(n_knots);
  auto space = SplineSpace<dim,range,rank>::create(p,grid);
  auto ref_basis = BSpline<dim,range,rank>::create(space) ;
  auto map = grid_functions::IdentityGridFunction<dim>::create(grid);
  auto domain = Domain<dim,0>::create(map);
  auto basis = PhysicalSpaceBasis<dim,range,rank,0>::create(ref_basis, domain);


  const int sdim = dim-1;
  const int s_id = 2;

  /*
  using SubGridElemMap = typename Grid<dim>::template SubGridMap<sdim>;
  SubGridElemMap sub_grid_elem_map;
  const std::shared_ptr<const Grid<sdim>> sub_grid = grid->template get_sub_grid<sdim>(s_id,sub_grid_elem_map);

  auto bndry_domain = domain->get_sub_domain(s_id,sub_grid_elem_map,sub_grid);
  auto f_at_bndry = TestBoundaryFunction<dim-1,range>::const_create(bndry_domain);
  //*/

  auto f = XProject<dim>::const_create(domain);

  const int n_qpoints = 4;
  auto quad = QGauss<sdim>::create(n_qpoints);

  const boundary_id dirichlet = 1;
  grid->set_boundary_id(s_id, dirichlet);
  SafeSTLSet<boundary_id> bdry_ids;
  bdry_ids.insert(dirichlet);

  std::map<Index,Real> boundary_values;
  space_tools::project_function_on_boundary<dim,0,range,rank>(
    *f, *basis, quad, bdry_ids,
    boundary_values);

  out << "basis index \t value" << endl;
  for (auto entry : boundary_values)
    out << entry.first << "\t" << entry.second << endl;

}


int main()
{
  {
    const int dim = 2;
    TensorSize<dim> n_knots { sequence<dim>(2)};
    do_test<dim, 1, 1>(2, n_knots);
  }
//    do_test<2,1,1>(3);
//    do_test<3,1,1>(2);

  return 0;
}

