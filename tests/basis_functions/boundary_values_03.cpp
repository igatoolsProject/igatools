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

#include <igatools/functions/grid_function_lib.h>

template<int dim,int codim = 0>
class XProject : public FormulaFunction<dim,codim>
{
private:
  using base_t = Function<dim,codim>;
  using parent_t = FormulaFunction<dim,codim>;
  using self_t = XProject<dim,codim>;
public:
  using typename parent_t::Point;
  using typename parent_t::Value;
  template <int order>
  using Derivative = typename parent_t::template Derivative<order>;
public:
  XProject(const SharedPtrConstnessHandler<Domain<dim,codim>> &domain)
    : FormulaFunction<dim,codim>(domain)
  {}

  static std::shared_ptr<self_t>
  const_create(const std::shared_ptr<Domain<dim,codim>> &domain)
  {
    return std::make_shared<self_t>(
             SharedPtrConstnessHandler<Domain<dim,codim>>(domain));
  }

  void evaluate_0(const ValueVector<Point> &points,
                  ValueVector<Value> &values) const override
  {
    for (int pt = 0; pt < points.size(); ++pt)
    {
      Point p = points[pt];
      values[pt][0] = p[0];
    }
  }

  void evaluate_1(const ValueVector<Point> &points,
                  ValueVector<Derivative<1>> &values) const override
  {}

  void evaluate_2(const ValueVector<Point> &points,
                  ValueVector<Derivative<2>> &values) const override
  {}

  void print_info(LogStream &out) const override
  {}
};



template<int dim , int range ,int rank>
void do_test(const int p, TensorSize<dim> n_knots)
{
  out << "Dimension: " << dim << endl;


  auto grid = Grid<dim>::create(n_knots);
  auto space = SplineSpace<dim,range,rank>::create(p,grid);
  auto ref_basis = BSpline<dim,range,rank>::create(space) ;
  auto map = grid_functions::IdentityGridFunction<dim>::create(grid);
  auto domain = Domain<dim,0>::create(map);
  auto basis = PhysicalBasis<dim,range,rank,0>::create(ref_basis, domain);


  const int sdim = dim-1;
  const int s_id = 2;


  typename Grid<dim>::template SubGridMap<sdim> elem_map;
  auto sub_grid = grid->template get_sub_grid<sdim>(s_id,elem_map);

  using SubDomain = Domain<sdim,1>;
  auto sub_domain = const_pointer_cast<SubDomain>(domain->template get_sub_domain<sdim>(s_id,elem_map,sub_grid));
  auto f = XProject<sdim,1>::const_create(sub_domain);

  const int n_qpoints = 4;
  auto quad = QGauss<sdim>::create(n_qpoints);

  std::map<Index,Real> boundary_values;


  std::map<int, std::shared_ptr<const Function<sdim,1,range,rank>>> bndry_funcs;
  bndry_funcs[s_id] = f;

  space_tools::project_boundary_values<dim,0,range,rank>(
    bndry_funcs,
    *basis,
    quad,
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

