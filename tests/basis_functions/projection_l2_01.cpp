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
 *  Test for the l2_projection function.
 *  Bspline spaces case
 *
 *  author: pauletti
 *  date: 2013-10-10
 */

#include "../tests.h"
#include "./common_functions.h"

#include <igatools/base/quadrature_lib.h>
#include <igatools/geometry/formula_grid_function.h>

#include <igatools/basis_functions/bspline.h>
#include <igatools/basis_functions/bspline_element.h>

#include <igatools/basis_functions/space_tools.h>





template<int dim,int range>
class TestGridFunction : public FormulaGridFunction<dim,range>
{
private:
  using base_t = GridFunction<dim,range>;
  using parent_t = FormulaGridFunction<dim,range>;
  using self_t = TestGridFunction<dim,range>;
public:
  using typename parent_t::Value;
  using typename parent_t::GridPoint;

  template <int order>
  using Derivative = typename parent_t::template Derivative<order>;

public:
  TestGridFunction(const SharedPtrConstnessHandler<Grid<dim>> &grid)
    : parent_t(grid)
  {}

  static std::shared_ptr<self_t>
  const_create(const std::shared_ptr<const Grid<dim>> &grid)
  {
    return std::shared_ptr<self_t>(new self_t(SharedPtrConstnessHandler<Grid<dim>>(grid)));
  }

  Real value(const Points<dim> &x) const
  {
    Real f = 1;
    for (int i = 0; i<dim; ++i)
      f = f * cos(2*PI*x[i]);
    return f;
  }

  void evaluate_0(const ValueVector<GridPoint> &points,
                  ValueVector<Value> &values) const override final
  {
    for (int i = 0; i<points.size(); ++i)
    {
      Points<dim> p = points[i];
      values[i][0] = this->value(p);
    }
  }
  void evaluate_1(const ValueVector<GridPoint> &points,
                  ValueVector<Derivative<1>> &values) const override final
  {
    Assert(false,ExcNotImplemented());
  }

  void evaluate_2(const ValueVector<GridPoint> &points,
                  ValueVector<Derivative<2>> &values) const override final
  {
    Assert(false,ExcNotImplemented());
  }

  void print_info(LogStream &out) const
  {
    Assert(false,ExcNotImplemented());
  }

  void rebuild_after_insert_knots(
    const iga::SafeSTLArray<iga::SafeSTLVector<double>, dim> &new_knots, const iga::Grid<dim> &g)
  {
    Assert(false,ExcNotImplemented());
  }
};


template<int dim , int range>
void project_l2(const int p, const int num_knots = 10)
{
  OUTSTART

  auto knots = Grid<dim>::const_create(num_knots);
  auto space = SplineSpace<dim,range>::const_create(p, knots) ;
  auto basis = BSpline<dim,range>::const_create(space) ;


  const int n_qpoints = 4;
  auto quad = QGauss<dim>::const_create(n_qpoints);

  auto f = TestGridFunction<dim,range>::const_create(knots);
  auto coeffs_func = space_tools::projection_l2_grid_function<dim,range>(*f,*basis,quad);

  auto proj_func = IgGridFunction<dim,range>::const_create(basis,coeffs_func);
  proj_func->print_info(out);

  OUTEND
}



int main()
{
//  project_l2<0,1>(1);
  project_l2<1,1>(3);
  project_l2<2,1>(3);
  project_l2<3,1>(1);

  return 0;
}

