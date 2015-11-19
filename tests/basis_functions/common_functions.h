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
 *  QA: The 0 dim case not returning appropriate value
 */


#include <igatools/base/quadrature_lib.h>
//#include <igatools/functions/identity_function.h>
#include <igatools/functions/formula_function.h>
#include <igatools/functions/function_lib.h>

#include <igatools/geometry/formula_grid_function.h>


// TODO (pauletti, Nov 13, 2014): delete this after cleaning and arranging test
using numbers::PI;




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

  static std::shared_ptr<const self_t>
  const_create(const std::shared_ptr<const Grid<dim>> &grid)
  {
    return std::shared_ptr<self_t>(new self_t(SharedPtrConstnessHandler<Grid<dim>>(grid)));
  }

  static std::shared_ptr<self_t>
  create(const std::shared_ptr<Grid<dim>> &grid)
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

  void print_info(LogStream &out) const override
  {
    Assert(false,ExcNotImplemented());
  }

  void rebuild_after_insert_knots(
    const iga::SafeSTLArray<iga::SafeSTLVector<double>, dim> &new_knots, const iga::Grid<dim> &g) override
  {
    Assert(false,ExcNotImplemented());
  }
};


template<int dim,int range>
class TestFunction : public FormulaFunction<dim,0,range,1>
{
private:
  using base_t = Function<dim,0,range,1>;
  using parent_t = FormulaFunction<dim,0,range,1>;
  using self_t = TestFunction<dim,range>;
public:
  using typename parent_t::Value;
  using typename parent_t::Point;

  template <int order>
  using Derivative = typename parent_t::template Derivative<order>;

public:
  TestFunction(const SharedPtrConstnessHandler<Domain<dim,0>> &domain)
    : parent_t(domain,"TestFunction")
  {}

  static std::shared_ptr<const self_t>
  const_create(const std::shared_ptr<const Domain<dim,0>> &domain)
  {
    return std::shared_ptr<self_t>(new self_t(SharedPtrConstnessHandler<Domain<dim,0>>(domain)));
  }

  static std::shared_ptr<self_t>
  create(const std::shared_ptr<Domain<dim,0>> &domain)
  {
    return std::shared_ptr<self_t>(new self_t(SharedPtrConstnessHandler<Domain<dim,0>>(domain)));
  }

  Real value(const Points<dim> &x) const
  {
    Real f = 1;
    for (int i = 0; i<dim; ++i)
      f = f * cos(2*PI*x[i]);
    return f;
  }

  void evaluate_0(const ValueVector<Point> &points,
                  ValueVector<Value> &values) const override final
  {
    for (int i = 0; i<points.size(); ++i)
      values[i][0] = this->value(points[i]);
  }

  void evaluate_1(const ValueVector<Point> &points,
                  ValueVector<Derivative<1>> &values) const override final
  {
    Assert(false,ExcNotImplemented());
  }

  void evaluate_2(const ValueVector<Point> &points,
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


template<int dim,int range>
class TestBoundaryFunction : public FormulaFunction<dim,1,range,1>
{
private:
  using base_t = Function<dim,1,range,1>;
  using parent_t = FormulaFunction<dim,1,range,1>;
  using self_t = TestBoundaryFunction<dim,range>;
public:
  using typename parent_t::Value;
  using typename parent_t::Point;

  template <int order>
  using Derivative = typename parent_t::template Derivative<order>;

public:
  TestBoundaryFunction(const SharedPtrConstnessHandler<Domain<dim,1>> &domain)
    : parent_t(domain,"TestBoundaryFunction")
  {}

  static std::shared_ptr<const self_t>
  const_create(const std::shared_ptr<const Domain<dim,1>> &domain)
  {
    return std::shared_ptr<self_t>(new self_t(SharedPtrConstnessHandler<Domain<dim,1>>(domain)));
  }

  static std::shared_ptr<self_t>
  create(const std::shared_ptr<Domain<dim,1>> &domain)
  {
    return std::shared_ptr<self_t>(new self_t(SharedPtrConstnessHandler<Domain<dim,1>>(domain)));
  }

  Real value(const Points<dim+1> &x) const
  {
    Real f = 1;
    for (int i = 0; i<dim+1; ++i)
      f = f * cos(2*PI*x[i]);
    return f;
  }

  void evaluate_0(const ValueVector<Point> &points,
                  ValueVector<Value> &values) const override final
  {
    for (int i = 0; i<points.size(); ++i)
      values[i][0] = this->value(points[i]);
  }

  void evaluate_1(const ValueVector<Point> &points,
                  ValueVector<Derivative<1>> &values) const override final
  {
    Assert(false,ExcNotImplemented());
  }

  void evaluate_2(const ValueVector<Point> &points,
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


/**
 * The p-norm of this function on the unit square is
 * (1/(p+1))^(n/p)
 */
template<int dim,int range = 1>
class ProductGridFunction : public FormulaGridFunction<dim,range>
{
private:
  using base_t = GridFunction<dim,range>;
  using parent_t = FormulaGridFunction<dim,range>;
  using self_t = ProductGridFunction<dim,range>;
public:
  using typename parent_t::GridPoint;
  using typename parent_t::Value;
  using typename parent_t::Gradient;
  using typename parent_t::ElementIterator;
  using typename parent_t::ElementAccessor;
  template <int order>
  using Derivative = typename parent_t::template Derivative<order>;

  using parent_t::FormulaGridFunction;

  static std::shared_ptr<const self_t>
  const_create(const std::shared_ptr<const Grid<dim>> &grid)
  {
    return std::shared_ptr<self_t>(new self_t(SharedPtrConstnessHandler<Grid<dim>>(grid)));
  }


private:
  void evaluate_0(const ValueVector<GridPoint> &points,
                  ValueVector<Value> &values) const
  {
    auto pt = points.begin();
    auto val = values.begin();

    for (; pt != points.end(); ++pt, ++val)
    {
      *val = 1.;
      for (int i=0; i<dim; ++i)
        (*val) *= (*pt)[i];
    }
  }

  void evaluate_1(const ValueVector<GridPoint> &points,
                  ValueVector<Derivative<1>> &values) const
  {
    Assert(false,ExcNotImplemented());
  }

  void evaluate_2(const ValueVector<GridPoint> &points,
                  ValueVector<Derivative<2>> &values) const
  {
    Assert(false,ExcNotImplemented());
  }

  void print_info(LogStream &out) const override
  {
    Assert(false,ExcNotImplemented());
  }

  void rebuild_after_insert_knots(
    const iga::SafeSTLArray<iga::SafeSTLVector<double>, dim> &new_knots, const iga::Grid<dim> &g) override
  {
    Assert(false,ExcNotImplemented());
  }

};


/**
 * Norm Function
 * F(x) = (sum x_i ^ p)^(1/p)
 */
template<int dim>
class NormFunction : public FormulaFunction<dim, 0, 1, 1>
{

public:
  using base_t = Function<dim, 0, 1, 1>;
  using parent_t = FormulaFunction<dim, 0, 1, 1>;
  using self_t = NormFunction<dim>;
  using GridType = Grid<dim>;
  using typename parent_t::DomainType;

  using typename parent_t::Point;
  using typename parent_t::Value;
  using typename parent_t::Gradient;
  using typename parent_t::ElementIterator;
  using typename parent_t::ElementAccessor;
  template <int order>
  using Derivative = typename parent_t::template Derivative<order>;

  static std::shared_ptr<const base_t>
  const_create(std::shared_ptr<const DomainType> &domain,
               const Real p=2.)
  {
    return std::shared_ptr<const self_t>(new self_t(SharedPtrConstnessHandler<DomainType>(domain), p));
  }

  NormFunction(const self_t &) = default;

protected:
  NormFunction(SharedPtrConstnessHandler<DomainType> domain,
               const Real p)
    :
    parent_t(domain),
    p_(p)
  {}



private:
  void evaluate_0(const ValueVector<Point> &points,
                  ValueVector<Value> &values) const override
  {
    auto pt = points.begin();
    for (auto &val : values)
    {
      val = std::pow(pt->norm_square(), p_/2.);
      ++pt;
    }
  }

  void evaluate_1(const ValueVector<Point> &points,
                  ValueVector<Derivative<1>> &values) const override
  {}

  void evaluate_2(const ValueVector<Point> &points,
                  ValueVector<Derivative<2>> &values) const override
  {}

private:
  const Real p_;
};





