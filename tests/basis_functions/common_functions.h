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


// TODO (pauletti, Nov 13, 2014): delete this after cleaning and arranging test
using numbers::PI;

template<int dim>
class BoundaryFunction : public FormulaFunction<dim>
{
private:
  using base_t = Function<dim>;
  using parent_t = FormulaFunction<dim>;
  using self_t = BoundaryFunction<dim>;
  using typename parent_t::DomainType;

public:
  using typename parent_t::Point;
  using typename parent_t::Value;
  template <int order>
  using Derivative = typename parent_t::template Derivative<order>;
public:
  BoundaryFunction(SharedPtrConstnessHandler<DomainType> domain)
    : FormulaFunction<dim>(domain)
  {}

  static std::shared_ptr<const base_t>
  const_create(std::shared_ptr<const DomainType> domain)
  {
    return std::shared_ptr<const self_t>(new self_t(SharedPtrConstnessHandler<DomainType>(domain)));
  }


  Real value(Points<dim> x) const
  {
    Real f = 1;
    for (int i = 0; i<dim; ++i)
      f = f * cos(2*PI*x[i]);
    return f;
  }

  void evaluate_0(const ValueVector<Point> &points,
                  ValueVector<Value> &values) const override
  {
    for (int i = 0; i<points.size(); ++i)
    {
      Points<dim> p = points[i];
      values[i][0] = this->value(p);
    }
  }
  void evaluate_1(const ValueVector<Point> &points,
                  ValueVector<Derivative<1>> &values) const override
  {}

  void evaluate_2(const ValueVector<Point> &points,
                  ValueVector<Derivative<2>> &values) const override
  {}
};



/**
 * The p-norm of this function on the unit square is
 * (1/(p+1))^(n/p)
 */
template<int dim, int codim=0, int range = 1, int rank = 1>
class ProductFunction : public FormulaFunction<dim>
{
private:
  using base_t = Function<dim, codim, range, rank>;
  using parent_t = FormulaFunction<dim, codim, range, rank>;
  using self_t = ProductFunction<dim, codim, range, rank>;
public:
  using typename parent_t::Point;
  using typename parent_t::Value;
  using typename parent_t::Gradient;
  using typename parent_t::ElementIterator;
  using typename parent_t::ElementAccessor;
  template <int order>
  using Derivative = typename parent_t::template Derivative<order>;

  using parent_t::FormulaFunction;



private:
  void evaluate_0(const ValueVector<Point> &points,
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

  void evaluate_1(const ValueVector<Point> &points,
                  ValueVector<Derivative<1>> &values) const
  {}

  void evaluate_2(const ValueVector<Point> &points,
                  ValueVector<Derivative<2>> &values) const
  {}

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





