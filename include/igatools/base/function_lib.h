//-+--------------------------------------------------------------------
// Igatools a general purpose Isogeometric analysis library.
// Copyright (C) 2012-2014  by the igatools authors (see authors.txt).
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

#ifndef __IGA_FUNCTION_LIB_H_
#define __IGA_FUNCTION_LIB_H_

#include <igatools/base/formula_function.h>

IGA_NAMESPACE_OPEN

/**
 * Collection of useful functions derived from the Function class.
 */
namespace functions
{

template<int dim, int codim, int range>
class LinearFunction : public FormulaFunction<dim, codim, range, 1>
{
public:
    using base_t = NewFunction<dim, codim, range, 1>;
    using parent_t = FormulaFunction<dim, codim, range, 1>;
    using self_t = LinearFunction<dim, codim, range>;
    using typename parent_t::Point;
    using typename parent_t::Value;
    using typename parent_t::Gradient;
    using typename parent_t::ElementIterator;
    using typename parent_t::ElementAccessor;
    template <int order>
    using Derivative = typename parent_t::template Derivative<order>;

    LinearFunction(std::shared_ptr<const CartesianGrid<dim>> grid,
                   const NewValueFlags &flag, const Quadrature<dim> &quad,
                   const Gradient &A, const Value &b);

    static std::shared_ptr<base_t>
    create(std::shared_ptr<const CartesianGrid<dim>> grid,
           const NewValueFlags &flag, const Quadrature<dim> &quad,
           const Gradient &A, const Value &b);
private:
    void evaluate_0(const ValueVector<Point> &points,
                    ValueVector<Value> &values) const;

    void evaluate_1(const ValueVector<Point> &points,
                    ValueVector<Derivative<1>> &values) const;

    void evaluate_2(const ValueVector<Point> &points,
                    ValueVector<Derivative<2>> &values) const;

private:
    const Gradient A_;
    const Value    b_;
};

#if 0
//------------------------------------------------------------------------------

/**
 * Constant scalar function.
 */
template<int dim, int range = 1, int rank = 1>
class ConstantFunction
    : public Function<dim, range, rank>
{
public:
    /**
     * Type for the input argument of the function.
     */
    using Point = typename Function<dim, range, rank>::Point;

    /**
     * Type for the return of the function.
     */
    using Value = typename Function<dim, range, rank>::Value;

    /**
     * Construct a constant function with the given value.
     */
    ConstantFunction(const Value value);

    /**
     * Destructor.
     */
    virtual ~ConstantFunction();

    /**
     * Compute the @p values of Function at some @p points.
     */
    void
    evaluate(const ValueVector<Point> &points,
             ValueVector<Value> &values) const;

private:
    /** Constant given value that defines the function. */
    const Value value_;

};

#endif

} // of namespace functions.

IGA_NAMESPACE_CLOSE

#endif /* __IGA_FUNCTION_LIB_H_ */
