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
template<int dim, int codim, int range, int rank = 1>
class ConstantFunction : public FormulaFunction<dim, codim, range, rank>
{
private:
    using base_t = NewFunction<dim, codim, range, rank>;
    using parent_t = FormulaFunction<dim, codim, range, rank>;
    using self_t = ConstantFunction<dim, codim, range, rank>;
    using typename base_t::GridType;
public:
    using typename parent_t::Point;
    using typename parent_t::Value;
    using typename parent_t::Gradient;
    using typename parent_t::ElementIterator;
    using typename parent_t::ElementAccessor;
    template <int order>
    using Derivative = typename parent_t::template Derivative<order>;

    static std::shared_ptr<base_t>
    create(std::shared_ptr<GridType> grid,
           const Value &b);

protected:
    ConstantFunction(std::shared_ptr<GridType> grid,
                     const Value &b);

private:
    void evaluate_0(const ValueVector<Point> &points,
                    ValueVector<Value> &values) const;

    void evaluate_1(const ValueVector<Point> &points,
                    ValueVector<Derivative<1>> &values) const;

    void evaluate_2(const ValueVector<Point> &points,
                    ValueVector<Derivative<2>> &values) const;

private:
    const Value b_;
};



//------------------------------------------------------------------------------

template<int dim, int codim, int range>
class LinearFunction : public FormulaFunction<dim, codim, range, 1>
{

public:
    using base_t = NewFunction<dim, codim, range, 1>;
    using parent_t = FormulaFunction<dim, codim, range, 1>;
    using self_t = LinearFunction<dim, codim, range>;
    using typename base_t::GridType;
    using typename parent_t::Point;
    using typename parent_t::Value;
    using typename parent_t::Gradient;
    using typename parent_t::ElementIterator;
    using typename parent_t::ElementAccessor;
    template <int order>
    using Derivative = typename parent_t::template Derivative<order>;

    static std::shared_ptr<base_t>
    create(std::shared_ptr<GridType> grid,
           const Gradient &A,
           const Value &b);

protected:
    LinearFunction(std::shared_ptr<GridType> grid,
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

//------------------------------------------------------------------------------
/**
 * Maps a hyper rectangle into a spherical ball sector using the
 * dim-dimensional spherical coordinates, maps a hyper-rectangle
 * r in [0,R], phi_1 in [0, 2 pi], and phi_2, phi_dim-1 in [0,pi]
 * such that
 * x1 = r cos (phi_1)
 * x2 = r sin (phi_1) cos (phi_2)
 * etc
 *
 */
template<int dim>
class BallFunction : public FormulaFunction<dim, 0, dim, 1>
{
private:
    using base_t = NewFunction<dim, 0, dim, 1>;
    using parent_t = FormulaFunction<dim, 0, dim, 1>;
    using self_t = BallFunction<dim>;
    using typename base_t::GridType;
public:
    using typename parent_t::Point;
    using typename parent_t::Value;
    using typename parent_t::Gradient;
    using typename parent_t::ElementIterator;
    using typename parent_t::ElementAccessor;
    template <int order>
    using Derivative = typename parent_t::template Derivative<order>;

    static std::shared_ptr<base_t>
    create(std::shared_ptr<GridType> grid);

protected:
    BallFunction(std::shared_ptr<GridType> grid);

private:
    template<int order>
    auto
    get_aux_vals(const ValueVector<Point> &points) const;

private:
    void evaluate_0(const ValueVector<Point> &points,
                    ValueVector<Value> &values) const;

    void evaluate_1(const ValueVector<Point> &points,
                    ValueVector<Derivative<1>> &values) const;

    void evaluate_2(const ValueVector<Point> &points,
                    ValueVector<Derivative<2>> &values) const;

private:
    const Value b_;
};



} // of namespace functions.

IGA_NAMESPACE_CLOSE

#endif
