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

#include <igatools/base/function_lib.h>
#include <igatools/base/function_element.h>

IGA_NAMESPACE_OPEN

namespace functions
{

template<int dim, int codim, int range>
LinearFunction<dim, codim, range>::
LinearFunction(std::shared_ptr<const CartesianGrid<dim>> grid,
               const ValueFlags &flag, const Quadrature<dim> &quad,
               const Gradient &A, const Value &b)
    :
    parent_t::FormulaFunction(grid, flag, quad),
    A_(A),
    b_(b)
{}

template<int dim, int codim, int range>
auto
LinearFunction<dim, codim, range>::
evaluate_0(const ValueVector<Point> &points,
           ValueVector<Value> &values) const -> void
{
    auto point = points.begin();
    for (auto &val : values)
    {
        val = action(A_, *point) + b_;
        ++point;
    }
}

template<int dim, int codim, int range>
auto
LinearFunction<dim, codim, range>::
evaluate_1(const ValueVector<Point> &points,
           ValueVector<Derivative<1>> &values) const -> void
{
    for (auto &val : values)
        val = A_;
}

template<int dim, int codim, int range>
auto
LinearFunction<dim, codim, range>::
evaluate_2(const ValueVector<Point> &points,
           ValueVector<Derivative<2>> &values) const -> void
{
    for (auto &val : values)
        val = 0.;
}



//------------------------------------------------------------------------------



template<int dim, int range, int rank>
ConstantFunction<dim, range, rank>::
ConstantFunction(const Value value)
    :value_ {value}
{}



template<int dim, int range, int rank>
ConstantFunction<dim, range, rank>::
~ConstantFunction()
{}



template<int dim, int range, int rank>
void
ConstantFunction<dim, range, rank>::
evaluate(
    const ValueVector<Point> &points,
    ValueVector<Value> &values) const
{
    Assert(points.size() == values.size(),
           ExcDimensionMismatch(points.size(), values.size())) ;

    for (auto &value : values)
        value = value_;
}

} // of namespace functions.


IGA_NAMESPACE_CLOSE

#include <igatools/base/function_lib.inst>
