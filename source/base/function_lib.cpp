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

template<int dim, int codim, int range, int rank>
ConstantFunction<dim, codim, range, rank>::
ConstantFunction(std::shared_ptr<const CartesianGrid<dim>> grid,
                 const Value &b,
                 const NewValueFlags &flag,
                 const Quadrature<dim> &quad)
    :
    parent_t::FormulaFunction(grid, flag, quad),
    b_(b)
{}



template<int dim, int codim, int range, int rank>
auto
ConstantFunction<dim, codim, range, rank>::
create(std::shared_ptr<const CartesianGrid<dim>> grid,
       const Value &b,
       const NewValueFlags &flag,
       const Quadrature<dim> &quad) ->  std::shared_ptr<base_t>
{
    return std::shared_ptr<base_t>(new self_t(grid, b, flag, quad));
}



template<int dim, int codim, int range, int rank>
auto
ConstantFunction<dim, codim, range, rank>::
evaluate_0(const ValueVector<Point> &points,
           ValueVector<Value> &values) const -> void
{
    for (auto &val : values)
        val =  b_;
}



template<int dim, int codim, int range, int rank>
auto
ConstantFunction<dim, codim, range, rank>::
evaluate_1(const ValueVector<Point> &points,
           ValueVector<Derivative<1>> &values) const -> void
{
    for (auto &val : values)
            val = 0.;
}



template<int dim, int codim, int range, int rank>
auto
ConstantFunction<dim, codim, range, rank>::
evaluate_2(const ValueVector<Point> &points,
           ValueVector<Derivative<2>> &values) const -> void
{
    for (auto &val : values)
        val = 0.;
}



//------------------------------------------------------------------------------
template<int dim, int codim, int range>
LinearFunction<dim, codim, range>::
LinearFunction(std::shared_ptr<const CartesianGrid<dim>> grid,
               const Gradient &A, const Value &b,
               const NewValueFlags &flag, const Quadrature<dim> &quad)
    :
    parent_t::FormulaFunction(grid, flag, quad),
    A_(A),
    b_(b)
{}



template<int dim, int codim, int range>
auto
LinearFunction<dim, codim, range>::
create(std::shared_ptr<const CartesianGrid<dim>> grid,
       const Gradient &A, const Value &b,
       const NewValueFlags &flag,
       const Quadrature<dim> &quad) ->  std::shared_ptr<base_t>
{
    return std::shared_ptr<base_t>(new self_t(grid, A, b, flag, quad));
}



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

} // of namespace functions.

IGA_NAMESPACE_CLOSE

#include <igatools/base/function_lib.inst>
