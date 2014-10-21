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

#ifndef FORMULA_FUNCTIONS_H
#define FORMULA_FUNCTIONS_H

#include <igatools/base/new_function.h>

IGA_NAMESPACE_OPEN

template<int dim, int codim=0, int range = 1, int rank = 1>
class FormulaFunction : public NewFunction<dim, codim, range, rank>
{
private:
    using parent_t = NewFunction<dim, codim, range, rank>;
public:
    using typename parent_t::Point;
    using typename parent_t::Value;
    using typename parent_t::Gradient;
    using typename parent_t::ElementIterator;
    using typename parent_t::ElementAccessor;
    using parent_t::space_dim;

    template <int order>
    using Derivative = typename parent_t::template Derivative<order>;

    FormulaFunction(std::shared_ptr<const CartesianGrid<dim>> grid,
                    const ValueFlags &flag,
                    const Quadrature<dim> &quad = Quadrature<dim>());

    void reset(const ValueFlags &flag, const Quadrature<dim> &quad);

    void init_elem(ElementAccessor &elem);

    void fill_elem(ElementAccessor &elem);

private:
    virtual void parametrization(const ValueVector<Points<dim>> &points_,
                                 ValueVector<Point> &values) const
    {
        const int num_points = points_.size();
        for (int i = 0; i<num_points; i++)
        {
            const auto &x = points_[i];
            for (int k = 0; k < dim; ++k)
            {
                values[i][k] = x[k];
            }
            for (int k = dim; k < codim; ++k)
            {
                values[i][k] = 0.;
            }
        }
    }

    virtual void evaluate_0(const ValueVector<Point> &points,
                            ValueVector<Value> &values) const = 0;

    virtual void evaluate_1(const ValueVector<Point> &points,
                            ValueVector<Derivative<1>> &values) const = 0;

    virtual void evaluate_2(const ValueVector<Point> &points,
                            ValueVector<Derivative<2>> &values) const = 0;

private:
    FunctionFlags flag_;
    Quadrature<dim> quad_;
};

IGA_NAMESPACE_CLOSE

#endif
