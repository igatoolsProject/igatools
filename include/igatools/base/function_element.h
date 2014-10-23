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

#ifndef FUNCTION_ELEMENT_H
#define FUNCTION_ELEMENT_H

#include <igatools/base/new_function.h>
#include <igatools/geometry/cartesian_grid_element.h>

IGA_NAMESPACE_OPEN

template<int dim, int codim, int range = 1, int rank = 1>
class FunctionElement : public CartesianGridElement<dim>
{
public:
    using Func = NewFunction<dim, codim, range, rank>;
    using Point = typename Func::Point;
    using Value = typename Func::Value;
    using Gradient = typename Func::Gradient;
    using Hessian  = typename Func::Hessian;
    using ContainerType = CartesianGrid<dim>;

private:
    template <int order>
    using Derivative = typename Func::template Derivative<order>;
public:
    using CartesianGridElement<dim>::CartesianGridElement;

    ValueVector<Point> const &get_points() const;

    ValueVector<Value> const &get_values() const;

private:
    template<int order>
    auto const &get_derivative() const;
public:
    ValueVector<Gradient> const &get_gradients() const;

    ValueVector<Hessian> const &get_hessians() const;

private:
    struct Cache : public CacheStatus
    {
        void resize(const FunctionFlags &flags_handler,
                    const int n_points)
        {
            //TODO(pauletti, Oct 11, 2014): missing all necesary clears
            flags_handler_ = flags_handler;

            if (flags_handler_.fill_points())
                points_.resize(n_points);

            if (flags_handler_.fill_values())
                values_.resize(n_points);

            if (flags_handler_.fill_gradients())
                std::get<1>(derivatives_).resize(n_points);

            if (flags_handler_.fill_hessians())
                std::get<2>(derivatives_).resize(n_points);

            set_initialized(true);
        }

        void print_info(LogStream &out) const
        {
            flags_handler_.print_info(out);
            values_.print_info(out);
            std::get<1>(derivatives_).print_info(out);
            std::get<2>(derivatives_).print_info(out);
        }

        ValueVector<Point> points_;
        ValueVector<Value> values_;
        std::tuple<ValueVector<Derivative<0>>,
            ValueVector<Derivative<1>>,
            ValueVector<Derivative<2>>> derivatives_;
        FunctionFlags flags_handler_;
    };

    std::shared_ptr<Cache> elem_cache_;
public:
    using CacheType = Cache;
private:
    template <typename Accessor> friend class GridForwardIterator;
    friend class NewFunction<dim, codim, range, rank>;
};

IGA_NAMESPACE_CLOSE

#endif
