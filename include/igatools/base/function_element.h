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
    using ContainerType = const CartesianGrid<dim>;

private:
    template <int order>
    using Derivative = typename Func::template Derivative<order>;
public:
    using CartesianGridElement<dim>::CartesianGridElement;

    template<int order, int k>
    auto
    get_values(const int j) const
    {
        const auto &cache = local_cache_->template get_value_cache<k>(j);
        Assert(cache.is_filled() == true, ExcCacheNotFilled());
        return std::get<order>(cache.values_);
    }

#if 0
    ValueVector<Point> const &get_points() const;
    ValueVector<Value> const &get_values() const;

private:
    template<int order>
    auto const &get_derivative() const;
public:
    ValueVector<Gradient> const &get_gradients() const;

    ValueVector<Hessian> const &get_hessians() const;
#endif

private:
    class ValuesCache : public CacheStatus
    {
    public:
        void resize(const FunctionFlags &flags_handler,
                    const int n_points)
        {
            //TODO(pauletti, Oct 11, 2014): missing all necesary clears
            flags_handler_ = flags_handler;

            if (flags_handler_.fill_points())
                points_.resize(n_points);

            if (flags_handler_.fill_values())
                std::get<0>(values_).resize(n_points);

            if (flags_handler_.fill_gradients())
                std::get<1>(values_).resize(n_points);

            if (flags_handler_.fill_hessians())
                std::get<2>(values_).resize(n_points);

            set_initialized(true);
        }

        void print_info(LogStream &out) const
        {
            flags_handler_.print_info(out);
            std::get<0>(values_).print_info(out);
            std::get<1>(values_).print_info(out);
            std::get<2>(values_).print_info(out);
        }

        FunctionFlags flags_handler_;

        ValueVector<Point> points_;
        std::tuple<ValueVector<Value>,
            ValueVector<Derivative<1>>,
            ValueVector<Derivative<2>>> values_;

    };

    class LocalCache
    {
    public:
        LocalCache() = default;

        LocalCache(const LocalCache &in) = default;

        LocalCache(LocalCache &&in) = default;

        ~LocalCache() = default;


        LocalCache &operator=(const LocalCache &in) = delete;

        LocalCache &operator=(LocalCache &&in) = delete;

        void print_info(LogStream &out) const;

        template <int k>
        ValuesCache &
        get_value_cache(const int j)
        {
            return std::get<k>(values_)[j];
        }

        template <int k>
        const ValuesCache &
        get_value_cache(const int j) const
        {
            return std::get<k>(values_)[j];
        }

        CacheList<ValuesCache, dim> values_;
    };

    std::shared_ptr<LocalCache> local_cache_;

public:
    using CacheType = LocalCache;
private:
    template <typename Accessor> friend class GridForwardIterator;
    friend class NewFunction<dim, codim, range, rank>;
};

IGA_NAMESPACE_CLOSE

#endif
