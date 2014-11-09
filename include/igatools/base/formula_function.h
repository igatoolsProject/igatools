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

/**
 *
 */
template<int dim, int codim=0, int range = 1, int rank = 1>
class FormulaFunction : public NewFunction<dim, codim, range, rank>
{
private:
    using parent_t = NewFunction<dim, codim, range, rank>;
    using self_t = FormulaFunction<dim, codim, range, rank>;
protected:
    using typename parent_t::GridType;
public:
    using typename parent_t::variant_1;
    using typename parent_t::variant_2;
    using typename parent_t::Point;
    using typename parent_t::Value;
    using typename parent_t::Gradient;
    using typename parent_t::ElementIterator;
    using typename parent_t::ElementAccessor;
    using parent_t::space_dim;

    template <int order>
    using Derivative = typename parent_t::template Derivative<order>;

    FormulaFunction(std::shared_ptr<GridType> grid);

//    void reset(const NewValueFlags &flag, const variant_1& quad) override;
//
//    void init_cache(ElementAccessor &elem, const variant_2& k) override;

    void fill_cache(ElementAccessor &elem, const int j, const variant_2& k) override;

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



    struct FillCacheDispatcher : boost::static_visitor<void>
    {
        template<class T>
        void operator()(const T& quad)
        {
            auto &local_cache = function->get_cache(*elem);
            auto &cache = local_cache->template get_value_cache<T::k>(j);
            auto &flags = cache.flags_handler_;

            if (!flags.fill_none())
            {
                const auto points =
                        elem->CartesianGridElement<dim>::template get_points<T::k>(j);

                if (flags.fill_points())
                    function->parametrization(points, cache.points_);
                if (flags.fill_values())
                    function->evaluate_0(cache.points_, std::get<0>(cache.values_));
                if (flags.fill_gradients())
                    function->evaluate_1(cache.points_, std::get<1>(cache.values_));
                if (flags.fill_hessians())
                    function->evaluate_2(cache.points_, std::get<2>(cache.values_));
            }

            cache.set_filled(true);
        }

        int j;
        self_t *function;
        ElementAccessor *elem;
        std::array<FunctionFlags, dim + 1> *flags_;
    };

    FillCacheDispatcher fill_cache_impl;
    friend class FillCacheDispatcher;
};

IGA_NAMESPACE_CLOSE

#endif
