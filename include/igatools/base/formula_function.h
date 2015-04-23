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

#ifndef FORMULA_FUNCTIONS_H
#define FORMULA_FUNCTIONS_H

#include <igatools/base/function.h>
#include <igatools/base/value_types.h>

IGA_NAMESPACE_OPEN

/**
 *
 */
template<int dim, int codim=0, int range = 1, int rank = 1>
class FormulaFunction :
    public Function<dim, codim, range, rank>,
    public std::enable_shared_from_this<FormulaFunction<dim,codim,range,rank>>
{
private:
    using parent_t = Function<dim, codim, range, rank>;
    using self_t = FormulaFunction<dim, codim, range, rank>;
protected:
    using typename parent_t::GridType;
public:
    using typename parent_t::topology_variant;
    using typename parent_t::eval_pts_variant;
    using typename parent_t::Point;
    using typename parent_t::Value;
    using typename parent_t::Gradient;
    using typename parent_t::ElementIterator;
    using typename parent_t::ElementAccessor;
    using parent_t::space_dim;

    template <int order>
    using Derivative = typename parent_t::template Derivative<order>;

    using Map = MapFunction<dim, space_dim>;

    FormulaFunction(std::shared_ptr<GridType> grid, std::shared_ptr<Map> map);

    FormulaFunction(const self_t &func)
        :
        parent_t::Function(func),
        mapping_(func.mapping_->clone()),
        map_elem_(func.mapping_->begin())
    {}

    void reset(const ValueFlags &flag, const eval_pts_variant &quad) override
    {
        parent_t::reset(flag, quad);
        mapping_->reset(ValueFlags::value|ValueFlags::point, quad);
    }

    void init_cache(ElementAccessor &elem, const topology_variant &k) override
    {
        parent_t::init_cache(elem, k);
        mapping_->init_cache(map_elem_, k);
    }

    void fill_cache(ElementAccessor &elem, const topology_variant &k, const int j) override;

protected:
    std::shared_ptr<const parent_t> shared_from_derived() const override final
    {
        return this->shared_from_this();
    }

private:

    virtual void evaluate_0(const ValueVector<Point> &points,
                            ValueVector<Value> &values) const = 0;

    virtual void evaluate_1(const ValueVector<Point> &points,
                            ValueVector<Derivative<1>> &values) const = 0;

    virtual void evaluate_2(const ValueVector<Point> &points,
                            ValueVector<Derivative<2>> &values) const = 0;

private:
    std::shared_ptr<Map> mapping_;
    typename Map::ElementIterator map_elem_;

    struct FillCacheDispatcher : boost::static_visitor<void>
    {
        template<int sub_elem_dim>
        void operator()(const Topology<sub_elem_dim> &sub_elem)
        {
            auto &local_cache = function->get_cache(*elem);
            auto &cache = local_cache->template get_sub_elem_cache<sub_elem_dim>(j);
            auto &flags = cache.flags_handler_;

            if (!flags.fill_none())
            {
                auto &cache_pts = cache.template get_data<_Point>();
                cache_pts = function->map_elem_->template get_values<_Value, sub_elem_dim>(j);
                if (flags.template fill<_Value>())
                {
                    function->evaluate_0(cache_pts, cache.template get_data<_Value>());
                    flags.template set_filled<_Value>(true);
                }
                if (flags.template fill<_Gradient>())
                {
                    function->evaluate_1(cache_pts, cache.template get_data<_Gradient>());
                    flags.template set_filled<_Gradient>(true);
                }
                if (flags.template fill<_Hessian>())
                {
                    function->evaluate_2(cache_pts, cache.template get_data<_Hessian>());
                    flags.template set_filled<_Hessian>(true);
                }
                if (flags.template fill<_Divergence>())
                    Assert(false,ExcNotImplemented());

            }

            cache.set_filled(true);
        }

        int j;
        self_t *function;
        ElementAccessor *elem;
//        std::array<FunctionFlags, dim + 1> *flags_;
    };

    FillCacheDispatcher fill_cache_impl;
    friend class FillCacheDispatcher;
};

IGA_NAMESPACE_CLOSE

#endif
