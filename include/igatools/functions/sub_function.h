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

#ifndef SUB_FUNCTION_H_
#define SUB_FUNCTION_H_

#include <igatools/functions/function.h>
#include <igatools/functions/function_element.h>
#include <boost/variant/get.hpp>
//#include<igatools/../../source/geometry/grid_forward_iterator.cpp>

IGA_NAMESPACE_OPEN

/**
 *
 * @author pauletti 2014
 */
template<int sub_dim, int dim, int codim, int range, int rank>
class SubFunction :
    public Function<sub_dim, codim + (dim-sub_dim), range, rank>,
    public std::enable_shared_from_this<SubFunction<sub_dim,dim,codim,range,rank>>
{
private:
    using self_t = SubFunction<sub_dim, dim, codim, range, rank>;


public:
    using base_t  = Function<sub_dim, codim + (dim-sub_dim), range, rank>;
    using SupFunc = Function<dim, codim, range, rank>;

    using typename base_t::GridType;

    using typename base_t::topology_variant;
    using typename base_t::eval_pts_variant;
    using typename base_t::ElementAccessor;

    using SuperGrid = typename SupFunc::GridType;
    template <int j>
    using InterGridMap = typename SuperGrid::template InterGridMap<j>;

private:
    std::shared_ptr<const base_t> shared_from_derived() const override final
    {
        return this->shared_from_this();
    }

public:

    SubFunction(std::shared_ptr<GridType> grid,
                std::shared_ptr<const SupFunc> func,
                const int s_id,
                InterGridMap<sub_dim> &elem_map)
        :
        base_t(grid),
        sup_func_(func->clone()),
        s_id_(s_id),
        elem_map_(elem_map),
        sup_elem_(sup_func_->begin())
    {}

    SubFunction(const self_t &sub_f)
        :
        base_t(sub_f),
        sup_func_(sub_f.sup_func_->clone()),
        s_id_(sub_f.s_id_),
        elem_map_(sub_f.elem_map_),
        sup_elem_(sub_f.sup_func_->begin())
    {}


    static std::shared_ptr<base_t>
    create(std::shared_ptr<GridType> grid,
           std::shared_ptr<const SupFunc> func,
           const int s_id,
           InterGridMap<sub_dim> &elem_map)
    {
        return std::make_shared<self_t>(grid, func, s_id, elem_map);
    }


    std::shared_ptr<base_t> clone() const override
    {

        return std::make_shared<self_t>(*this);
    }

    void reset(const ValueFlags &flag, const eval_pts_variant &eval_pts) override
    {
        base_t::reset(flag, eval_pts);
        auto q = boost::get<Quadrature<sub_dim>>(eval_pts);
        sup_func_->reset(flag, q);

    }


    void init_cache(ElementAccessor &elem, const topology_variant &k1) override
    {
        base_t::init_cache(elem, k1);
        sup_func_->init_cache(sup_elem_, Topology<sub_dim>());
    }

    void fill_cache(ElementAccessor &elem, const topology_variant &k1, const int j) override
    {
        Assert(j==0, ExcNotImplemented());
        using ElementIt = typename CartesianGrid<sub_dim>::ElementIterator;
        ElementIt el_it(elem.get_grid(),elem.get_flat_index(),ElementProperties::none);

        sup_elem_->move_to(elem_map_[el_it]->get_flat_index());

        base_t::fill_cache(elem,k1,j);
        sup_func_->fill_cache(sup_elem_, Topology<sub_dim>(),s_id_);
        auto &local_cache = this->get_cache(elem);
        auto &cache = local_cache->template get_sub_elem_cache<sub_dim>(j);
//        auto &flags = cache.flags_handler_;

        if (cache.template status_fill<_Value>())
        {
            cache.template get_data<_Value>() = sup_elem_->template get_values<_Value, sub_dim>(s_id_);
            cache.template set_status_filled<_Value>(true);
        }
        if (cache.template status_fill<_Gradient>())
        {
            auto active = UnitElement<dim>::template get_elem<sub_dim>(s_id_).active_directions;
            auto DSupF  = sup_elem_->template get_values<_Gradient, sub_dim>(s_id_);
            auto &DSubF = cache.template get_data<_Gradient>();

            const auto n_points = DSupF.get_num_points();
            for (int pt = 0; pt < n_points; ++pt)
            {
                int j = 0;
                for (auto &i : active)
                {
                    DSubF[pt][j] = DSupF[pt][i];
                    ++j;
                }
            }
            cache.template set_status_filled<_Gradient>(true);
        }

        if (cache.template status_fill<_Hessian>())
        {
            Assert(false, ExcNotImplemented());
//      std::get<2>(cache.values_) = sup_elem_->template get_values<2, sub_dim>(j);
        }
        cache.set_filled(true);

    }

private:
    std::shared_ptr<SupFunc> sup_func_;
    const int s_id_;
    InterGridMap<sub_dim> &elem_map_;

    typename SupFunc::ElementIterator sup_elem_;

};



template<int sub_dim, int dim, int space_dim>
class SubMapFunction :
    public MapFunction<sub_dim, space_dim>,
    public std::enable_shared_from_this<SubMapFunction<sub_dim,dim,space_dim>>
{
private:
    using self_t = SubMapFunction<sub_dim, dim, space_dim>;
public:
    using base_t  = MapFunction<sub_dim, space_dim>;
    using SupFunc = MapFunction<dim, space_dim>;

    using typename base_t::GridType;

    using typename base_t::topology_variant;
    using typename base_t::eval_pts_variant;
    using typename base_t::ElementAccessor;

    using SuperGrid = typename SupFunc::GridType;
    template <int j>
    using InterGridMap = typename SuperGrid::template InterGridMap<j>;

private:
    std::shared_ptr<const base_t> shared_from_derived() const override final
    {
        return this->shared_from_this();
    }


public:

    SubMapFunction(std::shared_ptr<GridType> grid,
                   std::shared_ptr<const SupFunc> func,
                   const int s_id,
                   InterGridMap<sub_dim> &elem_map)
        :
        base_t(grid),
        sup_func_(func->clone()),
        s_id_(s_id),
        elem_map_(elem_map),
        sup_elem_(sup_func_->begin())
    {}

    SubMapFunction(std::shared_ptr<const SupFunc> func,
                   const int s_id)
        :
        base_t(func->get_grid()->template get_sub_grid<sub_dim>(s_id, elem_map_)),
        sup_func_(func->clone()),
        s_id_(s_id),
        sup_elem_(sup_func_->begin())
    {}

    SubMapFunction(const self_t &sub_f)
        :
        base_t(sub_f),
        sup_func_(sub_f.sup_func_->clone()),
        s_id_(sub_f.s_id_),
        elem_map_(sub_f.elem_map_),
        sup_elem_(sub_f.sup_func_->begin())
    {}

    std::shared_ptr<base_t> clone() const override
    {

        return std::make_shared<self_t>(self_t(*this));
    }


    static std::shared_ptr<base_t>
    create(std::shared_ptr<const SupFunc> func,
           const int s_id)
    {
        return std::make_shared<self_t>(func, s_id);
    }

    static std::shared_ptr<base_t>
    create(std::shared_ptr<GridType> grid,
           std::shared_ptr<const SupFunc> func,
           const int s_id,
           InterGridMap<sub_dim> &elem_map)
    {
        return std::make_shared<self_t>(grid, func, s_id, elem_map);
    }

    void reset(const ValueFlags &flag, const eval_pts_variant &eval_pts) override
    {
        base_t::reset(flag, eval_pts);
        auto q = boost::get<Quadrature<sub_dim>>(eval_pts);
        sup_func_->reset(flag, q);

    }


    void init_cache(ElementAccessor &elem, const topology_variant &k1) override
    {
        base_t::init_cache(elem, k1);
        sup_func_->init_cache(sup_elem_, Topology<sub_dim>());
    }

    void fill_cache(ElementAccessor &elem, const topology_variant &k1, const int j) override
    {
        Assert(j==0, ExcNotImplemented());
//        typename CartesianGrid<sub_dim>::ElementIterator el_it(elem);
        using ElementIt = typename CartesianGrid<sub_dim>::ElementIterator;
        ElementIt el_it(elem.get_grid()->create_element(elem.get_flat_index()),ElementProperties::none);

        sup_elem_->move_to(elem_map_[el_it]->get_flat_index());

        base_t::fill_cache(elem, k1, j);
        sup_func_->fill_cache(sup_elem_,Topology<sub_dim>(),s_id_);
        auto &local_cache = this->get_cache(elem);
        auto &cache = local_cache->template get_sub_elem_cache<sub_dim>(j);
//        auto &flags = cache.flags_handler_;

        if (cache.template status_fill<_Value>())
        {
            cache.template get_data<_Value>() = sup_elem_->template get_values<_Value, sub_dim>(s_id_);

            cache.template set_status_filled<_Value>(true);
        }
        if (cache.template status_fill<_Gradient>())
        {
            auto active = UnitElement<dim>::template get_elem<sub_dim>(s_id_).active_directions;
            auto DSupF  = sup_elem_->template get_values<_Gradient, sub_dim>(s_id_);
            auto &DSubF = cache.template get_data<_Gradient>();

            const auto n_points = DSupF.get_num_points();
            for (int pt = 0; pt<n_points; ++pt)
            {
                int j = 0;
                for (auto &i : active)
                {
                    DSubF[pt][j] = DSupF[pt][i];
                    ++j;
                }
            }

            cache.template set_status_filled<_Gradient>(true);
        }
        if (cache.template status_fill<_Hessian>())
        {
            Assert(false, ExcNotImplemented());
            //      std::get<2>(cache.values_) = sup_elem_->template get_values<2, sub_dim>(j);
        }

        cache.set_filled(true);

    }

private:
    std::shared_ptr<SupFunc> sup_func_;
    const int s_id_;
    InterGridMap<sub_dim> &elem_map_;

    typename SupFunc::ElementIterator sup_elem_;

};

IGA_NAMESPACE_CLOSE

#endif

