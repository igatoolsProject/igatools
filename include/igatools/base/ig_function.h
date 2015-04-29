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

#ifndef __IG_FUNCTION_H
#define __IG_FUNCTION_H

#include <igatools/base/value_types.h>
#include <igatools/base/function.h>
#include <igatools/basis_functions/spline_space.h>
#include <igatools/linear_algebra/epetra_vector.h>

IGA_NAMESPACE_OPEN

using IgCoefficients = EpetraTools::Vector;

template<class Space>
class IgFunction :
    public Function<Space::dim, Space::codim, Space::range, Space::rank>,
    public std::enable_shared_from_this<IgFunction<Space>>
{
public:
    static const int dim = Space::dim;
    static const int codim = Space::codim;
    static const int range = Space::range;
    static const int rank = Space::rank;

    using CoeffType = IgCoefficients;


private:
    using base_t = Function<dim, codim, range, rank>;
    using parent_t = Function<dim, codim, range, rank>;
    using self_t = IgFunction<Space>;

    std::shared_ptr<const base_t> shared_from_derived() const override final
    {
        return this->shared_from_this();
    }

public:
    //TODO (pauletti, Mar 23, 2015): should we make this private?
    IgFunction(std::shared_ptr<const Space> space,
               const CoeffType &coeff,
               const std::string &property = DofProperties::active);

    IgFunction(const self_t &);

    virtual ~IgFunction() = default;

    using typename parent_t::topology_variant;
    using typename parent_t::eval_pts_variant;
    using typename parent_t::Point;
    using typename parent_t::Value;
    using typename parent_t::Gradient;
    using typename parent_t::ElementIterator;
    using typename parent_t::ElementAccessor;
    template <int order>
    using Derivative = typename parent_t::template Derivative<order>;



public:
    static std::shared_ptr<self_t>
    create(std::shared_ptr<const Space> space, const CoeffType &coeff,
           const std::string &property = DofProperties::active);


    std::shared_ptr<base_t> clone() const override final
    {
        return std::make_shared<self_t>(self_t(*this));
    }


    void reset(const ValueFlags &flag, const eval_pts_variant &eval_pts) override;

    void reset_selected_elements(const ValueFlags &flag,
                                 const eval_pts_variant &eval_pts,
                                 const vector<Index> &elements_flat_id);

    void init_cache(ElementAccessor &elem, const topology_variant &k) override;

    void fill_cache(ElementAccessor &elem, const topology_variant &k, const int j) override;

    std::shared_ptr<const Space> get_ig_space() const;

    const CoeffType &get_coefficients() const;

    const std::string &get_property() const
    {
        return property_;
    }

    self_t &operator +=(const self_t &fun);

    void print_info(LogStream &out) const;

private:
    std::shared_ptr<const Space> space_;

    CoeffType coeff_;

    const std::string property_;

    typename Space::ElementIterator elem_;

    std::shared_ptr<typename Space::ElementHandler> space_filler_;

private:
    struct ResetDispatcher : boost::static_visitor<void>
    {
        template<int sub_elem_dim>
        void operator()(const Quadrature<sub_elem_dim> &quad)
        {
            (*flags_)[sub_elem_dim] = flag;
            Assert(space_handler_ != nullptr, ExcNullPtr());
            space_handler_->reset_selected_elements(flag, quad,*elements_flat_id_);
        }

        ValueFlags flag;
        typename Space::ElementHandler *space_handler_;
        std::array<ValueFlags, dim + 1> *flags_;

        /**
         * Elements to reset.
         */
        const vector<Index> *elements_flat_id_;
    };


    struct InitCacheDispatcher : boost::static_visitor<void>
    {
        template<int sub_elem_dim>
        void operator()(const Topology<sub_elem_dim> &sub_elem)
        {
            Assert(space_handler_ != nullptr, ExcNullPtr());
            space_handler_->template init_cache<sub_elem_dim>(*space_elem);
        }

        typename Space::ElementHandler  *space_handler_;
        typename Space::ElementAccessor *space_elem;
    };


    struct FillCacheDispatcher : boost::static_visitor<void>
    {
        template<int sub_elem_dim>
        void operator()(const Topology<sub_elem_dim> &sub_elem)
        {
            Assert(space_handler_ != nullptr, ExcNullPtr());
            space_handler_->template fill_cache<sub_elem_dim>(*space_elem,j);

            auto &local_cache = function->get_cache(*func_elem);
            auto &cache = local_cache->template get_sub_elem_cache<sub_elem_dim>(j);

#if 0
            boost::fusion::for_each(cache.get_values(),
                                    [&](auto & type_and_value) -> void
            {
                using ValueType_ValueContainer = typename std::remove_reference<decltype(type_and_value)>::type;
                using ValueType = typename ValueType_ValueContainer::first_type;
                auto &value = type_and_value.second;

                if (value.status_fill())
                {
                    value = space_elem->template linear_combination<ValueType,sub_elem_dim>(*loc_coeff,j, *property);
                    value.set_status_filled(true);
                }
            } // end lambda function
                                                           );
#endif

            //TODO (martinelli Mar 27,2015): bad style. Use the ValueType mechanism in order to avoid the if-switch
            if (cache.template status_fill<_Value>())
            {
                cache.template get_data<_Value>() =
                    space_elem->template linear_combination<_Value,sub_elem_dim>(*loc_coeff,j, *property);
                cache.template set_status_filled<_Value>(true);
            }
            if (cache.template status_fill<_Gradient>())
            {
                cache.template get_data<_Gradient>() =
                    space_elem->template linear_combination<_Gradient,sub_elem_dim>(*loc_coeff,j, *property);
                cache.template set_status_filled<_Gradient>(true);
            }
            if (cache.template status_fill<_Hessian>())
            {
                cache.template get_data<_Hessian>() =
                    space_elem->template linear_combination<_Hessian,sub_elem_dim>(*loc_coeff,j, *property);
                cache.template set_status_filled<_Hessian>(true);
            }
            if (cache.template status_fill<_Divergence>())
            {
                cache.template get_data<_Divergence>() =
                    space_elem->template linear_combination<_Divergence,sub_elem_dim>(*loc_coeff,j, *property);
                cache.template set_status_filled<_Divergence>(true);
            }
//#endif
            cache.set_filled(true);
        }

        int j;
        self_t *function;
        typename Space::ElementHandler *space_handler_;
        ElementAccessor *func_elem;
        typename Space::ElementAccessor *space_elem;
        vector<Real> *loc_coeff;
        std::string const *property;
    };

    ResetDispatcher reset_impl;
    InitCacheDispatcher init_cache_impl;
    FillCacheDispatcher fill_cache_impl;

#ifdef REFINE
    void create_connection_for_insert_knots(std::shared_ptr<self_t> ig_function);

    void rebuild_after_insert_knots(
        const SafeSTLArray<vector<Real>,dim> &knots_to_insert,
        const CartesianGrid<dim> &old_grid);
#endif

};

IGA_NAMESPACE_CLOSE

#endif
