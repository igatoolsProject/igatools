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

#ifndef NEW_FUNCTIONS_H
#define NEW_FUNCTIONS_H

#include <igatools/base/config.h>
#include <igatools/base/tensor.h>
#include <igatools/utils/value_vector.h>
#include <igatools/geometry/grid_wrapper.h>
#include <igatools/geometry/cartesian_grid_iterator.h>
#include <igatools/geometry/grid_element_handler.h>
#include <igatools/base/quadrature.h>

#include <boost/variant.hpp>
#include <boost/mpl/vector.hpp>
IGA_NAMESPACE_OPEN

template<int k_>
struct Int
{
    static const int k = k_;
};

template<template<int> class Q, int start, std::size_t N>
struct seq;

template<template<int> class Q, int start>
struct seq<Q, start, start>
{
    using type = boost::mpl::vector<Q<start>>;
};

template<template<int> class Q, int start, std::size_t N>
struct seq
{
    using v1 = typename seq<Q, start, N-1>::type;
    using type = typename boost::mpl::push_back<v1, Q<N>>::type;
};

template <int, int, int, int> class FunctionElement;



/**
 * Function Class
 */
template<int dim_, int codim_ = 0, int range_ = 1, int rank_ = 1>
class Function
    : public GridElementHandler<dim_>
{
private:
    using base_t = Function<dim_, codim_, range_, rank_>;
    using self_t = Function<dim_, codim_, range_, rank_>;
    using parent_t = GridElementHandler<dim_>;


    virtual std::shared_ptr<const self_t> shared_from_derived() const = 0;

public:
    using Topology = UnitElement<dim_>;
    using GridType = CartesianGrid<dim_>;

public:
    static const int l = iga::max(0, dim_-num_sub_elem);

    using v2 = typename seq<Int, l, dim_>::type;
    using topology_variant = typename boost::make_variant_over<v2>::type;

    using v3 = typename seq<Quadrature, l, dim_>::type;
    using eval_pts_variant = typename boost::make_variant_over<v3>::type;

public:
    using ElementAccessor = FunctionElement<dim_, codim_, range_, rank_>;
    using ElementIterator = CartesianGridIterator<ElementAccessor>;

    static const int space_dim = dim_ + codim_;
    static const int dim       = dim_;
    static const int codim     = codim_;
    static const int range     = range_;
    static const int rank      = rank_;

    /** Types for the input/output evaluation arguments */
    ///@{
    /**
     * Type for the input argument of the function.
     */
    using Point = Points<space_dim>;

    /**
     * Type for the return of the function.
     */
    using Value = Values<space_dim, range_, rank_>;

    /**
     * Type for the derivative of the function.
     */
    template <int order>
    using Derivative = Derivatives<space_dim, range_, rank_, order>;

    /**
     * Type for the gradient of the function.
     */
    using Gradient = Derivative<1>;

    /**
     * Type for the hessian of the function.
     */
    using Hessian = Derivative<2>;

    /**
     * Type for the divergence of function.
     */
    using Div = Values<space_dim, range_, rank_-1>;
    ///@}

    /** @name Constructors and destructor. */
    ///@{
protected:
    /** Constructor */
    Function(std::shared_ptr<GridType> grid);


    Function(const self_t &) = default;

public:
    /** Destructor */
    virtual ~Function() = default;
    ///@}


    virtual std::shared_ptr<base_t> clone() const = 0;

#if 0
    virtual std::shared_ptr<base_t> clone() const
    {
        Assert(false, ExcNotImplemented());
        return std::make_shared<self_t>(self_t(*this));
    }
#endif

    virtual void reset(const ValueFlags &flag, const eval_pts_variant &quad)
    {
        reset_impl.flag = flag;
        reset_impl.grid_handler = this;
        reset_impl.flags_ = &flags_;
        boost::apply_visitor(reset_impl, quad);
    }

    void reset_one_element(
        const ValueFlags &flag,
        const eval_pts_variant &eval_pts,
        const Index elem_id)
    {
#ifndef NDEBUG
        const auto grid = this->get_grid();
//        Assert(grid->is_element_active(elem_id),
//               ExcMessage("The element " + std::to_string(elem_id) + " is not active."));
#endif

        this->reset(flag,eval_pts);
    }

    virtual void init_cache(ElementAccessor &elem, const topology_variant &k)
    {
        init_cache_impl.grid_handler = this;
        init_cache_impl.elem = &elem;
        init_cache_impl.flags_ = &flags_;
        init_cache_impl.quad_ = &(this->quad_);
        boost::apply_visitor(init_cache_impl, k);
    }

    void init_cache(ElementIterator &elem, const topology_variant &k)
    {
        init_cache(*elem, k);
    }

    template <int k>
    void init_cache(ElementAccessor &elem)
    {
        const auto topology = Int<k>();
        this->init_cache(elem, topology);
    }

    void init_element_cache(ElementAccessor &elem)
    {
        this->init_cache(elem, Int<dim_>());
    }

    template <int k>
    void init_cache(ElementIterator &elem)
    {
        this->template init_cache<k>(*elem);
    }

    void init_element_cache(ElementIterator &elem)
    {
        this->init_cache(*elem, Int<dim_>());
    }

    virtual void fill_cache(ElementAccessor &elem, const topology_variant &k,const int j)
    {
        fill_cache_impl.j = j;
        fill_cache_impl.grid_handler = this;
        fill_cache_impl.elem = &elem;
        boost::apply_visitor(fill_cache_impl, k);
    }

    void fill_cache(ElementIterator &elem, const topology_variant &k, const int j)
    {
        fill_cache(*elem, k, j);
    }

    template <int k>
    void fill_cache(ElementAccessor &elem, const int j)
    {
        const auto topology = Int<k>();
        this->fill_cache(elem, topology,j);
    }

    void fill_element_cache(ElementAccessor &elem)
    {
        this->fill_cache(elem, Int<dim_>(),0);
    }

    template <int k>
    void fill_cache(ElementIterator &elem, const int j)
    {
        this->template fill_cache<k>(*elem,j);
    }

    void fill_element_cache(ElementIterator &elem)
    {
        this->fill_cache(*elem, Int<dim_>(),0);
    }

    std::shared_ptr<ElementAccessor> create_element(const Index flat_index) const
    {
        auto elem = std::shared_ptr<ElementAccessor>(
                        new ElementAccessor(this->shared_from_derived(),flat_index));
        Assert(elem != nullptr,ExcNullPtr());

        return elem;
    }

    auto begin()  const -> ElementIterator
    {
        return ElementIterator(this->create_element(0),ElementProperties::none);
    }

    auto end() const -> ElementIterator
    {
        return ElementIterator(this->create_element(IteratorState::pass_the_end),ElementProperties::none);
    }



    virtual void print_info(LogStream &out) const
    {
        using std::to_string;
        out.begin_item("Function<" + to_string(dim_) + "," +
                       to_string(codim_) + "," +
                       to_string(range_) + "," +
                       to_string(rank_) + ">");
        parent_t::print_info(out);
        out.end_item();
    }

private:


    struct ResetDispatcher : boost::static_visitor<void>
    {
        template<class T>
        void operator()(const T &quad)
        {
            (*flags_)[T::dim] = flag;

            grid_handler->template reset<T::dim>(FunctionFlags::to_grid_flags(flag), quad);
        }

        ValueFlags flag;
        parent_t *grid_handler;
        std::array<FunctionFlags, dim_ + 1> *flags_;
    };

    struct FillCacheDispatcher : boost::static_visitor<void>
    {
        template<class T>
        void operator()(const T &quad)
        {
            grid_handler->template fill_cache<T::k>(*elem, j);
        }

        int j;
        parent_t *grid_handler;
        ElementAccessor *elem;
    };

    struct InitCacheDispatcher : boost::static_visitor<void>
    {
        template<class T>
        void operator()(const T &quad)
        {
            grid_handler->template init_cache<T::k>(*elem);

            auto &cache = elem->local_cache_;
            if (cache == nullptr)
            {
                using Cache = typename ElementAccessor::LocalCache;
                cache = std::shared_ptr<Cache>(new Cache);
            }

            for (auto &s_id: UnitElement<dim_>::template elems_ids<T::k>())
            {
                auto &s_cache = cache->template get_value_cache<T::k>(s_id);
                auto &quad = std::get<T::k>(*quad_);
                s_cache.resize((*flags_)[T::k], quad.get_num_points());
            }
        }

        parent_t *grid_handler;
        ElementAccessor *elem;
        std::array<FunctionFlags, dim_ + 1> *flags_;
        EvalPtsList<dim_> *quad_;
    };

    ResetDispatcher reset_impl;
    FillCacheDispatcher fill_cache_impl;
    InitCacheDispatcher init_cache_impl;


public:
    std::shared_ptr<typename ElementAccessor::CacheType>
    &get_cache(ElementAccessor &elem);

protected:
    std::array<FunctionFlags, dim_ + 1> flags_;


    /**
     * This member is used to handle the knots-refinement.
     */
    GridWrapper<CartesianGrid<dim_> > functions_knots_refinement_;

public:

    /**
     * Perform the h-refinement of the grid in all the directions.
     * Each interval in the unrefined grid is uniformly divided in @p n_subdivisions sub-intervals.
     */
    void refine_h(const Size n_subdivisions=2)
    {
        functions_knots_refinement_.refine_h(n_subdivisions);
    }
};


template<int dim, int space_dim>
using MapFunction = Function<dim, 0, space_dim, 1>;

IGA_NAMESPACE_CLOSE

#endif
