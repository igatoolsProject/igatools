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

#ifndef NEW_FUNCTIONS_H
#define NEW_FUNCTIONS_H

#include <igatools/base/config.h>
#include <igatools/base/tensor.h>
#include <igatools/utils/value_vector.h>
#include <igatools/geometry/grid_forward_iterator.h>
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
template<int dim, int codim = 0, int range = 1, int rank = 1>
class NewFunction : public GridElementHandler<dim>
{
private:
    using base_t = NewFunction<dim, codim, range, rank>;
    using self_t = NewFunction<dim, codim, range, rank>;
    using parent_t = GridElementHandler<dim>;

public:
    using typename parent_t::GridType;

public:
    static const int l= iga::max(0, dim-num_sub_elem);
    using v1 = typename seq<Quadrature, l, dim>::type;
    using variant_1 = typename boost::make_variant_over<v1>::type;

    using v2 = typename seq<Int, l, dim>::type;
    using variant_2 = typename boost::make_variant_over<v2>::type;

public:
    using ElementAccessor = FunctionElement<dim, codim, range, rank>;
    using ElementIterator = GridForwardIterator<ElementAccessor>;

    static const int space_dim = dim + codim;

    /** Types for the input/output evaluation arguments */
    ///@{
    /**
     * Type for the input argument of the function.
     */
    using Point = Points<space_dim>;

    /**
     * Type for the return of the function.
     */
    using Value = Values<space_dim, range, rank>;

    /**
     * Type for the derivative of the function.
     */
    template <int order>
    using Derivative = Derivatives<space_dim, range, rank, order>;

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
    using Div = Values<space_dim, range, rank-1>;
    ///@}

    /** @name Constructors and destructor. */
    ///@{
    /** Constructor */
    NewFunction(std::shared_ptr<GridType> grid);

    /** Destructor */
    virtual ~NewFunction() = default;
    ///@}

    NewFunction(const self_t &) = default;

    virtual std::shared_ptr<base_t> clone() const
    {
        Assert(false, ExcNotImplemented());
        return std::make_shared<self_t>(self_t(*this));
    }

    virtual void reset(const NewValueFlags &flag, const variant_1 &quad)
    {
        reset_impl.flag = flag;
        reset_impl.grid_handler = this;
        reset_impl.flags_ = &flags_;
        boost::apply_visitor(reset_impl, quad);
    }


    virtual void init_cache(ElementAccessor &elem, const variant_2 &k)
    {
        init_cache_impl.grid_handler = this;
        init_cache_impl.elem = &elem;
        init_cache_impl.flags_ = &flags_;
        init_cache_impl.quad_ = &(this->quad_);
        boost::apply_visitor(init_cache_impl, k);
    }


    virtual void fill_cache(ElementAccessor &elem, const int j, const variant_2 &k)
    {
        fill_cache_impl.j = j;
        fill_cache_impl.grid_handler = this;
        fill_cache_impl.elem = &elem;
        boost::apply_visitor(fill_cache_impl, k);
    }

    void fill_cache(ElementIterator &elem, const int j, const variant_2 &k)
    {
        fill_cache(elem.get_accessor(), j, k);
    }

    void init_cache(ElementIterator &elem, const variant_2 &k)
    {
        init_cache(elem.get_accessor(), k);
    }

    auto begin()  const -> ElementIterator
    {
        return ElementIterator(this->get_grid(), 0);
    }

    auto end() const -> ElementIterator
    {
        return ElementIterator(this->get_grid(),
                               IteratorState::pass_the_end);
    }

private:
    struct ResetDispatcher : boost::static_visitor<void>
    {
        template<class T>
        void operator()(const T &quad)
        {
            (*flags_)[T::dim] = flag;
            grid_handler->template reset<T::dim>(flag, quad);
        }

        NewValueFlags flag;
        parent_t *grid_handler;
        std::array<FunctionFlags, dim + 1> *flags_;
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

            for (auto &s_id: UnitElement<dim>::template elems_ids<T::k>())
            {
                auto &s_cache = cache->template get_value_cache<T::k>(s_id);
                auto &quad = std::get<T::k>(*quad_);
                s_cache.resize((*flags_)[T::k], quad.get_num_points());
            }
        }

        parent_t *grid_handler;
        ElementAccessor *elem;
        std::array<FunctionFlags, dim + 1> *flags_;
        QuadList<dim> *quad_;
    };

    ResetDispatcher reset_impl;
    FillCacheDispatcher fill_cache_impl;
    InitCacheDispatcher init_cache_impl;



//    virtual void init_elem(ElementAccessor &elem) = 0;
//
//    virtual void fill_elem(ElementAccessor &elem) = 0;
//
//    virtual void init_elem(ElementIterator &elem)
//    {
//        this->init_elem(elem.get_accessor());
//    }
//
//    virtual void fill_elem(ElementIterator &elem)
//    {
//        this->fill_elem(elem.get_accessor());
//    }

//protected:
public:
    std::shared_ptr<typename ElementAccessor::CacheType>
    &get_cache(ElementAccessor &elem);

protected:
    std::array<FunctionFlags, dim + 1> flags_;
};


template<int dim, int space_dim>
using MapFunction = NewFunction<dim, 0, space_dim, 1>;

IGA_NAMESPACE_CLOSE

#endif
