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

#ifndef FUNCTION_ELEMENT_H
#define FUNCTION_ELEMENT_H

#include <igatools/base/function.h>
#include <igatools/geometry/cartesian_grid_element.h>
#include <igatools/base/value_types.h>

IGA_NAMESPACE_OPEN

template<int dim, int codim, int range = 1, int rank = 1>
class FunctionElement : public CartesianGridElement<dim>
{
public:
    using Func = Function<dim, codim, range, rank>;
    using Point = typename Func::Point;
    using Value = typename Func::Value;
    using Gradient = typename Func::Gradient;
    using Hessian  = typename Func::Hessian;
//    using ContainerType = const CartesianGrid<dim>;
    using ContainerType = const Func;

private:
    template <int order>
    using Derivative = typename Func::template Derivative<order>;

public:

    /** @name Constructors */
    ///@{
    /**
     * Default constructor.
     */
    FunctionElement() = delete;

    /**
     * Construct an accessor pointing to the element with
     * flat index @p elem_index of the CartesianGrid @p grid.
     */
    FunctionElement(const std::shared_ptr<const Func> func,
                    const Index elem_index)
        :
        CartesianGridElement<dim>(func->get_grid(),elem_index),
        func_(std::const_pointer_cast<Func>(func))
    {
        Assert(func_ != nullptr ,ExcNullPtr());
    }

    /**
     * Copy constructor.
     * It can be used with different copy policies
     * (i.e. deep copy or shallow copy).
     * The default behaviour (i.e. using the proper interface of a
     * classic copy constructor)
     * uses the deep copy.
     */
    FunctionElement(const FunctionElement<dim,codim,range,rank> &elem,
                    const CopyPolicy &copy_policy = CopyPolicy::deep);

    /**
     * Move constructor.
     */
    FunctionElement(FunctionElement<dim,codim,range,rank> &&elem) = default;

    /**
     * Destructor.
     */
    ~FunctionElement() = default;
    ///@}

    /**
     * @name Functions for performing different kind of copy.
     */
    ///@{
    /**
     * Performs a deep copy of the input @p element,
     * i.e. a new local cache is built using the copy constructor on the local cache of @p element.
     *
     * @note In DEBUG mode, an assertion will be raised if the input local cache is not allocated.
     */
    void deep_copy_from(const FunctionElement<dim,codim,range,rank> &element)
    {
        Assert(false,ExcNotImplemented());
    }

    /**
     * Performs a shallow copy of the input @p element. The current object will contain a pointer to the
     * local cache used by the input @p element.
     */
    void shallow_copy_from(const FunctionElement<dim,codim,range,rank> &element)
    {
        Assert(false,ExcNotImplemented());
    }

    ///@}


    /** @name Assignment operators */
    ///@{
    /**
     * Copy assignment operator. Performs a <b>shallow copy</b> of the input @p element.
     *
     * @note Internally it uses the function shallow_copy_from().
     */
    FunctionElement<dim,codim,range,rank> &operator=(const FunctionElement<dim,codim,range,rank> &element)
    {
        shallow_copy_from(element);
        return *this;
    }

    /**
     * Move assignment operator.
     */
    FunctionElement<dim,codim,range,rank> &operator=(FunctionElement<dim,codim,range,rank> &&elem) = default;
    ///@}



    template<int order, int k>
    auto
    get_values(const int j) const
    {
        Assert(local_cache_ != nullptr,ExcNullPtr());
        const auto &cache = local_cache_->template get_value_cache<k>(j);
        Assert(cache.is_filled() == true, ExcCacheNotFilled());
        return std::get<order>(cache.values_);
    }



    /**
     * @name Methods for the for the evaluations of Functions's derivatives
     *  without the use of the cache.
     */
    ///@{
    /**
     * Returns a ValueTable with the <tt>deriv_order</tt>-th derivatives of all local basis function
     * at each point (in the unit domain) specified by the input argument <tt>points</tt>.
     * @note This function does not use the cache and therefore can be called any time without
     * needing to pre-call init_cache()/fill_cache().
     * @warning The evaluation <tt>points</tt> must belong to the unit hypercube
     * \f$ [0,1]^{\text{dim}} \f$ otherwise, in Debug mode, an assertion will be raised.
     */
    template <int deriv_order>
    ValueVector<
    Conditional< deriv_order==0,
                 Value,
                 Derivative<deriv_order> > >
                 evaluate_derivatives_at_points(const Quadrature<dim> &points)
    {
        ValueFlags flags;
        if (deriv_order == 0)
            flags = ValueFlags::value;
        else if (deriv_order == 1)
            flags = ValueFlags::gradient;
        else if (deriv_order == 2)
            flags = ValueFlags::hessian;
        else
        {
            Assert(false,ExcNotImplemented());
        }

        func_->reset_one_element(flags,points,this->get_flat_index());
        const auto topology = Int<dim>();
        func_->init_cache(*this,topology);
        func_->fill_cache(*this,topology,0);

        return this->template get_values<deriv_order,dim>(0);


//        Assert(false,ExcNotImplemented());
//        return dummy;
    }

    ValueVector<Value>
    evaluate_values_at_points(const Quadrature<dim> &points)
    {
        return this->template evaluate_derivatives_at_points<0>(points);
    }

    ValueVector<Derivative<1> >
    evaluate_gradients_at_points(const Quadrature<dim> &points)
    {
        return this->template evaluate_derivatives_at_points<1>(points);
    }

    ValueVector<Derivative<2> >
    evaluate_hessians_at_points(const Quadrature<dim> &points)
    {
        return this->template evaluate_derivatives_at_points<2>(points);
    }
    ///@}

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

            if (flags_handler_.fill<_Value>())
                std::get<0>(values_).resize(n_points);

            if (flags_handler_.fill<_Gradient>())
                std::get<1>(values_).resize(n_points);

            if (flags_handler_.fill<_Hessian>())
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
    std::shared_ptr<Func> func_;

    template <class Accessor> friend class CartesianGridIteratorBase;
    friend class Function<dim, codim, range, rank>;

    /**
     * Creates a new object performing a deep copy of the current object using the FunctionElement
     * copy constructor.
     */
    std::shared_ptr<FunctionElement<dim,codim,range,rank> > clone() const
    {
        auto elem = std::shared_ptr<FunctionElement<dim,codim,range,rank> >(
                        new FunctionElement(*this,CopyPolicy::deep));
        Assert(elem != nullptr, ExcNullPtr());
        return elem;
    }
};

IGA_NAMESPACE_CLOSE

#endif
