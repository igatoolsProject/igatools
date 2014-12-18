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

#ifndef BSPLINE_ELEMENT_HANDLER_H_
#define BSPLINE_ELEMENT_HANDLER_H_

#include <igatools/base/config.h>
#include <igatools/base/cache_status.h>
#include <igatools/base/flags_handler.h>
#include <igatools/base/quadrature.h>

//TODO(pauletti, Sep 9, 2014): should we instantiate the cartesian product instead
#include <igatools/utils/cartesian_product_array-template.h>

#include <igatools/utils/value_table.h>
#include <igatools/geometry/grid_element_handler.h>
#include <igatools/basis_functions/bspline_space.h>
#include <igatools/basis_functions/bspline_element_scalar_evaluator.h>


IGA_NAMESPACE_OPEN

template<int dim_, int range_ = 1, int rank_ = 1>
class ReferenceElementHandler
//      : protected GridElementHandler<dim_>
{
public:
//    using base_t = GridElementHandler<dim_>;
    using Space = ReferenceSpace<dim_,range_,rank_>;
    using ElementIterator = typename Space::ElementIterator;
    using ElementAccessor = typename Space::ElementAccessor;

//    using base_t::get_num_points;

    static const int l = iga::max(0, dim_-num_sub_elem);
    using v1 = typename seq<Quadrature, l, dim_>::type;
    using quadrature_variant = typename boost::make_variant_over<v1>::type;

    using v2 = typename seq<Int, l, dim_>::type;
    using topology_variant = typename boost::make_variant_over<v2>::type;

    //Allocates and fill the (global) cache
    ReferenceElementHandler(std::shared_ptr<const Space> space)
        :
//        base_t(space->get_grid()),
        grid_handler_(space->get_grid()),
        space_(space)
    {
        Assert(space != nullptr, ExcNullPtr());
    };

    virtual ~ReferenceElementHandler() = default;

    ReferenceElementHandler(const ReferenceElementHandler<dim_,range_,rank_> &elem_handler) = delete;
    ReferenceElementHandler(ReferenceElementHandler<dim_,range_,rank_> &&elem_handler) = delete;

    virtual void reset(const ValueFlags &flag, const quadrature_variant &quad) = 0;

    virtual void init_cache(ElementAccessor &elem, const topology_variant &topology) = 0;

    void init_cache(ElementIterator &elem, const topology_variant &topology)
    {
        init_cache(*elem,topology);
    }

    //Allocates the ElementIterator element_cache
    void init_element_cache(ElementIterator &elem)
    {
        init_cache(*elem,Int<dim_>());
    }

    virtual void fill_cache(ElementAccessor &elem, const topology_variant &topology, const int j) = 0;

    void fill_cache(ElementIterator &elem, const topology_variant &topology, const int j)
    {
        fill_cache(*elem,topology,j);
    }

    //Fill the ElementIterator element_cache
    void fill_element_cache(ElementIterator &elem)
    {
        fill_cache(*elem,Int<dim_>(),0);
    }


    template <int k>
    void fill_cache(ElementAccessor &elem, const int j)
    {
        Assert(false,ExcNotImplemented());
    }

    virtual void print_info(LogStream &out) const = 0;

    template <int k = dim_>
    Size get_num_points() const
    {
        return grid_handler_.template get_num_points<k>();
    }

protected:
    GridElementHandler<dim_> grid_handler_;

private:
    std::shared_ptr<const Space> space_;

public:
    std::shared_ptr<const Space> get_space() const
    {
        Assert(space_ != nullptr,ExcNullPtr());
        return space_;
    }
};


/**
 * Global BSplineSpace uniform quadrature
 * computational optimization cache, storing the interval length
 * in each direction
 */
template<int dim_, int range_ = 1, int rank_ = 1>
class BSplineElementHandler : public ReferenceElementHandler<dim_,range_,rank_>
{
    using base_t = ReferenceElementHandler<dim_,range_,rank_>;
    using self_t = BSplineElementHandler<dim_,range_,rank_>;
    using Space = BSplineSpace<dim_,range_,rank_>;
    static const Size n_components =  Space::n_components;


    template<class T>
    using ComponentContainer = typename Space::template ComponentContainer<T>;

    template<class T>
    using ComponentDirectionTable = ComponentContainer<CartesianProductArray<T,dim_>>;

    template<class T>
    using ComponentDirectionContainer = ComponentContainer<std::array<T,dim_>>;

    using SpaceDimensionTable = typename Space::SpaceDimensionTable;

    template <int order>
    using Derivative = typename Space::template Derivative<order>;

    using Value = typename Space::Value;

protected:

    using BaseSpace = ReferenceSpace<dim_,range_,rank_>;
    using RefElementIterator = typename BaseSpace::ElementIterator;
    using RefElementAccessor = typename BaseSpace::ElementAccessor;

    using ElementIterator = typename Space::ElementIterator;
    using ElementAccessor = typename Space::ElementAccessor;



public:
    static const int dim = dim_;

protected:
    //Allocates and fill the (global) cache
    BSplineElementHandler(std::shared_ptr<const Space> space);

public:
    virtual ~BSplineElementHandler() = default;

    static std::shared_ptr<self_t> create(std::shared_ptr<const Space> space)
    {
        return std::shared_ptr<self_t>(new self_t(space));
    }

    using quadrature_variant = typename base_t::quadrature_variant;
    using topology_variant = typename base_t::topology_variant;

    virtual void reset(const ValueFlags &flag, const quadrature_variant &quad) override final;

    virtual void init_cache(RefElementAccessor &elem, const topology_variant &topology) override final;

    void init_cache(ElementIterator &elem, const topology_variant &topology)
    {
        init_cache(*elem,topology);
    }

    //Allocates the ElementIterator element_cache
    void init_element_cache(ElementIterator &elem)
    {
        init_cache(*elem,Int<dim_>());
    }


    virtual void fill_cache(RefElementAccessor &elem, const topology_variant &topology, const int j) override final;


    void fill_cache(ElementIterator &elem, const topology_variant &topology, const int j)
    {
        fill_cache(*elem,topology,j);
    }

    //Fill the ElementIterator element_cache
    void fill_element_cache(ElementIterator &elem)
    {
        fill_cache(*elem,Int<dim_>(),0);
    }

public:

    virtual void print_info(LogStream &out) const override final ;

    const GridElementHandler<dim_> &get_grid_handler() const
    {
        return this->grid_handler_;
    }

private:
//    std::shared_ptr<const Space> space_;

//    ComponentContainer<Size> comp_offset_;

    std::array<FunctionFlags, dim + 1> flags_;

    template <class T>
    using DirectionTable = CartesianProductArray<T, dim_>;
    using BasisValues = ComponentContainer<BasisValues1d>;
    /**
     * B-splines values and derivatives at quadrature points.
     * The values are stored in the un tensor product way.
     *
     * splines1d_[dir][interval][comp][order][function][point]
     */
    class GlobalCache : public DirectionTable<BasisValues>
    {
    public:
        using DirectionTable<BasisValues>::DirectionTable;
        // TODO (pauletti, Sep 23, 2014): document and split definition
        auto get_element_values(const TensorIndex<dim> &id) const
        {
            ComponentContainer<TensorProductFunctionEvaluator<dim> >
            result((this->entry(0,0)).get_comp_map());
            for (auto c : result.get_active_components_id())
            {
                for (int i = 0; i < dim; ++i)
                    result[c][i] =
                        BasisValues1dConstView((this->entry(i, id[i]))[c]);
                result[c].update_size();

            }
            return result;
        }
    };

    CacheList<GlobalCache, dim> splines1d_;


    struct ResetDispatcher : boost::static_visitor<void>
    {
        template<class T>
        void operator()(const T &quad);

        GridElementHandler<dim_> *grid_handler_;
        ValueFlags flag_;
        std::array<FunctionFlags, dim + 1> *flags_;
        CacheList<GlobalCache, dim> *splines1d_;
        const Space *space_;
    };

    ResetDispatcher reset_impl_;


    struct InitCacheDispatcher : boost::static_visitor<void>
    {
        template<class T>
        void operator()(const T &quad);

        GridElementHandler<dim_> *grid_handler_;
        ReferenceElement<dim_,range_,rank_> *elem_;
        std::array<FunctionFlags, dim + 1> *flags_;
    };

    InitCacheDispatcher init_cache_impl_;


    struct FillCacheDispatcher : boost::static_visitor<void>
    {
        template<class T>
        void operator()(const T &quad);

        /**
         * Computes the values (i.e. the 0-th order derivative) of the non-zero B-spline basis
         * functions over the current element,
         *   at the evaluation points pre-allocated in the cache.
         *
         * \warning If the output result @p D_phi is not correctly pre-allocated,
         * an exception will be raised.
         */
        void evaluate_bspline_values(
            const  ComponentContainer<TensorProductFunctionEvaluator<dim>> &elem_values,
            ValueTable<Value> &D_phi) const;

        /**
         * Computes the k-th order derivative of the non-zero B-spline basis
         * functions over the current element,
         *   at the evaluation points pre-allocated in the cache.
         *
         * \warning If the output result @p D_phi is not correctly pre-allocated,
         * an exception will be raised.
         */
        template <int order>
        void evaluate_bspline_derivatives(
            const  ComponentContainer<TensorProductFunctionEvaluator<dim>> &elem_values,
            ValueTable<Derivative<order>> &D_phi) const;



        GridElementHandler<dim_> *grid_handler_;
        int j_;
        const CacheList<GlobalCache, dim> *splines1d_;
        ReferenceElement<dim_,range_,rank_> *elem_;

    private:
        void
        copy_to_inactive_components_values(const vector<Index> &inactive_comp,
                                           const std::array<Index, n_components> &active_map,
                                           ValueTable<Value> &D_phi) const;

        template <int order>
        void
        copy_to_inactive_components(const vector<Index> &inactive_comp,
                                    const std::array<Index, n_components> &active_map,
                                    ValueTable<Derivative<order>> &D_phi) const;

    };

    FillCacheDispatcher fill_cache_impl_;

    /**
     * Returns the BSplineSpace used to define the BSplineElementHandler object.
     */
    std::shared_ptr<const Space> get_bspline_space() const;

};

IGA_NAMESPACE_CLOSE


#endif
