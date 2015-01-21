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

#ifndef BSPLINE_ELEMENT_HANDLER_H_
#define BSPLINE_ELEMENT_HANDLER_H_

#include <igatools/base/config.h>

//TODO(pauletti, Sep 9, 2014): should we instantiate the cartesian product instead
#include <igatools/utils/cartesian_product_array-template.h>

#include <igatools/basis_functions/reference_element_handler.h>
#include <igatools/basis_functions/bspline_space.h>
#include <igatools/basis_functions/bspline_element_scalar_evaluator.h>
#include <igatools/basis_functions/bernstein_basis.h>


IGA_NAMESPACE_OPEN

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
    static const Size n_components =  SplineSpace<dim_,range_,rank_>::n_components;


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

//    using ElementIterator = typename Space::ElementIterator;
//    using ElementAccessor = typename Space::ElementAccessor;



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
    using eval_pts_variant = typename base_t::eval_pts_variant;

    virtual void reset(const ValueFlags &flag, const quadrature_variant &quad) override final;

    virtual void reset_one_element(const ValueFlags &flag, const eval_pts_variant &eval_points, const int elem_flat_id) override final;

    virtual void init_cache(RefElementAccessor &elem, const topology_variant &topology) override final;

    template <int k>
    void init_cache(RefElementAccessor &elem)
    {
        this->init_cache(elem,Int<k>());
    }

    template <int k>
    void init_cache(RefElementIterator &elem)
    {
        this->template init_cache<k>(*elem);
    }


    virtual void fill_cache(RefElementAccessor &elem, const topology_variant &topology, const int j) override final;

    template<int k>
    void fill_cache(RefElementAccessor &elem, const int j)
    {
        this->fill_cache(elem,Int<k>(),j);
    }

    template<int k>
    void fill_cache(RefElementIterator &elem, const int j)
    {
        this->template fill_cache<k>(*elem,j);
    }

    virtual void print_info(LogStream &out) const override final ;


private:

    static void
    fill_interval_values(const Real one_len,
                         const BernsteinOperator &oper,
                         const BasisValues1d &bernstein_vals,
                         BasisValues1d &spline_vals)
    {
        for (int order = 0; order < max_der; ++order)
        {
            auto &spline = spline_vals.get_derivative(order);
            const auto &berns = bernstein_vals.get_derivative(order);
            spline = oper.scale_action(std::pow(one_len, order), berns);
        }
    }


    static void
    resize_and_fill_bernstein_values(
        const int deg,
        const vector<Real> &pt_coords,
        BasisValues1d &bernstein_values)
    {
        bernstein_values.resize(max_der, deg+1, pt_coords.size());
        for (int order = 0; order < max_der; ++order)
            bernstein_values.get_derivative(order) =
                BernsteinBasis::derivative(order, deg, pt_coords);
    }


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

        /**
         * id of the intervals that must be processed
         */
        std::array<vector<int>,dim> intervals_id_directions_;
    };

    ResetDispatcher reset_impl_;


    struct ResetDispatcherOneElem : boost::static_visitor<void>
    {
        template<class T>
        void operator()(const T &eval_pts);

        GridElementHandler<dim_> *grid_handler_;
        ValueFlags flag_;
        std::array<FunctionFlags, dim + 1> *flags_;
        CacheList<GlobalCache, dim> *splines1d_;
        const Space *space_;

        int elem_flat_id_;
    };
    ResetDispatcherOneElem reset_one_elem_impl_;

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
