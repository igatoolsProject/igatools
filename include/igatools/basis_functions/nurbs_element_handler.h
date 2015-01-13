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


#ifndef NURBS_ELEMENT_HANDLER_H_
#define NURBS_ELEMENT_HANDLER_H_

#include <igatools/base/config.h>

#ifdef NURBS

#include <igatools/base/cache_status.h>
#include <igatools/base/flags_handler.h>
#include <igatools/base/quadrature.h>

//TODO(pauletti, Sep 9, 2014): should we instantiate the cartesian product instead
//#include <igatools/utils/cartesian_product_array-template.h>
//#include <igatools/basis_functions/nurbs_space.h>

#include <igatools/utils/value_table.h>
#include <igatools/geometry/grid_element_handler.h>
#include <igatools/basis_functions/bspline_element_handler.h>




IGA_NAMESPACE_OPEN


template<int,int,int> class NURBSSpace;
template<int,int,int> class NURBSElement;

/**
 * Global NURBSSpace uniform quadrature
 * computational optimization cache.
 */
template<int dim_, int range_ = 1, int rank_ = 1>
class NURBSElementHandler : public ReferenceElementHandler<dim_,range_,rank_>
{
    using base_t = ReferenceElementHandler<dim_,range_,rank_>;
    using self_t = NURBSElementHandler<dim_,range_,rank_>;
    using Space = NURBSSpace<dim_,range_,rank_>;
    static const Size n_components =  Space::n_components;

    template<class T>
    using ComponentContainer = typename Space::template ComponentContainer<T>;

    template <int order>
    using Derivative = typename Space::template Derivative<order>;

    using Value = typename Space::Value;

protected:
    using ElementIterator = typename Space::ElementIterator;
    using ElementAccessor = typename Space::ElementAccessor;

    using BaseSpace = ReferenceSpace<dim_,range_,rank_>;
    using RefElementIterator = typename BaseSpace::ElementIterator;
    using RefElementAccessor = typename BaseSpace::ElementAccessor;


public:
    static const int dim = dim_;



private:
    NURBSElementHandler(std::shared_ptr<const Space> space);

public:
    virtual ~NURBSElementHandler() = default;


    static std::shared_ptr<self_t> create(std::shared_ptr<const Space> space)
    {
        return std::shared_ptr<self_t>(new self_t(space));
    }

    using quadrature_variant = typename base_t::quadrature_variant;
    using topology_variant = typename base_t::topology_variant;

    virtual void reset(const ValueFlags &flag, const quadrature_variant &quad) override final;

    virtual void reset(const ValueFlags &flag, const ValueVector<typename Space::RefPoint> &points) override final
    {
        Assert(false,ExcNotImplemented());
        AssertThrow(false,ExcNotImplemented());
    }

    virtual void init_cache(RefElementAccessor &elem, const topology_variant &topology) override final;

    virtual void fill_cache(RefElementAccessor &elem, const topology_variant &topology, const int j) override final;


    virtual void print_info(LogStream &out) const override final;

private:
//    std::shared_ptr<const Space> space_;
    std::shared_ptr<BSplineElementHandler<dim_,range_,rank_>> bspline_handler_;

    std::array<FunctionFlags, dim + 1> flags_;


    using WeightElem = typename Space::WeightFunction::ElementAccessor;
    using WeightElemTable = typename Space::template ComponentContainer<WeightElem>;



#if 0
    /**
     * Returns the active components id for the NURBS values and derivatives.
     *
     * @note The active components id is the union of the active components for the numerator
     * (basis function belonging to a BSplineSpace) and the active components for the denominator
     * (a ComponentTable of scalar IgFunction(s)).
     */
    vector<int> get_active_components_id() const;


    /**
     * Returns the active components id for the NURBS values and derivatives.
     *
     * @note The incative components id are the complement of the active components id with respect
     * to the sequence 0,1,...,n_components-1
     *
     * @see get_active_components_id()
     */
    vector<int> get_inactive_components_id() const;
#endif


    struct ResetDispatcher : boost::static_visitor<void>
    {
        template<class T>
        void operator()(const T &quad);

        GridElementHandler<dim_> *grid_handler_;
        ValueFlags flag_;
        std::array<FunctionFlags, dim + 1> *flags_;
    };

    ResetDispatcher reset_impl_;

    struct InitCacheDispatcher : boost::static_visitor<void>
    {
        template<class T>
        void operator()(const T &quad);

        GridElementHandler<dim_> *grid_handler_;
        ReferenceElement<dim_,range_,rank_> *elem_;
        int n_points_;
        std::array<FunctionFlags, dim + 1> *flags_;

    };

    InitCacheDispatcher init_cache_impl_;

    struct FillCacheDispatcher : boost::static_visitor<void>
    {
        template<class T>
        void operator()(const T &quad);

        /**
         * Computes the value of the non-zero NURBS basis
         * functions over the current element,
         *   at the evaluation points pre-allocated in the cache.
         *
         * \warning If the output result @p derivatives_phi_hat is not correctly pre-allocated,
         * an exception will be raised.
         */
        void evaluate_nurbs_values_from_bspline(
            const typename Space::SpSpace::ElementAccessor &bspline_elem,
            const WeightElemTable &weight_elem_table,
            ValueTable<Value> &phi) const;

        /**
         * Computes the 1st order derivative of the non-zero NURBS basis
         * functions over the current element,
         *   at the evaluation points pre-allocated in the cache.
         *
         * \warning If the output result @p derivatives_phi_hat is not correctly pre-allocated,
         * an exception will be raised.
         */
        void evaluate_nurbs_gradients_from_bspline(
            const typename Space::SpSpace::ElementAccessor &bspline_elem,
            const WeightElemTable &weight_elem_table,
            ValueTable<Derivative<1>> &D1_phi) const;

        /**
         * Computes the 2nd order derivative of the non-zero NURBS basis
         * functions over the current element,
         *   at the evaluation points pre-allocated in the cache.
         *
         * \warning If the output result @p derivatives_phi_hat is not correctly pre-allocated,
         * an exception will be raised.
         */
        void evaluate_nurbs_hessians_from_bspline(
            const typename Space::SpSpace::ElementAccessor &bspline_elem,
            const WeightElemTable &weight_elem_table,
            ValueTable<Derivative<2>> &D2_phi) const;

        int j_;
        NURBSElement<dim_,range_,rank_> *nrb_elem_;
    };

    FillCacheDispatcher fill_cache_impl_;



    /**
     * Returns the NURBSSpace used to define the NURBSElementHandler object.
     */
    std::shared_ptr<const Space> get_nurbs_space() const;
};




IGA_NAMESPACE_CLOSE

#endif // #ifdef NURBS

#endif // #ifndef NURBS_ELEMENT_HANDLER_H_

