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
 *
 * @ingroup serializable
 */
template<int dim_, int range_, int rank_>
class BSplineElementHandler
    : public ReferenceElementHandler<dim_,range_,rank_>
{
    using base_t = ReferenceElementHandler<dim_,range_,rank_>;
    using self_t = BSplineElementHandler<dim_,range_,rank_>;
    using Space = BSplineSpace<dim_,range_,rank_>;
    static const Size n_components =  SplineSpace<dim_,range_,rank_>::n_components;

    using IndexType = typename CartesianGrid<dim_>::IndexType;

    template<class T>
    using ComponentContainer = typename Space::template ComponentContainer<T>;

    template<class T>
    using ComponentDirectionTable = ComponentContainer<CartesianProductArray<T,dim_>>;

    template<class T>
    using ComponentDirectionContainer = ComponentContainer<SafeSTLArray<T,dim_>>;

    using TensorSizeTable = typename Space::TensorSizeTable;

    template <int order>
    using Derivative = typename Space::template Derivative<order>;

    using Value = typename Space::Value;


protected:

    using BaseSpace = ReferenceSpace<dim_,range_,rank_>;
    using RefElementIterator = typename BaseSpace::ElementIterator;
    using RefElementAccessor = typename BaseSpace::ElementAccessor;


private:
    /**
     * Assignment operators.
     */
    ///@{
    /**
     * Copy assignment operator. Not allowed to be used.
     */
    self_t &operator=(const self_t &) = delete;

    /**
     * Move assignment operator. Not allowed to be used.
     */
    self_t &operator=(self_t &&) = delete;
    ///@}

    /**
     * @name Constructors.
     */
    ///@{

    /**
     * Default constructor. It does nothing but it is needed for the
     * <a href="http://www.boost.org/doc/libs/release/libs/serialization/">boost::serialization</a>
     * mechanism.
     */
    BSplineElementHandler() = default;

    BSplineElementHandler(std::shared_ptr<const Space> space);

    /**
     * Copy constructor. Not allowed to be used.
     */
    BSplineElementHandler(const self_t &) = delete;

    /**
     * Move constructor. Not allowed to be used.
     */
    BSplineElementHandler(self_t &&) = delete;
    ///@}

public:
    static const int dim = dim_;

    /**
     * Destructor.
     */
    virtual ~BSplineElementHandler() = default;

    static std::shared_ptr<self_t> create(std::shared_ptr<const Space> space);

    using topology_variant = typename base_t::topology_variant;
    using eval_pts_variant = typename base_t::eval_pts_variant;

#if 0
    /**
     * @name Reset functions
     */
    ///@{
    virtual void reset_selected_elements(
        const ValueFlags &flag,
        const eval_pts_variant &eval_points,
        const SafeSTLVector<IndexType> &elements_id) override final;
    ///@}

private:
    virtual void init_ref_elem_cache(RefElementAccessor &elem,
                                     const topology_variant &topology) override final;

    virtual void fill_ref_elem_cache(RefElementAccessor &elem,
                                     const topology_variant &topology,
                                     const int sub_elem_id) override final;
#endif

public:
    virtual void print_info(LogStream &out) const override final ;


private:

    static void
    fill_interval_values(const Real one_len,
                         const BernsteinOperator &oper,
                         const BasisValues1d &bernstein_vals,
                         BasisValues1d &spline_vals);


    static void
    resize_and_fill_bernstein_values(
        const int deg,
        const SafeSTLVector<Real> &pt_coords,
        BasisValues1d &bernstein_values);




    /**
     * One-dimensional B-splines values and derivatives at quadrature points.
     * The values are stored with the following index ordering:
     *
     * splines1d_[comp][dir][interval][order][function][point]
     *
     * @ingroup serializable
     */
    class GlobalCache
    {
    private:
        using BasisValues1dTable = ComponentContainer<SafeSTLArray<std::map<Index,BasisValues1d>,dim>>;


        /**
         * Quadrature points used for the 1D basis evaluation.
         */
        std::shared_ptr<const Quadrature<dim>> quad_;

        /**
         * Values (and derivatives) of 1D basis precomputed in the initalized
         * interval of a given direction.
         *
         * @note The map's key is the interval id. In Debug mode, it will be
         * raised an assertion if
         * the requested values are not initialized for the interval.
         *
         */
        BasisValues1dTable basis_values_1d_table_;


    public:
        using ComponentMap = typename BasisValues1dTable::ComponentMap;

        GlobalCache() = default;

        GlobalCache(const std::shared_ptr<const Quadrature<dim>> &quad, const ComponentMap &component_map);


        ComponentContainer<std::unique_ptr<const TensorProductFunctionEvaluator<dim>> >
                get_element_values(const TensorIndex<dim> &elem_tensor_id) const;

        BasisValues1d &entry(const int comp, const int dir, const Index interval_id);

        void print_info(LogStream &out) const;

    private:
#ifdef SERIALIZATION
        /**
         * @name Functions needed for boost::serialization
         * @see <a href="http://www.boost.org/doc/libs/release/libs/serialization/">boost::serialization</a>
         */
        ///@{
        friend class boost::serialization::access;

        template<class Archive>
        void
        serialize(Archive &ar, const unsigned int version)
        {
            ar &boost::serialization::make_nvp("quad_",quad_);
            ar &boost::serialization::make_nvp("basis_values_1d_table_",basis_values_1d_table_);
        }
        ///@}
#endif // SERIALIZATION
    };


#if 0
    struct ResetDispatcher : boost::static_visitor<void>
    {
        ResetDispatcher(const Space &space,
                        const ValueFlags flag_in,
                        GridElementHandler<dim_> &grid_handler,
                        SafeSTLArray<ValueFlags, dim+1> &flags,
                        CacheList<GlobalCache, dim> &splines1d)
            :
            space_(space),
            flag_(flag_in),
            grid_handler_(grid_handler),
            flags_(flags),
            splines1d_(splines1d)
        {}

        template<int sub_elem_dim>
        void operator()(const std::shared_ptr<const Quadrature<sub_elem_dim>> &quad);

        const Space &space_;
        const ValueFlags flag_;
        GridElementHandler<dim_> &grid_handler_;
        SafeSTLArray<ValueFlags, dim+1> &flags_;
        CacheList<GlobalCache, dim> &splines1d_;

        /**
         * id of the intervals that must be processed
         */
        SafeSTLArray<SafeSTLVector<int>,dim> intervals_id_directions_;
    };


    struct InitCacheDispatcher : boost::static_visitor<void>
    {
        InitCacheDispatcher(GridElementHandler<dim_> &grid_handler,
                            ReferenceElement<dim_,range_,rank_> &elem,
                            SafeSTLArray<ValueFlags, dim+1> &flags)
            :
            grid_handler_(grid_handler),
            elem_(elem),
            flags_(flags)
        {}

        template<int sub_elem_dim>
        void operator()(const Topology<sub_elem_dim> &sub_elem);

        GridElementHandler<dim_> &grid_handler_;
        ReferenceElement<dim_,range_,rank_> &elem_;
        SafeSTLArray<ValueFlags, dim+1> &flags_;
    };



    struct FillCacheDispatcher : boost::static_visitor<void>
    {
        FillCacheDispatcher(const int sub_elem_id,
                            const CacheList<GlobalCache, dim> &splines1d,
                            GridElementHandler<dim_> &grid_handler,
                            ReferenceElement<dim_,range_,rank_> &elem)
            :
            j_(sub_elem_id),
            splines1d_(splines1d),
            grid_handler_(grid_handler),
            elem_(elem)
        {}

        template<int sub_elem_dim>
        void operator()(const Topology<sub_elem_dim> &sub_elem);

        /**
         * Computes the values (i.e. the 0-th order derivative) of the non-zero
         *  B-spline basis
         * functions over the current element,
         *   at the evaluation points pre-allocated in the cache.
         *
         * \warning If the output result @p D_phi is not correctly
         * pre-allocated,
         * an exception will be raised.
         */
        void evaluate_bspline_values(
            const ComponentContainer<std::unique_ptr<const TensorProductFunctionEvaluator<dim>>> &elem_values,
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
            const ComponentContainer<std::unique_ptr<const TensorProductFunctionEvaluator<dim>>> &elem_values,
            ValueTable<Derivative<order>> &D_phi) const;



        const int j_;
        const CacheList<GlobalCache, dim> &splines1d_;
        GridElementHandler<dim_> &grid_handler_;
        ReferenceElement<dim_,range_,rank_> &elem_;

    private:
        void
        copy_to_inactive_components_values(const SafeSTLVector<Index> &inactive_comp,
                                           const SafeSTLArray<Index, n_components> &active_map,
                                           ValueTable<Value> &D_phi) const;

        template <int order>
        void
        copy_to_inactive_components(const SafeSTLVector<Index> &inactive_comp,
                                    const SafeSTLArray<Index, n_components> &active_map,
                                    ValueTable<Derivative<order>> &D_phi) const;

    };
#endif

    /**
     * Returns the BSplineSpace used to define the BSplineElementHandler object.
     */
    std::shared_ptr<const Space> get_bspline_space() const;


private:

    SafeSTLArray<ValueFlags, dim_ + 1> flags_;

    CacheList<GlobalCache, dim> splines1d_;

#ifdef SERIALIZATION
    /**
     * @name Functions needed for boost::serialization
     * @see <a href="http://www.boost.org/doc/libs/release/libs/serialization/">boost::serialization</a>
     */
    ///@{
    friend class boost::serialization::access;

    template<class Archive>
    void
    serialize(Archive &ar, const unsigned int version);
    ///@}
#endif
};

IGA_NAMESPACE_CLOSE


#endif
