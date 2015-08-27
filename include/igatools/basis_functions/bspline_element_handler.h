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

    virtual void set_flags_impl(const typename space_element::Flags &flag, const topology_variant &topology) override final
    {
        auto set_flag_dispatcher = SetFlagDispatcher(flag,this->grid_handler_,flags_);
        boost::apply_visitor(set_flag_dispatcher,topology);
    }

    struct SetFlagDispatcher : boost::static_visitor<void>
    {
        SetFlagDispatcher(const typename space_element::Flags flag_in,
                          GridElementHandler<dim_> &grid_handler,
                          SafeSTLArray<typename space_element::Flags, dim+1> &flags)
            :
            flag_in_(flag_in),
            grid_handler_(grid_handler),
            flags_(flags)
        {}

        template<int sdim>
        void operator()(const Topology<sdim> &topology)
        {
            using GridFlags = grid_element::Flags;
            //TODO (martinelli, Aug 27, 2015): select the proper grid flags depending on the BSpline element flags
            const auto grid_flags = GridFlags::point |
                                    GridFlags::w_measure;
            grid_handler_.template set_flags<sdim>(grid_flags);

            flags_[sdim] = flag_in_;
        }

        const typename space_element::Flags flag_in_;
        GridElementHandler<dim_> &grid_handler_;
        SafeSTLArray<typename space_element::Flags, dim+1> &flags_;
    };


    using BaseElem = SpaceElement<dim_,0,range_,rank_,Transformation::h_grad>;

    virtual void init_cache_impl(BaseElem &elem,
                                 const eval_pts_variant &quad) const override final
    {
        auto init_cache_dispatcher = InitCacheDispatcher(this->grid_handler_,flags_,elem);
        boost::apply_visitor(init_cache_dispatcher,quad);
    }

    struct InitCacheDispatcher : boost::static_visitor<void>
    {
        InitCacheDispatcher(const GridElementHandler<dim_> &grid_handler,
                            const SafeSTLArray<typename space_element::Flags, dim+1> &flags,
                            BaseElem &elem)
            :
            grid_handler_(grid_handler),
            flags_(flags),
            elem_(elem)
        {}

        template<int sdim>
        void operator()(const std::shared_ptr<const Quadrature<sdim>> &quad)
        {
            grid_handler_.template init_cache<sdim>(elem_.get_grid_element(),quad);

            auto &cache = elem_.get_all_sub_elems_cache();
            if (cache == nullptr)
            {
                using VCache = typename BSplineElement<dim_,range_,rank_>::parent_t::Cache;

                using Cache = AllSubElementsCache<VCache>;
                cache = std::make_shared<Cache>();
            }

            const auto n_basis = elem_.get_max_num_basis();
            const auto n_points = quad->get_num_points();
            const auto flag = flags_[sdim];

            for (auto &s_id: UnitElement<dim_>::template elems_ids<sdim>())
            {
                auto &s_cache = cache->template get_sub_elem_cache<sdim>(s_id);
                s_cache.resize(flag, n_points, n_basis);
            }
        }

        const GridElementHandler<dim_> &grid_handler_;
        const SafeSTLArray<typename space_element::Flags, dim+1> &flags_;
        BaseElem &elem_;
    };


    virtual void fill_cache_impl(BaseElem &elem,
                                 const topology_variant &topology,
                                 const int s_id) const override final
    {
        auto fill_cache_dispatcher = FillCacheDispatcherNoGlobalCache(s_id,this->grid_handler_,elem);
        boost::apply_visitor(fill_cache_dispatcher,topology);
    }


    struct FillCacheDispatcherNoGlobalCache : boost::static_visitor<void>
    {
        FillCacheDispatcherNoGlobalCache(const int s_id,
                                         const GridElementHandler<dim_> &grid_handler,
                                         BaseElem &elem)
            :
            s_id_(s_id),
            grid_handler_(grid_handler),
            elem_(elem)
        {}

        template<int sdim>
        void operator()(const Topology<sdim> &topology)
        {
            auto &grid_elem = elem_.get_grid_element();
            grid_handler_.template fill_cache<sdim>(grid_elem,s_id_);



            //--------------------------------------------------------------------------------------
            // filling the 1D cache --- begin

            const auto &grid = *grid_elem.get_grid();
            const auto n_inter = grid.get_num_intervals();

            const auto elem_size = grid_elem.template get_side_lengths<dim>(0);
            const auto elem_tensor_id = grid_elem.get_index();

            const auto &bsp_space = dynamic_cast<const Space &>(*elem_.get_space());

            const auto &space_data = *bsp_space.space_data_;

            const auto &degree = bsp_space.get_degree_table();

            const auto &active_components_id = space_data.get_active_components_id();

            const auto quad_in_cache = grid_elem.template get_quadrature<sdim>();

            const auto quad = extend_sub_elem_quad<sdim,dim>(*quad_in_cache, s_id_);

            using BasisValues1dTable = ComponentContainer<SafeSTLArray<BasisValues1d,dim>>;
            BasisValues1dTable splines_derivatives_1D_table(space_data.get_components_map());

            const auto &n_coords = quad.get_num_coords_direction();

            /*
             * For each direction, interval and component we compute the 1D bspline
             * basis evaluate at the 1D component of the tensor product quadrature
             */
            const auto &bezier_op   = bsp_space.operators_;
            const auto &end_interval = bsp_space.end_interval_;

            using BasisValues = ComponentContainer<BasisValues1d>;
            const auto &deg_comp_map = degree.get_comp_map();
            BasisValues bernstein_values(deg_comp_map);

            for (const int dir : UnitElement<dim>::active_directions)
            {
                const auto &pt_coords_internal = quad.get_coords_direction(dir);

                const auto len = elem_size[dir];

                const auto interval_id = elem_tensor_id[dir];

                Real alpha;

                const int n_pts_1D = n_coords[dir];

                SafeSTLVector<Real> pt_coords_boundary(n_pts_1D);

                const SafeSTLVector<Real> *pt_coords_ptr = nullptr;

                for (auto comp : active_components_id)
                {
                    const int deg = degree[comp][dir];

                    auto &splines_derivatives_1D = splines_derivatives_1D_table[comp][dir];
                    splines_derivatives_1D.resize(MAX_NUM_DERIVATIVES,deg+1,n_pts_1D);

                    if (interval_id == 0) // processing the leftmost interval
                    {
                        // first interval (i.e. left-most interval)

                        alpha = end_interval[comp][dir].first;
                        const Real one_minus_alpha = 1. - alpha;

                        for (int ipt = 0 ; ipt < n_pts_1D ; ++ipt)
                            pt_coords_boundary[ipt] = one_minus_alpha +
                                                      pt_coords_internal[ipt] * alpha;

                        pt_coords_ptr = &pt_coords_boundary;
                    } // end process_interval_left
                    else if (interval_id == n_inter[dir]-1) // processing the rightmost interval
                    {
                        // last interval (i.e. right-most interval)

                        alpha = end_interval[comp][dir].second;

                        for (int ipt = 0 ; ipt < n_pts_1D ; ++ipt)
                            pt_coords_boundary[ipt] = pt_coords_internal[ipt] *
                                                      alpha;

                        pt_coords_ptr = &pt_coords_boundary;
                    } // end process_interval_right
                    else
                    {
                        // internal interval

                        alpha = 1.0;

                        pt_coords_ptr = &pt_coords_internal;
                    } // end process_interval_internal


                    const Real alpha_div_interval_length = alpha / len;

                    const auto &oper = bezier_op.get_operator(dir,interval_id,comp);

                    //------------------------------------------------------------
                    //resize_and_fill_bernstein_values
                    bernstein_values[comp].resize(MAX_NUM_DERIVATIVES,deg+1,n_pts_1D);
                    for (int order = 0; order < MAX_NUM_DERIVATIVES; ++order)
                    {
                        auto &berns = bernstein_values[comp].get_derivative(order);
                        berns = BernsteinBasis::derivative(order, deg,*pt_coords_ptr);

                        auto &splines = splines_derivatives_1D.get_derivative(order);
                        splines = oper.scale_action(
                                      std::pow(alpha_div_interval_length, order),
                                      berns);
                    }
                    //------------------------------------------------------------

                } // end loop comp

            } // end loop dir
            //
            // filling the 1D cache --- end
            //-------------------------------------------------------------------------------



            //-------------------------------------------------------------------------------
            // Multi-variate spline evaluation from 1D values --- begin
            using TPFE = const TensorProductFunctionEvaluator<dim>;
            ComponentContainer<std::unique_ptr<TPFE>> val_1d(splines_derivatives_1D_table.get_comp_map());

            SafeSTLArray<BasisValues1dConstView, dim> values_1D;
            for (auto c : val_1d.get_active_components_id())
            {
                const auto &value = splines_derivatives_1D_table[c];
                for (int i = 0 ; i < dim_ ; ++i)
                    values_1D[i] = BasisValues1dConstView(value[i]);

                val_1d[c] = std::make_unique<TPFE>(quad,values_1D);
            }
            // Multi-variate spline evaluation from 1D values --- end
            //-------------------------------------------------------------------------------



            auto &all_sub_elems_cache = elem_.get_all_sub_elems_cache();
            Assert(all_sub_elems_cache != nullptr, ExcNullPtr());
            auto &sub_elem_cache = all_sub_elems_cache->template get_sub_elem_cache<sdim>(s_id_);
//#if 0
//            const auto val_1d = g_cache.get_element_values();

            using Elem = SpaceElement<dim_,0,range_,rank_,Transformation::h_grad>;
            using _Value      = typename Elem::_Value;
            using _Gradient   = typename Elem::_Gradient;
            using _Hessian    = typename Elem::_Hessian;
            using _Divergence = typename Elem::_Divergence;

            if (sub_elem_cache.template status_fill<_Value>())
            {
                auto &values = sub_elem_cache.template get_data<_Value>();
                evaluate_bspline_values(val_1d, values);
                sub_elem_cache.template set_status_filled<_Value>(true);
            }
            if (sub_elem_cache.template status_fill<_Gradient>())
            {
                auto &values = sub_elem_cache.template get_data<_Gradient>();
                evaluate_bspline_derivatives<1>(val_1d, values);
                sub_elem_cache.template set_status_filled<_Gradient>(true);
            }
            if (sub_elem_cache.template status_fill<_Hessian>())
            {
                auto &values = sub_elem_cache.template get_data<_Hessian>();
                evaluate_bspline_derivatives<2>(val_1d, values);
                sub_elem_cache.template set_status_filled<_Hessian>(true);
            }
            if (sub_elem_cache.template status_fill<_Divergence>())
            {
                eval_divergences_from_gradients(
                    sub_elem_cache.template get_data<_Gradient>(),
                    sub_elem_cache.template get_data<_Divergence>());
                sub_elem_cache.template set_status_filled<_Divergence>(true);
            }
//#endif
            sub_elem_cache.set_filled(true);
        }

//#if 0
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
//#endif

        const int s_id_;
        const GridElementHandler<dim_> &grid_handler_;
        BaseElem &elem_;

    private:

//#if 0
        void
        copy_to_inactive_components_values(const SafeSTLVector<Index> &inactive_comp,
                                           const SafeSTLArray<Index, n_components> &active_map,
                                           ValueTable<Value> &D_phi) const;

        template <int order>
        void
        copy_to_inactive_components(const SafeSTLVector<Index> &inactive_comp,
                                    const SafeSTLArray<Index, n_components> &active_map,
                                    ValueTable<Derivative<order>> &D_phi) const;
//#endif
    };


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



#if 0
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


        /**
         * Quadrature points used for the 1D basis evaluation.
         */
        std::shared_ptr<const Quadrature<dim>> quad_;


        using BasisValues1dTable = ComponentContainer<SafeSTLArray<std::map<Index,BasisValues1d>,dim>>;

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
#endif


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

    SafeSTLArray<typename space_element::Flags, dim_ + 1> flags_;

#if 0
    CacheList<GlobalCache, dim> splines1d_;
#endif

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
