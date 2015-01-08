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

#include <igatools/basis_functions/nurbs_element_handler.h>
#include <igatools/basis_functions/nurbs_element.h>

#include <algorithm>

#ifdef NURBS
using std::shared_ptr;
IGA_NAMESPACE_OPEN

template<int dim_, int range_ , int rank_>
NURBSElementHandler<dim_, range_, rank_>::
NURBSElementHandler(shared_ptr<const Space> space)
    :
    base_t(space)
//  ,
//    space_(space)
//  ,
//    bspline_handler_(BSplineElementHandler<dim_,range_,rank_>::create(space->get_spline_space()))
{
    const auto bsp_space = space->get_spline_space();
    bspline_handler_ = BSplineElementHandler<dim_,range_,rank_>::create(bsp_space);
}


#if 0
template<int dim_, int range_ , int rank_>
NURBSElementHandler<dim_, range_, rank_>::
NURBSElementHandler(shared_ptr<const Space> space)
    :
    base_t(space->get_grid()),
    space_(space),
    n_basis_(space_->get_num_all_element_basis())
{

    // Compute the component offsets
    comp_offset_[0] = 0;
    for (int j = 1; j < Space::n_components; ++j)
        comp_offset_[j] = comp_offset_[j-1] + n_basis_.comp_dimension[j-1];
}
#endif


#if 0
template<int dim_, int range_ , int rank_>
template<int k>
void
NURBSElementHandler<dim_, range_, rank_>::
reset(const ValueFlags flag,
      const Quadrature<k> &quad1)
{
    //--------------------------------------
    // resetting the BSplineElementHandler (for the numerator)
    bspline_handler_->reset(flag, quad1);
    //--------------------------------------


    //--------------------------------------------------
    // resetting the Function for the weight (for the denominator)
    int max_deriv_order = -1;
    if (contains(flag, ValueFlags::point) ||
        contains(flag, ValueFlags::value))
        max_deriv_order = 0;

    if (contains(flag, ValueFlags::measure) ||
        contains(flag, ValueFlags::w_measure) ||
        contains(flag, ValueFlags::boundary_normal) ||
        contains(flag, ValueFlags::outer_normal) ||
        contains(flag, ValueFlags::gradient) ||
        contains(flag, ValueFlags::inv_gradient))
        max_deriv_order = 1;


    if (contains(flag, ValueFlags::curvature) ||
        contains(flag, ValueFlags::hessian) ||
        contains(flag, ValueFlags::inv_hessian))
        max_deriv_order = 2;


    ValueFlags weight_flag;
    if (max_deriv_order == 0)
        weight_flag = ValueFlags::value;
    else if (max_deriv_order == 1)
        weight_flag = ValueFlags::value | ValueFlags::gradient;
    else if (max_deriv_order == 2)
        weight_flag = ValueFlags::value | ValueFlags::gradient | ValueFlags::hessian;
    else
        Assert(false,ExcMessage("Not a right value flag."));

    for (const auto &comp_id : space_->weight_func_table_.get_active_components_id())
        space_->weight_func_table_[comp_id]->reset(weight_flag,quad1);
    //--------------------------------------------------



    flags_[k] = flag;
}
#endif

template<int dim_, int range_ , int rank_>
template<class T>
void
NURBSElementHandler<dim_, range_, rank_>::
ResetDispatcher::
operator()(const T &quad1)
{
    (*flags_)[T::dim] = flag_;
}


template<int dim_, int range_ , int rank_>
void
NURBSElementHandler<dim_, range_, rank_>::
reset(const ValueFlags &flag, const quadrature_variant &quad)
{
    //--------------------------------------
    // resetting the BSplineElementHandler (for the numerator)
    Assert(bspline_handler_ != nullptr, ExcNullPtr());
    bspline_handler_->reset(flag, quad);
    //--------------------------------------


    //--------------------------------------------------
    // resetting the Function for the weight (for the denominator)
    int max_deriv_order = -1;
    if (contains(flag, ValueFlags::point) ||
        contains(flag, ValueFlags::value))
        max_deriv_order = 0;

    if (contains(flag, ValueFlags::measure) ||
        contains(flag, ValueFlags::w_measure) ||
        contains(flag, ValueFlags::boundary_normal) ||
        contains(flag, ValueFlags::outer_normal) ||
        contains(flag, ValueFlags::gradient) ||
        contains(flag, ValueFlags::inv_gradient))
        max_deriv_order = 1;


    if (contains(flag, ValueFlags::curvature) ||
        contains(flag, ValueFlags::hessian) ||
        contains(flag, ValueFlags::inv_hessian))
        max_deriv_order = 2;


    ValueFlags weight_flag;
    if (max_deriv_order == 0)
        weight_flag = ValueFlags::value;
    else if (max_deriv_order == 1)
        weight_flag = ValueFlags::value | ValueFlags::gradient;
    else if (max_deriv_order == 2)
        weight_flag = ValueFlags::value | ValueFlags::gradient | ValueFlags::hessian;
    else
        Assert(false,ExcMessage("Not a right value flag."));


    const auto nrb_space = this->get_nurbs_space();
    for (const auto &comp_id : nrb_space->weight_func_table_.get_active_components_id())
        nrb_space->weight_func_table_[comp_id]->reset(weight_flag,quad);
    //--------------------------------------------------


    //--------------------------------------------------
//    reset_impl_.grid_handler_ = &(this->grid_handler_);
    reset_impl_.grid_handler_ =
        const_cast<GridElementHandler<dim_> *>(&(bspline_handler_->get_grid_handler()));
    reset_impl_.flag_ = flag;
    reset_impl_.flags_ = &flags_;

    boost::apply_visitor(reset_impl_, quad);
    //--------------------------------------------------
}


template<int dim_, int range_ , int rank_>
template<class T>
void
NURBSElementHandler<dim_, range_, rank_>::
InitCacheDispatcher::
operator()(const T &quad1)
{
    Assert(grid_handler_ != nullptr,ExcNullPtr());
    Assert(elem_ != nullptr,ExcNullPtr());
    grid_handler_->template init_cache<T::k>(*elem_);

    auto &cache = elem_->get_local_cache();
    if (cache == nullptr)
    {
        using Cache = typename ElementAccessor::LocalCache;
        cache = shared_ptr<Cache>(new Cache);
    }

    const auto n_basis = elem_->get_num_basis();
    const auto n_points = grid_handler_->template get_num_points<T::k>();
    const auto flag = (*flags_)[T::k];

    for (auto &s_id: UnitElement<dim>::template elems_ids<T::k>())
    {
        auto &s_cache = cache->template get_value_cache<T::k>(s_id);
        s_cache.resize(flag, n_points, n_basis);
    }
}

template<int dim_, int range_ , int rank_>
void
NURBSElementHandler<dim_, range_, rank_>::
init_cache(RefElementAccessor &elem, const topology_variant &topology)
{
    Assert(!elem.get_space()->is_bspline(),ExcMessage("Not a NURBSElement."));

    auto nrb_elem = dynamic_cast<NURBSElement<dim_,range_,rank_>*>(&elem);
    auto &bsp_elem = nrb_elem->bspline_elem_;
    bspline_handler_->init_cache(bsp_elem,topology);

    const auto nrb_space = this->get_nurbs_space();
    for (const auto &comp_id : nrb_space->weight_func_table_.get_active_components_id())
        nrb_space->weight_func_table_[comp_id]->init_cache(nrb_elem->weight_elem_table_[comp_id],topology);


    //-------------------------------------
    init_cache_impl_.grid_handler_ =
        const_cast<GridElementHandler<dim_> *>(&(bspline_handler_->get_grid_handler()));
    init_cache_impl_.elem_ = nrb_elem;
    init_cache_impl_.flags_ = &flags_;

    boost::apply_visitor(init_cache_impl_,topology);
    //-------------------------------------
}



template<int dim_, int range_ , int rank_>
template<class T>
void
NURBSElementHandler<dim_, range_, rank_>::
FillCacheDispatcher::
operator()(const T &quad1)
{
    Assert(nrb_elem_ != nullptr, ExcNullPtr());

    Assert(nrb_elem_->local_cache_ != nullptr, ExcNullPtr());
    auto &cache = nrb_elem_->local_cache_->template get_value_cache<T::k>(j_);

    const auto &bsp_elem = nrb_elem_->bspline_elem_;
    const auto &wght_table = nrb_elem_->weight_elem_table_;

    auto &flags = cache.flags_handler_;
    if (flags.fill_values())
    {
        auto &values = cache.template get_der<0>();
        evaluate_nurbs_values_from_bspline(bsp_elem, wght_table, values);
        flags.set_values_filled(true);
    }
    if (flags.fill_gradients())
    {
        auto &gradients = cache.template get_der<1>();
        evaluate_nurbs_gradients_from_bspline(bsp_elem, wght_table, gradients);
        flags.set_gradients_filled(true);
    }
    if (flags.fill_hessians())
    {
        auto &hessians = cache.template get_der<2>();
        evaluate_nurbs_hessians_from_bspline(bsp_elem, wght_table, hessians);
        flags.set_hessians_filled(true);
    }
    cache.set_filled(true);
}


template<int dim_, int range_ , int rank_>
void
NURBSElementHandler<dim_, range_, rank_>::
fill_cache(RefElementAccessor &elem, const topology_variant &topology, const int j)
{
    using Elem = NURBSElement<dim_,range_,rank_>;
    const auto nrb_elem = dynamic_cast<Elem *>(&elem);
    Assert(nrb_elem != nullptr, ExcNullPtr());

    const auto nrb_space = this->get_nurbs_space();
    Assert(nrb_space == nrb_elem->get_nurbs_space(),
           ExcMessage("The element accessor and the element handler cannot have different spaces."));

    bspline_handler_->fill_cache(nrb_elem->bspline_elem_,topology,j);

    const auto &weight_func_table = nrb_space->weight_func_table_;
    for (const auto &comp_id : weight_func_table.get_active_components_id())
        weight_func_table[comp_id]->fill_cache(nrb_elem->weight_elem_table_[comp_id],topology,j);


    //-----------------------------------------
    fill_cache_impl_.j_ = j;
    fill_cache_impl_.nrb_elem_ = nrb_elem;
    boost::apply_visitor(fill_cache_impl_,topology);
    //-----------------------------------------
}




template<int dim_, int range_ , int rank_>
auto
NURBSElementHandler<dim_, range_, rank_>::
get_nurbs_space() const -> std::shared_ptr<const Space>
{
    auto nrb_space = std::dynamic_pointer_cast<const Space>(this->get_space());
    Assert(nrb_space != nullptr,ExcNullPtr());
    return nrb_space;
}

template<int dim_, int range_ , int rank_>
void
NURBSElementHandler<dim_, range_, rank_>::
print_info(LogStream &out) const
{
    Assert(false,ExcNotImplemented());
    AssertThrow(false,ExcNotImplemented());
#if 0

    out.begin_item("Grid Cache:");
    base_t::print_info(out);
    out.end_item();

    cacheutils::print_caches(splines1d_, out);
    //out.begin_item("One dimensional splines cache:");
    //splines1d_.print_info(out);
    //out.end_item();
#endif
}

template<int dim_, int range_ , int rank_>
void
NURBSElementHandler<dim_, range_, rank_>::
FillCacheDispatcher::
evaluate_nurbs_values_from_bspline(
    const typename Space::SpSpace::ElementAccessor &bspline_elem,
    const WeightElemTable &weight_elem_table,
    ValueTable<Value> &phi) const
{
    /*
     * This function evaluates the values of the n+1 NURBS basis function R_0,...,R_n
     * from the set of BSpline basis function N_0,...,N_n
     * where the i-th NURBS basis function is defined as
     *
     *         P_i
     * R_i = -------
     *          Q
     *
     * with P_i a basis function of a BSplineSpace
     * and Q an IgFunction built over a scalar BSpline space
     *
     */

    Assert(!phi.empty(), ExcEmptyObject());

    const auto &P = bspline_elem.template get_values<0,dim>(0);
    const auto n_pts = P.get_num_points();

    const auto bsp_local_to_patch = bspline_elem.get_local_to_patch();

    const auto nrb_space = nrb_elem_->get_nurbs_space();
    const auto comp_offset = nrb_space->sp_space_->get_basis_offset();

    int bsp_fn_id = 0;
    for (int comp = 0 ; comp < n_components ; ++comp)
    {
        const auto &weight_elem = weight_elem_table[comp];

        const auto &Q = weight_elem.template get_values<0,dim>(0);

        Assert(n_pts == Q.get_num_points(),
               ExcDimensionMismatch(n_pts,Q.get_num_points()));

        vector<Real> invQ(n_pts);
        for (int pt = 0 ; pt < n_pts ; ++pt)
            invQ[pt] = 1.0 / Q[pt](0);

        const int n_funcs_comp = bspline_elem.get_num_basis(comp);

        const auto &w_coefs = nrb_space->weight_func_table_[comp]->get_coefficients();
        Assert(nrb_space->get_num_basis(comp) == w_coefs.size(),
               ExcDimensionMismatch(nrb_space->get_num_basis(comp),w_coefs.size()));

        const auto offset = comp_offset[comp];
        for (int w_fn_id = 0 ; w_fn_id < n_funcs_comp ; ++w_fn_id, ++bsp_fn_id)
        {
            const auto &P_fn = P.get_function_view(bsp_fn_id);

            auto R_fn = phi.get_function_view(bsp_fn_id);

            const Real w = w_coefs[bsp_local_to_patch[bsp_fn_id]-offset];

            for (int pt = 0 ; pt < n_pts ; ++pt)
                R_fn[pt](comp) = P_fn[pt](comp) * invQ[pt] * w ;
        } // end loop w_fn_id
    } // end loop comp
}

template<int dim_, int range_ , int rank_>
void
NURBSElementHandler<dim_, range_, rank_>::
FillCacheDispatcher::
evaluate_nurbs_gradients_from_bspline(
    const typename Space::SpSpace::ElementAccessor &bspline_elem,
    const WeightElemTable &weight_elem_table,
    ValueTable<Derivative<1>> &D1_phi) const
{
    /*
     * This function evaluates the gradients of the n+1 NURBS basis function R_0,...,R_n
     * from the set of BSpline basis function N_0,...,N_n
     * where the i-th NURBS basis function is defined as
     *
     *         P_i
     * R_i = -------
     *          Q
     *
     * with P_i a basis function of a BSplineSpace
     * and Q an IgFunction built over a scalar space.
     *
     * Then the gradient dR_i is:
     *
     *         dP_i       P_i * dQ
     * dR_i = ------- -  ----------
     *           Q           Q^2
     *
     */

    Assert(!D1_phi.empty(), ExcEmptyObject());

    const auto &P  = bspline_elem.template get_values<0,dim>(0);
    const auto &dP = bspline_elem.template get_values<1,dim>(0);

    const auto n_pts = P.get_num_points();

    const auto nrb_space = nrb_elem_->get_nurbs_space();
    const auto bsp_local_to_patch = bspline_elem.get_local_to_patch();
    const auto comp_offset = nrb_space->sp_space_->get_basis_offset();

    int bsp_fn_id = 0;
    for (int comp = 0 ; comp < n_components ; ++comp)
    {
        const auto &weight_elem = weight_elem_table[comp];

        const auto &Q  = weight_elem.template get_values<0,dim>(0);
        const auto &dQ = weight_elem.template get_values<1,dim>(0);

        Assert(n_pts == Q.get_num_points(),
               ExcDimensionMismatch(n_pts,Q.get_num_points()));


        vector<Real> invQ(n_pts);
        vector<array<Real,dim>> dQ_invQ2(n_pts);
        for (int pt = 0 ; pt < n_pts ; ++pt)
        {
            invQ[pt] = 1.0 / Q[pt](0);

            const auto &dQ_pt = dQ[pt];
            auto &dQ_invQ2_pt = dQ_invQ2[pt];

            for (int i = 0 ; i < dim ; ++i)
                dQ_invQ2_pt[i] = invQ[pt] * invQ[pt] * dQ_pt(i)(0);
        }

        const int n_funcs_comp = bspline_elem.get_num_basis(comp);

        const auto &w_coefs = nrb_space->weight_func_table_[comp]->get_coefficients();
        Assert(nrb_space->get_num_basis(comp) == w_coefs.size(),
               ExcDimensionMismatch(nrb_space->get_num_basis(comp),w_coefs.size()));

        const auto offset = comp_offset[comp];
        for (int w_fn_id = 0 ; w_fn_id < n_funcs_comp ; ++w_fn_id, ++bsp_fn_id)
        {
            const auto &P_fn  =  P.get_function_view(bsp_fn_id);
            const auto &dP_fn = dP.get_function_view(bsp_fn_id);

            auto dR_fn = D1_phi.get_function_view(bsp_fn_id);

            const Real w = w_coefs[bsp_local_to_patch[bsp_fn_id]-offset];

            for (int pt = 0 ; pt < n_pts ; ++pt)
            {
                auto &dR_fn_pt = dR_fn[pt];

                const auto   &P_fn_pt =  P_fn[pt];
                const auto &dP_fn_pt = dP_fn[pt];

                const Real invQ_pt = invQ[pt];
                const auto &dQ_invQ2_pt = dQ_invQ2[pt];

                for (int i = 0 ; i < dim ; ++i)
                    dR_fn_pt(i)(comp) = (dP_fn_pt(i)(comp)*invQ_pt - P_fn_pt(comp)*dQ_invQ2_pt[i]) * w;
            } // end loop pt
        } // end loop w_fn_id
    } // end loop comp
}

template<int dim_, int range_ , int rank_>
void
NURBSElementHandler<dim_, range_, rank_>::
FillCacheDispatcher::
evaluate_nurbs_hessians_from_bspline(
    const typename Space::SpSpace::ElementAccessor &bspline_elem,
    const WeightElemTable &weight_elem_table,
    ValueTable<Derivative<2>> &D2_phi) const
{
    /*
     * This function evaluates the gradients of the n+1 NURBS basis function R_0,...,R_n
     * from the set of BSpline basis function N_0,...,N_n
     * where the k-th NURBS basis function is defined as
     *
     *        Pk
     * Rk = -------
     *         Q
     *
     * with Pk a basis function of a BSplineSpace
     * and Q an IgFunction built over a scalar space.
     *
     * Then the gradient dRk is defined by the partial derivatives:
     *
     *          dPk_i     Pk * dQ_i
     * dRk_i = ------- - -----------
     *            Q          Q^2
     *
     * And the hessian d2Rk is:
     *                                _                                         _
     *            d2Pk_ij      1     |                                           |    2 * Pk
     * d2Rk_ij = --------- - ----- * | dPk_i * dQ_j + dPk_j * dQ_i + Pk * d2Q_ij | + -------- * dQ_i * dQ_j
     *               Q        Q^2    |_                                         _|      Q^3
     *
     */
    Assert(!D2_phi.empty(), ExcEmptyObject());

    const auto &P   = bspline_elem.template get_values<0,dim>(0);
    const auto &dP  = bspline_elem.template get_values<1,dim>(0);
    const auto &d2P = bspline_elem.template get_values<2,dim>(0);

    const auto n_pts = P.get_num_points();

    const auto nrb_space = nrb_elem_->get_nurbs_space();
    const auto bsp_local_to_patch = bspline_elem.get_local_to_patch();
    const auto comp_offset = nrb_space->sp_space_->get_basis_offset();

    int bsp_fn_id = 0;
    for (int comp = 0 ; comp < n_components ; ++comp)
    {
        const auto &weight_elem = weight_elem_table[comp];

        const auto &Q   = weight_elem.template get_values<0,dim>(0);
        const auto &dQ  = weight_elem.template get_values<1,dim>(0);
        const auto &d2Q = weight_elem.template get_values<2,dim>(0);

        Assert(n_pts == Q.get_num_points(),
               ExcDimensionMismatch(n_pts,Q.get_num_points()));


        vector<Real> invQ(n_pts);
        vector<Real> invQ2(n_pts);
        vector<array<Real,dim> > dQ_invQ2(n_pts);
        vector<array<array<Real,dim>,dim> > Q_terms_2nd_order(n_pts);
        for (int pt = 0 ; pt < n_pts ; ++pt)
        {
            auto &invQ_pt  = invQ[pt];
            auto &invQ2_pt = invQ2[pt];

            const auto &dQ_pt = dQ[pt];
            const auto &d2Q_pt = d2Q[pt];

            auto &dQ_invQ2_pt = dQ_invQ2[pt];
            auto &Q_terms_2nd_order_pt = Q_terms_2nd_order[pt];

            invQ_pt = 1.0 / Q[pt](0);
            invQ2_pt = invQ_pt * invQ_pt;

            const Real two_invQ_pt = 2.0 * invQ_pt;

            int hessian_entry_fid = 0;
            for (int i = 0 ; i < dim ; ++i)
            {
                dQ_invQ2_pt[i] = invQ2_pt * dQ_pt(i)(0);

                for (int j = 0 ; j < dim ; ++j, ++hessian_entry_fid)
                    Q_terms_2nd_order_pt[i][j] = dQ_invQ2_pt[i] * dQ_pt(j)(0) * two_invQ_pt
                                                 - d2Q_pt(hessian_entry_fid)(0) * invQ2_pt;
            } // end loop i

        } // end loop pt


        const int n_funcs_comp = bspline_elem.get_num_basis(comp);

        const auto &w_coefs = nrb_space->weight_func_table_[comp]->get_coefficients();
        Assert(nrb_space->get_num_basis(comp) == w_coefs.size(),
               ExcDimensionMismatch(nrb_space->get_num_basis(comp),w_coefs.size()));

        const auto offset = comp_offset[comp];
        for (int w_fn_id = 0 ; w_fn_id < n_funcs_comp ; ++w_fn_id, ++bsp_fn_id)
        {
            const auto &P_fn   =   P.get_function_view(bsp_fn_id);
            const auto &dP_fn  =  dP.get_function_view(bsp_fn_id);
            const auto &d2P_fn = d2P.get_function_view(bsp_fn_id);

            auto d2R_fn = D2_phi.get_function_view(bsp_fn_id);

            const Real w = w_coefs[bsp_local_to_patch[bsp_fn_id]-offset];

            for (int pt = 0 ; pt < n_pts ; ++pt)
            {
                auto &d2R_fn_pt = d2R_fn[pt];

                const auto   &P_fn_pt_comp = P_fn[pt](comp);
                const auto  &dP_fn_pt =  dP_fn[pt];
                const auto &d2P_fn_pt = d2P_fn[pt];

                const Real invQ_pt = invQ[pt];
                const auto &dQ_invQ2_pt = dQ_invQ2[pt];
                const auto &Q_terms_2nd_order_pt = Q_terms_2nd_order[pt];

                int hessian_entry_fid = 0;
                for (int i = 0 ; i < dim ; ++i)
                {
                    const auto &dP_fn_pt_i_comp = dP_fn_pt(i)(comp);
                    const auto &dQ_invQ2_pt_i = dQ_invQ2_pt[i];
                    const auto &Q_terms_2nd_order_pt_i = Q_terms_2nd_order_pt[i];

                    for (int j = 0 ; j < dim ; ++j, ++hessian_entry_fid)
                    {
                        const auto &dP_fn_pt_j_comp = dP_fn_pt(j)(comp);
                        const auto &dQ_invQ2_pt_j = dQ_invQ2_pt[j];
                        const auto &Q_terms_2nd_order_pt_i_j = Q_terms_2nd_order_pt_i[j];

                        auto &d2R_fn_pt_i_j_comp = d2R_fn_pt(hessian_entry_fid)(comp);
                        const auto &d2P_fn_pt_i_j_comp = d2P_fn_pt(hessian_entry_fid)(comp);

                        d2R_fn_pt_i_j_comp = (d2P_fn_pt_i_j_comp * invQ_pt
                                              - dP_fn_pt_i_comp  * dQ_invQ2_pt_j
                                              - dP_fn_pt_j_comp * dQ_invQ2_pt_i
                                              +  P_fn_pt_comp * Q_terms_2nd_order_pt_i_j) * w;
                    } // end loop j
                } // end loop i
            } // end loop pt
        } // end loop w_fn_id
    } // end loop comp
}


IGA_NAMESPACE_CLOSE

#include <igatools/basis_functions/nurbs_element_handler.inst>

#endif // #ifdef NURBS
