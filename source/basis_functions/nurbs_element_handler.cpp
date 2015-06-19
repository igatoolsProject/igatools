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
{
    const auto bsp_space = space->get_spline_space();
    bspline_handler_ = BSplineElementHandler<dim_,range_,rank_>::create(bsp_space);
}


template<int dim_, int range_ , int rank_>
auto
NURBSElementHandler<dim_, range_, rank_>::
create(std::shared_ptr<const Space> space) -> std::shared_ptr<base_t>
{
    return std::shared_ptr<self_t>(new self_t(space));
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


template<int dim_, int range_ , int rank_>
template<int sub_elem_dim>
void
NURBSElementHandler<dim_, range_, rank_>::
ResetDispatcher::
operator()(const Quadrature<sub_elem_dim> &quad)
{
    flags_[sub_elem_dim] = flag_;
}


template<int dim_, int range_ , int rank_>
void
NURBSElementHandler<dim_, range_, rank_>::
reset_selected_elements(
    const ValueFlags &flag,
    const eval_pts_variant &quad,
    const SafeSTLVector<int> &elements_flat_id)
{
    //--------------------------------------------------
    // resetting the Function for the bspline (numerator and weight function at the denominator)
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


    ValueFlags bspline_flag;
    if (max_deriv_order == 0)
        bspline_flag = ValueFlags::value;
    else if (max_deriv_order == 1)
        bspline_flag = ValueFlags::value | ValueFlags::gradient;
    else if (max_deriv_order == 2)
        bspline_flag = ValueFlags::value | ValueFlags::gradient | ValueFlags::hessian;
    else
        Assert(false,ExcMessage("Not a right value flag."));


    if (contains(flag, ValueFlags::point))
        bspline_flag |= ValueFlags::point;
    if (contains(flag, ValueFlags::w_measure))
        bspline_flag |= ValueFlags::w_measure;
    //--------------------------------------------------



    //--------------------------------------
    // resetting the BSplineElementHandler (for the numerator)
    Assert(bspline_handler_ != nullptr, ExcNullPtr());
    bspline_handler_->reset_selected_elements(bspline_flag, quad, elements_flat_id);
    //--------------------------------------


    //--------------------------------------------------
    const auto nrb_space = this->get_nurbs_space();
    nrb_space->weight_func_->reset_selected_elements(bspline_flag,quad,elements_flat_id);
    //--------------------------------------------------


    //--------------------------------------------------
    auto reset_dispatcher = ResetDispatcher(flag,flags_);
    boost::apply_visitor(reset_dispatcher, quad);
    //--------------------------------------------------
}


template<int dim_, int range_ , int rank_>
template<int sub_elem_dim>
void
NURBSElementHandler<dim_, range_, rank_>::
InitCacheDispatcher::
operator()(const Topology<sub_elem_dim> &sub_elem)
{
    grid_handler_.template init_cache<sub_elem_dim>(elem_.as_cartesian_grid_element_accessor());

    auto &cache = elem_.get_all_sub_elems_cache();
    if (cache == nullptr)
    {
        using VCache = typename NURBSElement<dim_,range_,rank_>::parent_t::Cache;

        using Cache = AllSubElementsCache<VCache>;

        cache = shared_ptr<Cache>(new Cache);
    }

    const auto n_basis = elem_.get_max_num_basis();//elem_.get_num_basis(DofProperties::active);
    const auto n_points = grid_handler_.template get_num_points<sub_elem_dim>();
    const auto flag = flags_[sub_elem_dim];

    for (auto &s_id: UnitElement<dim>::template elems_ids<sub_elem_dim>())
    {
        auto &s_cache = cache->template get_sub_elem_cache<sub_elem_dim>(s_id);
        s_cache.resize(flag, n_points, n_basis);
    }
}

template<int dim_, int range_ , int rank_>
void
NURBSElementHandler<dim_, range_, rank_>::
init_ref_elem_cache(RefElementAccessor &elem, const topology_variant &topology)
{
    Assert(!elem.get_space()->is_bspline(),ExcMessage("Not a NURBSElement."));

    auto &nrb_elem = dynamic_cast<NURBSElement<dim_,range_,rank_>&>(elem);
    auto &bsp_elem = nrb_elem.bspline_elem_;
    bspline_handler_->init_cache(bsp_elem,topology);

    const auto nrb_space = this->get_nurbs_space();
    nrb_space->weight_func_->init_cache(nrb_elem.weight_elem_,topology);


    //-------------------------------------
    auto init_cache_dispatcher = InitCacheDispatcher(
                                     const_cast<GridElementHandler<dim_> &>(bspline_handler_->get_grid_handler()),
                                     nrb_elem,flags_);
    boost::apply_visitor(init_cache_dispatcher,topology);
    //-------------------------------------
}



template<int dim_, int range_ , int rank_>
template<int sub_elem_dim>
void
NURBSElementHandler<dim_, range_, rank_>::
FillCacheDispatcher::
operator()(const Topology<sub_elem_dim> &sub_elem)
{
    grid_handler_.template fill_cache<sub_elem_dim>(
        nrb_elem_.as_cartesian_grid_element_accessor(),sub_elem_id_);

    Assert(nrb_elem_.all_sub_elems_cache_ != nullptr, ExcNullPtr());
    auto &sub_elem_cache = nrb_elem_.all_sub_elems_cache_->template get_sub_elem_cache<sub_elem_dim>(sub_elem_id_);

    const auto &bsp_elem = nrb_elem_.bspline_elem_;
    const auto &weight_elem = nrb_elem_.weight_elem_;

    if (sub_elem_cache.template status_fill<_Value>())
    {
        auto &values = sub_elem_cache.template get_data<_Value>();
        evaluate_nurbs_values_from_bspline(bsp_elem, weight_elem, values);
        sub_elem_cache.template set_status_filled<_Value>(true);
    }
    if (sub_elem_cache.template status_fill<_Gradient>())
    {
        auto &gradients = sub_elem_cache.template get_data<_Gradient>();
        evaluate_nurbs_gradients_from_bspline(bsp_elem, weight_elem, gradients);
        sub_elem_cache.template set_status_filled<_Gradient>(true);
    }
    if (sub_elem_cache.template status_fill<_Hessian>())
    {
        auto &hessians = sub_elem_cache.template get_data<_Hessian>();
        evaluate_nurbs_hessians_from_bspline(bsp_elem, weight_elem, hessians);
        sub_elem_cache.template set_status_filled<_Hessian>(true);
    }
    if (sub_elem_cache.template status_fill<_Divergence>())
    {
        eval_divergences_from_gradients(
            sub_elem_cache.template get_data<_Gradient>(),
            sub_elem_cache.template get_data<_Divergence>());
        sub_elem_cache.template set_status_filled<_Divergence>(true);
    }

    sub_elem_cache.set_filled(true);
}


template<int dim_, int range_ , int rank_>
void
NURBSElementHandler<dim_, range_, rank_>::
fill_ref_elem_cache(RefElementAccessor &elem, const topology_variant &topology, const int sub_elem_id)
{
    using Elem = NURBSElement<dim_,range_,rank_>;
    auto &nrb_elem = dynamic_cast<Elem &>(elem);
    auto &bsp_elem = nrb_elem.bspline_elem_;

    const auto nrb_space = this->get_nurbs_space();
    Assert(nrb_space == nrb_elem.get_nurbs_space(),
           ExcMessage("The element accessor and the element handler cannot have different spaces."));

    bspline_handler_->fill_cache(bsp_elem,topology,sub_elem_id);

    nrb_space->weight_func_->fill_cache(nrb_elem.weight_elem_,topology,sub_elem_id);


    //-----------------------------------------
    auto fill_cache_dispatcher =
        FillCacheDispatcher(bspline_handler_->get_grid_handler(),sub_elem_id,nrb_elem);
    boost::apply_visitor(fill_cache_dispatcher,topology);
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
    const WeightElem &weight_elem,
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

    const auto &P = bspline_elem.template get_basis<_Value,dim>(0,DofProperties::active);
    const auto n_pts = P.get_num_points();

    const auto bsp_local_to_patch = bspline_elem.get_local_to_patch(DofProperties::active);

    const auto nrb_space = nrb_elem_.get_nurbs_space();
    const auto comp_offset = nrb_space->get_dof_distribution()->get_dofs_offset();

    const auto &w_coefs = nrb_space->weight_func_->get_coefficients();

    int bsp_fn_id = 0;
    for (int comp = 0 ; comp < n_components ; ++comp)
    {
        Assert(nrb_space->get_num_basis(comp) == w_coefs.size(),
               ExcDimensionMismatch(nrb_space->get_num_basis(comp),w_coefs.size()));

        const auto &Q = weight_elem.template get_values<_Value,dim>(0);

        Assert(n_pts == Q.get_num_points(),
               ExcDimensionMismatch(n_pts,Q.get_num_points()));

        SafeSTLVector<Real> invQ(n_pts);
        for (int pt = 0 ; pt < n_pts ; ++pt)
            invQ[pt] = 1.0 / Q[pt](0);

        const int n_funcs_comp = bspline_elem.get_num_basis_comp(comp);


        const auto offset = comp_offset[comp];
        for (int w_fn_id = 0 ; w_fn_id < n_funcs_comp ; ++w_fn_id, ++bsp_fn_id)
        {
            const auto &P_fn = P.get_function_view(bsp_fn_id);

            auto R_fn = phi.get_function_view(bsp_fn_id);

            //TODO (pauletti, Mar 22, 2015): study the following line, why is it like this?
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
    const WeightElem &weight_elem,
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

    const auto &P  = bspline_elem.template get_basis<   _Value,dim>(0,DofProperties::active);
    const auto &dP = bspline_elem.template get_basis<_Gradient,dim>(0,DofProperties::active);

    const auto n_pts = P.get_num_points();

    const auto nrb_space = nrb_elem_.get_nurbs_space();
    const auto bsp_local_to_patch = bspline_elem.get_local_to_patch(DofProperties::active);
    const auto comp_offset = nrb_space->get_dof_distribution()->get_dofs_offset();

    const auto &w_coefs = nrb_space->weight_func_->get_coefficients();

    int bsp_fn_id = 0;
    for (int comp = 0 ; comp < n_components ; ++comp)
    {
        Assert(nrb_space->get_num_basis(comp) == w_coefs.size(),
               ExcDimensionMismatch(nrb_space->get_num_basis(comp),w_coefs.size()));

        const auto &Q  = weight_elem.template get_values<_Value,dim>(0);
        const auto &dQ = weight_elem.template get_values<_Gradient,dim>(0);

        Assert(n_pts == Q.get_num_points(),
               ExcDimensionMismatch(n_pts,Q.get_num_points()));


        SafeSTLVector<Real> invQ(n_pts);
        SafeSTLVector<SafeSTLArray<Real,dim>> dQ_invQ2(n_pts);
        for (int pt = 0 ; pt < n_pts ; ++pt)
        {
            invQ[pt] = 1.0 / Q[pt](0);

            const auto &dQ_pt = dQ[pt];
            auto &dQ_invQ2_pt = dQ_invQ2[pt];

            for (int i = 0 ; i < dim ; ++i)
                dQ_invQ2_pt[i] = invQ[pt] * invQ[pt] * dQ_pt(i)(0);
        }

        const int n_funcs_comp = bspline_elem.get_num_basis_comp(comp);

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
    const WeightElem &weight_elem,
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

    const auto &P   = bspline_elem.template get_basis<   _Value,dim>(0,DofProperties::active);
    const auto &dP  = bspline_elem.template get_basis<_Gradient,dim>(0,DofProperties::active);
    const auto &d2P = bspline_elem.template get_basis< _Hessian,dim>(0,DofProperties::active);

    const auto n_pts = P.get_num_points();

    const auto nrb_space = nrb_elem_.get_nurbs_space();
    const auto bsp_local_to_patch = bspline_elem.get_local_to_patch(DofProperties::active);
    const auto comp_offset = nrb_space->get_dof_distribution()->get_dofs_offset();

    const auto &w_coefs = nrb_space->weight_func_->get_coefficients();

    int bsp_fn_id = 0;
    for (int comp = 0 ; comp < n_components ; ++comp)
    {
        Assert(nrb_space->get_num_basis(comp) == w_coefs.size(),
               ExcDimensionMismatch(nrb_space->get_num_basis(comp),w_coefs.size()));


        const auto &Q   = weight_elem.template get_values<_Value,dim>(0);
        const auto &dQ  = weight_elem.template get_values<_Gradient,dim>(0);
        const auto &d2Q = weight_elem.template get_values<_Hessian,dim>(0);

        Assert(n_pts == Q.get_num_points(),
               ExcDimensionMismatch(n_pts,Q.get_num_points()));


        SafeSTLVector<Real> invQ(n_pts);
        SafeSTLVector<Real> invQ2(n_pts);
        SafeSTLVector<SafeSTLArray<Real,dim> > dQ_invQ2(n_pts);
        SafeSTLVector<SafeSTLArray<SafeSTLArray<Real,dim>,dim> > Q_terms_2nd_order(n_pts);
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


        const int n_funcs_comp = bspline_elem.get_num_basis_comp(comp);

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


#ifdef SERIALIZATION
template<int dim_, int range_ , int rank_>
template<class Archive>
void
NURBSElementHandler<dim_, range_, rank_>::
serialize(Archive &ar, const unsigned int version)
{
    ar &boost::serialization::make_nvp("NURBSElementHandler_base_t",
                                       boost::serialization::base_object<base_t>(*this));

    ar &boost::serialization::make_nvp("bspline_handler_",bspline_handler_);
    ar &boost::serialization::make_nvp("flags_",flags_);
}
///@}
#endif // SERIALIZATION


IGA_NAMESPACE_CLOSE

#include <igatools/basis_functions/nurbs_element_handler.inst>

#endif // #ifdef NURBS
