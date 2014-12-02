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

#include <igatools/basis_functions/nurbs_element_handler.h>
#include <igatools/basis_functions/nurbs_element.h>

using std::shared_ptr;

#ifdef NURBS

IGA_NAMESPACE_OPEN

template<int dim_, int range_ , int rank_>
NURBSElementHandler<dim_, range_, rank_>::
NURBSElementHandler(shared_ptr<const Space> space)
    :
    base_t(space->get_spline_space()),
    space_(space)
{}


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
template<int k>
void
NURBSElementHandler<dim_, range_, rank_>::
reset(const NewValueFlags flag,
      const Quadrature<k> &quad1)
{
    //--------------------------------------
    // resetting the BSplineElementHandler (for the numerator)
    base_t::template reset<k>(flag, quad1);
    //--------------------------------------


    //--------------------------------------------------
    // resetting the Function for the weight (for the denominator)
    int max_deriv_order = -1;
    if (contains(flag, NewValueFlags::point) ||
        contains(flag, NewValueFlags::value))
        max_deriv_order = 0;

    if (contains(flag, NewValueFlags::measure) ||
        contains(flag, NewValueFlags::w_measure) ||
        contains(flag, NewValueFlags::boundary_normal) ||
        contains(flag, NewValueFlags::outer_normal) ||
        contains(flag, NewValueFlags::gradient) ||
        contains(flag, NewValueFlags::inv_gradient))
        max_deriv_order = 1;


    if (contains(flag, NewValueFlags::curvature) ||
        contains(flag, NewValueFlags::hessian) ||
        contains(flag, NewValueFlags::inv_hessian))
        max_deriv_order = 2;


    NewValueFlags weight_flag;
    if (max_deriv_order == 0)
        weight_flag = NewValueFlags::value;
    else if (max_deriv_order == 1)
        weight_flag = NewValueFlags::value | NewValueFlags::gradient;
    else if (max_deriv_order == 2)
        weight_flag = NewValueFlags::value | NewValueFlags::gradient | NewValueFlags::hessian;
    else
        Assert(false,ExcMessage("Not a right value flag."));

    space_->weight_func_->reset(weight_flag,quad1);
    //--------------------------------------------------



    flags_[k] = flag;
}



template<int dim_, int range_ , int rank_>
template<int k>
void
NURBSElementHandler<dim_, range_, rank_>::
init_cache(ElementAccessor &elem)
{
    base_t::template init_cache<k>(elem.bspline_elem_);

    const auto topology = Int<k>();
    space_->weight_func_->init_cache(elem.weight_elem_,topology);



    auto &cache = elem.local_cache_;
    if (cache == nullptr)
    {
        using Cache = typename ElementAccessor::LocalCache;
        cache = shared_ptr<Cache>(new Cache);
    }

    const auto n_basis = elem.get_num_basis();

    for (auto &s_id: UnitElement<dim>::template elems_ids<k>())
    {
        auto &s_cache = cache->template get_value_cache<k>(s_id);
        const auto n_points = this->template get_num_points<k>();
        s_cache.resize(flags_[k], n_points, n_basis);
    }
}





template<int dim_, int range_ , int rank_>
void
NURBSElementHandler<dim_, range_, rank_>::
init_element_cache(ElementIterator &elem)
{
    init_cache<dim>(elem.get_accessor());
}


#if 0
template <int dim, int range, int rank>
void
NURBSElementHandler<dim, range, rank>::
copy_to_inactive_components_values(const vector<Index> &inactive_comp,
                                   const std::array<Index, n_components> &active_map,
                                   ValueTable<Value> &D_phi) const
{
    const Size n_points = D_phi.get_num_points();
    for (int comp : inactive_comp)
    {
        const auto act_comp = active_map[comp];
        const auto n_basis = n_basis_.comp_dimension[comp];
        const Size act_offset = comp_offset_[act_comp];
        const Size offset     = comp_offset_[comp];
        for (Size basis_i = 0; basis_i < n_basis;  ++basis_i)
        {
            const auto act_D_phi = D_phi.get_function_view(act_offset+basis_i);
            auto inact_D_phi = D_phi.get_function_view(offset+basis_i);
            for (int qp = 0; qp < n_points; ++qp)
                inact_D_phi[qp](comp) = act_D_phi[qp](act_comp);
        }
    }
}



template <int dim, int range, int rank>
template <int order>
void
NURBSElementHandler<dim, range, rank>::
copy_to_inactive_components(const vector<Index> &inactive_comp,
                            const std::array<Index, n_components> &active_map,
                            ValueTable<Derivative<order>> &D_phi) const
{
    const Size n_points = D_phi.get_num_points();
    const Size n_ders = Derivative<order>::size;
    for (int comp : inactive_comp)
    {
        const auto act_comp = active_map[comp];
        const auto n_basis = n_basis_.comp_dimension[comp];
        const Size act_offset = comp_offset_[act_comp];
        const Size offset     = comp_offset_[comp];
        for (Size basis_i = 0; basis_i < n_basis;  ++basis_i)
        {
            const auto act_D_phi = D_phi.get_function_view(act_offset+basis_i);
            auto     inact_D_phi = D_phi.get_function_view(offset+basis_i);
            for (int qp = 0; qp < n_points; ++qp)
                for (int der = 0; der < n_ders; ++der)
                    inact_D_phi[qp](der)(comp) = act_D_phi[qp](der)(act_comp);
        }
    }
}


template <int dim, int range, int rank>
void
NURBSElementHandler<dim, range, rank>::
evaluate_bspline_values(
    const ComponentContainer<TensorProductFunctionEvaluator<dim>> &elem_values,
    ValueTable<Value> &D_phi) const
{
    const Size n_points = D_phi.get_num_points();
    const TensorIndex<dim> der_tensor_id; // [0,0,..,0] tensor index
    for (int comp : elem_values.get_active_components_id())
    {
        auto &values = elem_values[comp];
        const int total_n_basis = n_basis_.comp_dimension[comp];
        const Size offset = comp_offset_[comp];

        for (int func_id = 0; func_id < total_n_basis; ++func_id)
        {
            auto D_phi_i = D_phi.get_function_view(offset + func_id);
            auto const &func = values.func_flat_to_tensor(func_id);
            for (int point_id = 0; point_id < n_points; ++point_id)
            {
                auto const &pts  = values.points_flat_to_tensor(point_id);
                D_phi_i[point_id](comp) = values.evaluate(der_tensor_id, func, pts);
            }
        } // end func_id loop
    } // end comp loop

    copy_to_inactive_components_values(elem_values.get_inactive_components_id(),
                                       elem_values.get_comp_map(), D_phi);
}



template <int dim, int range, int rank>
template <int order>
void
NURBSElementHandler<dim, range, rank>::
evaluate_bspline_derivatives(
    const ComponentContainer<TensorProductFunctionEvaluator<dim>> &elem_values,
    ValueTable<Derivative<order>> &D_phi) const
{
    /*
     * This code computes any order of derivatives for a multivariate
     * B-spline on the current element
     * We use the formula
     * \partial_(\alpha_1,...,\alpha_n) B(qp) = \Pi d^{\alpha_i} B_i(qp_i)
     */

    Assert(D_phi.size() > 0, ExcEmptyObject());
    //  Assert(D_phi.get_num_functions() == this->get_num_basis(),
    //           ExcDimensionMismatch(D_phi.get_num_functions(),this->get_num_basis()));
    const Size n_points = D_phi.get_num_points();


    TensorFunctionDerivativesSymmetry<dim,order> sym;
    const auto n_der =  TensorFunctionDerivativesSymmetry<dim,order>::num_entries_eval;

    const auto &univariate_order = sym.univariate_order ;
    const auto &copy_indices = sym.copy_indices;

    for (int comp : elem_values.get_active_components_id())
    {
        auto &values = elem_values[comp];
        const int total_n_basis = n_basis_.comp_dimension[comp];
        const Size offset = comp_offset_[comp];

        for (int func_id = 0; func_id < total_n_basis; ++func_id)
        {
            auto D_phi_i = D_phi.get_function_view(offset + func_id);
            auto const &func = values.func_flat_to_tensor(func_id);
            for (int der_id = 0; der_id<n_der; ++der_id)
            {
                const auto &copy_indices_der = copy_indices[der_id];
                const auto copy_indices_der_size = copy_indices_der.size();

                auto const &der_tensor_id = univariate_order[der_id];
                for (int point_id = 0; point_id < n_points; ++point_id)
                {
                    auto const &pts  = values.points_flat_to_tensor(point_id);
                    auto &der = D_phi_i[point_id];
                    der(copy_indices_der[0])(comp) = values.evaluate(der_tensor_id, func, pts);
                    for (int k=1; k<copy_indices_der_size; ++k)
                        der(copy_indices_der[k])(comp) = der(copy_indices_der[0])(comp);
                }
            }

        }
    } // end comp loop

    copy_to_inactive_components<order>(elem_values.get_inactive_components_id(),
                                       elem_values.get_comp_map(), D_phi);
}
#endif



template<int dim_, int range_ , int rank_>
template <int k>
void
NURBSElementHandler<dim_, range_, rank_>::
fill_cache(ElementAccessor &elem, const int j)
{
    base_t::template fill_cache<k>(elem.bspline_elem_, j);

    const auto topology = Int<k>();
    space_->weight_func_->fill_cache(elem.weight_elem_,j,topology);

    Assert(elem.local_cache_ != nullptr, ExcNullPtr());
    auto &cache = elem.local_cache_->template get_value_cache<k>(j);

    auto &flags = cache.flags_handler_;
    if (flags.fill_values())
    {
        auto &values = cache.template get_der<0>();
        evaluate_nurbs_values_from_bspline(elem.bspline_elem_,elem.weight_elem_,values);
        flags.set_values_filled(true);
    }
    if (flags.fill_gradients())
    {
        auto &gradients = cache.template get_der<1>();
        evaluate_nurbs_gradients_from_bspline(elem.bspline_elem_,elem.weight_elem_,gradients);
        flags.set_gradients_filled(true);
    }
    if (flags.fill_hessians())
    {
        auto &hessians = cache.template get_der<2>();
        evaluate_nurbs_hessians_from_bspline(elem.bspline_elem_,elem.weight_elem_,hessians);
        flags.set_hessians_filled(true);
    }

    cache.set_filled(true);
}



template<int dim_, int range_ , int rank_>
void
NURBSElementHandler<dim_, range_, rank_>::
fill_element_cache(ElementIterator &elem)
{
    fill_cache<dim>(elem.get_accessor(), 0);
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
evaluate_nurbs_values_from_bspline(
    const typename Space::SpSpace::ElementAccessor &bspline_elem,
    const typename Space::WeightFunction::ElementAccessor &weight_elem,
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
    const auto &Q = weight_elem.template get_values<0,dim>(0);

    /*
        LogStream out;
        out.begin_item("Denominator values");
        Q.print_info(out);
        out.end_item();
    //*/

    Assert(P.get_num_points() == Q.get_num_points(),
           ExcDimensionMismatch(P.get_num_points(),Q.get_num_points()));
    const auto n_pts = P.get_num_points();
    const auto n_funcs = P.get_num_functions();


    vector<Real> invQ(n_pts);
    for (int pt = 0 ; pt < n_pts ; ++pt)
        invQ[pt] = 1.0 / Q[pt](0);

    const auto local_to_patch = bspline_elem.get_local_to_patch();

    for (int fn = 0 ; fn < n_funcs ; ++fn)
    {
        const auto &P_fn = P.get_function_view(fn);

        auto R_fn = phi.get_function_view(fn);

        const Real w = space_->get_weight_coef_from_basis_id(local_to_patch[fn]);

        for (int pt = 0 ; pt < n_pts ; ++pt)
            R_fn[pt] = P_fn[pt] * invQ[pt] * w ;
    }
}

template<int dim_, int range_ , int rank_>
void
NURBSElementHandler<dim_, range_, rank_>::
evaluate_nurbs_gradients_from_bspline(
    const typename Space::SpSpace::ElementAccessor &bspline_elem,
    const typename Space::WeightFunction::ElementAccessor &weight_elem,
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

    const auto &Q  = weight_elem.template get_values<0,dim>(0);
    const auto &dQ = weight_elem.template get_values<1,dim>(0);

    Assert(P.get_num_points() == Q.get_num_points(),
           ExcDimensionMismatch(P.get_num_points(),Q.get_num_points()));
    const auto n_pts = P.get_num_points();
    const auto n_funcs = P.get_num_functions();

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


    const auto local_to_patch = bspline_elem.get_local_to_patch();

    for (int fn = 0 ; fn < n_funcs ; ++fn)
    {
        const auto &P_fn =  P.get_function_view(fn);
        const auto &dP_fn = dP.get_function_view(fn);

        auto dR_fn = D1_phi.get_function_view(fn);

        const Real w = space_->get_weight_coef_from_basis_id(local_to_patch[fn]);

        for (int pt = 0 ; pt < n_pts ; ++pt)
        {
            auto &dR_fn_pt = dR_fn[pt];

            const auto   &P_fn_pt =  P_fn[pt];
            const auto &dP_fn_pt = dP_fn[pt];

            const Real invQ_pt = invQ[pt];
            const auto &dQ_invQ2_pt = dQ_invQ2[pt];

            for (int i = 0 ; i < dim ; ++i)
            {
                const auto &dP_fn_pt_i = dP_fn_pt(i);
                const auto &dQ_invQ2_pt_i = dQ_invQ2_pt[i];

                auto &dR_fn_pt_i = dR_fn_pt(i);
                for (int comp = 0 ; comp < n_components ; ++comp)
                    dR_fn_pt_i(comp) = (dP_fn_pt_i(comp) * invQ_pt - P_fn_pt(comp) * dQ_invQ2_pt_i) * w;
            } // end loop i
        } // end loop pt
    } // end loop fn
}

template<int dim_, int range_ , int rank_>
void
NURBSElementHandler<dim_, range_, rank_>::
evaluate_nurbs_hessians_from_bspline(
    const typename Space::SpSpace::ElementAccessor &bspline_elem,
    const typename Space::WeightFunction::ElementAccessor &weight_elem,
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

    const auto &Q   = weight_elem.template get_values<0,dim>(0);
    const auto &dQ  = weight_elem.template get_values<1,dim>(0);
    const auto &d2Q = weight_elem.template get_values<2,dim>(0);

    Assert(P.get_num_points() == Q.get_num_points(),
           ExcDimensionMismatch(P.get_num_points(),Q.get_num_points()));
    const auto n_pts = P.get_num_points();
    const auto n_funcs = P.get_num_functions();

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


    const auto local_to_patch = bspline_elem.get_local_to_patch();

    for (int fn = 0 ; fn < n_funcs ; ++fn)
    {
        const auto &P_fn   =   P.get_function_view(fn);
        const auto &dP_fn  =  dP.get_function_view(fn);
        const auto &d2P_fn = d2P.get_function_view(fn);

        auto d2R_fn = D2_phi.get_function_view(fn);

        const Real w = space_->get_weight_coef_from_basis_id(local_to_patch[fn]);

        for (int pt = 0 ; pt < n_pts ; ++pt)
        {
            auto &d2R_fn_pt = d2R_fn[pt];

            const auto   &P_fn_pt =   P_fn[pt];
            const auto  &dP_fn_pt =  dP_fn[pt];
            const auto &d2P_fn_pt = d2P_fn[pt];

            const Real invQ_pt = invQ[pt];
            const auto &dQ_invQ2_pt = dQ_invQ2[pt];
            const auto &Q_terms_2nd_order_pt = Q_terms_2nd_order[pt];

            int hessian_entry_fid = 0;
            for (int i = 0 ; i < dim ; ++i)
            {

                const auto &dP_fn_pt_i = dP_fn_pt(i);
                const auto &dQ_invQ2_pt_i = dQ_invQ2_pt[i];
                const auto &Q_terms_2nd_order_pt_i = Q_terms_2nd_order_pt[i];

                for (int j = 0 ; j < dim ; ++j, ++hessian_entry_fid)
                {
                    const auto &dP_fn_pt_j = dP_fn_pt(j);
                    const auto &dQ_invQ2_pt_j = dQ_invQ2_pt[j];
                    const auto &Q_terms_2nd_order_pt_i_j = Q_terms_2nd_order_pt_i[j];

                    auto &d2R_fn_pt_i_j = d2R_fn_pt(hessian_entry_fid);
                    const auto &d2P_fn_pt_i_j = d2P_fn_pt(hessian_entry_fid);

                    for (int comp = 0 ; comp < n_components ; ++comp)
                        d2R_fn_pt_i_j(comp) =
                            (d2P_fn_pt_i_j(comp) * invQ_pt
                             - dP_fn_pt_i(comp) * dQ_invQ2_pt_j
                             - dP_fn_pt_j(comp) * dQ_invQ2_pt_i
                             +  P_fn_pt(comp) * Q_terms_2nd_order_pt_i_j) * w;
                } // end loop j
            } // end loop i
        } // end loop pt
    } // end loop fn
}


IGA_NAMESPACE_CLOSE

#include <igatools/basis_functions/nurbs_element_handler.inst>

#endif // #ifdef NURBS
