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

#include <igatools/basis_functions/nurbs_uniform_quad_cache.h>
#include <igatools/utils/multi_array_utils.h>

#include <algorithm>
using std::shared_ptr;

IGA_NAMESPACE_OPEN

namespace
{
ValueFlags
nurbs_to_bspline_flags(const ValueFlags nurbs_flags)
{
    int max_der_order = -1;
    ValueFlags bspline_flags = nurbs_flags;
    if (contains(nurbs_flags, ValueFlags::value))
    {
        max_der_order=std::max(max_der_order,0);
        bspline_flags |= ValueFlags::value;
    }

    if (contains(nurbs_flags, ValueFlags::face_value))
    {
        max_der_order=std::max(max_der_order,0);
        bspline_flags |= ValueFlags::face_value;
    }

    if (contains(nurbs_flags, ValueFlags::gradient))
    {
        max_der_order=std::max(max_der_order,1);
        bspline_flags |= ValueFlags::value |
                         ValueFlags::gradient;
    }

    if (contains(nurbs_flags, ValueFlags::face_gradient))
    {
        max_der_order=std::max(max_der_order,1);
        bspline_flags |= ValueFlags::face_value |
                         ValueFlags::face_gradient;
    }

    if (contains(nurbs_flags, ValueFlags::hessian))
    {
        max_der_order=std::max(max_der_order,2);
        bspline_flags |= ValueFlags::value |
                         ValueFlags::gradient |
                         ValueFlags::hessian;
    }

    if (contains(nurbs_flags, ValueFlags::face_hessian))
    {
        max_der_order=std::max(max_der_order,2);
        bspline_flags |= ValueFlags::face_value |
                         ValueFlags::face_gradient |
                         ValueFlags::face_hessian;
    }
    Assert(max_der_order>=0, ExcMessage("Not a right ValueFlag"));

    return bspline_flags;
}

};



template<int dim_, int range_ , int rank_>
NURBSUniformQuadCache<dim_, range_, rank_>::
NURBSUniformQuadCache(shared_ptr<const Space> space,
                      const ValueFlags flag,
                      const Quadrature<dim> &quad)
    :
    GridUniformQuadCache<dim_>(space->get_grid(), flag, quad),
    space_(space),
//    n_basis_(space_->get_num_all_element_basis()),
    flags_(flag),
    face_flags_(flag),
//    quad_(quad),
    bspline_uniform_quad_cache_(
        space->get_spline_space(),
        nurbs_to_bspline_flags(flag),
        quad)
{
    Assert(contains(flag, ValueFlags::none),
           ExcMessage("Nothing to reset"));
}

#if 0
template<int dim_, int range_ , int rank_>
BSplineUniformQuadCache<dim_, range_, rank_>::
BSplineUniformQuadCache(shared_ptr<const Space> space,
                        const ValueFlags flag,
                        const Quadrature<dim> &quad)
    :
    GridUniformQuadCache<dim_>(space->get_grid(), flag, quad),
    space_(space),
    n_basis_(space_->get_num_all_element_basis()),
    flags_(flag),
    face_flags_(flag),
    quad_(quad),
    splines1d_(space->get_grid()->get_num_intervals(),
               BasisValues(space->get_components_map()))
{


    // Compute the component offsets
    comp_offset_[0] = 0;
    for (int j = 1; j < Space::n_components; ++j)
        comp_offset_[j] = comp_offset_[j-1] + n_basis_.comp_dimension[j-1];


    const auto &n_inter = space->get_grid()->get_num_intervals();
    const auto &n_points = quad.get_num_points_direction();

    // Allocate space for the BasisValues1D
    for (int dir = 0 ; dir < dim ; ++dir)
    {
        const auto &n_pts = n_points[dir];
        for (int j = 0 ; j < n_inter[dir] ; ++j)
        {
            auto &splines1d = splines1d_.entry(dir, j);
            for (auto comp : splines1d.get_active_components_id())
                splines1d[comp].resize(n_derivatives, n_basis_[comp][dir], n_pts);
        }
    }

    /*
     * For each direction, interval and component we compute the 1D bspline
     * basis evaluate at the 1D component of the tensor product quadrature
     */
    const auto &degree      = space->get_degree();
    const auto &bezier_op   = space_->operators_;
    const auto &points      = quad_.get_points();
    const auto &lengths = this->lengths_;

    BasisValues bernstein_values(n_basis_.get_comp_map());

    for (int dir = 0 ; dir < dim ; ++dir)
    {
        // fill values and derivatives of the Bernstein's polynomials at
        // quad points in [0,1]
        for (auto comp : bernstein_values.get_active_components_id())
        {
            const int deg = degree[comp][dir];
            bernstein_values[comp].resize(n_derivatives, deg+1, n_points[dir]);
            const auto &pt_coords = points.get_data_direction(dir);
            for (int order = 0; order < n_derivatives; ++order)
                bernstein_values[comp].get_derivative(order) =
                    BernsteinBasis::derivative(order, deg, pt_coords);
        }

        const auto &inter_lengths = lengths.get_data_direction(dir);
        for (int j = 0 ; j < n_inter[dir] ; ++j)
        {
            auto &splines1d = splines1d_.entry(dir, j);
            for (auto comp : splines1d.get_active_components_id())
            {
                const auto &berns_values = bernstein_values[comp];
                auto &basis = splines1d[comp];
                const auto &oper = bezier_op.get_operator(comp,dir)[j];
                const Real one_div_size = 1.0 / inter_lengths[j];
                for (int order = 0; order < n_derivatives; ++order)
                {
                    const Real scale = std::pow(one_div_size, order);
                    const auto &b_values = berns_values.get_derivative(order);
                    basis.get_derivative(order) =
                        scale * prec_prod(oper, b_values);
                }
            }
        }

    }

//        for (auto comp : splines1d_.get_active_components_id())
//    {
//        auto &splines1d = splines1d_[comp];
//
//        {
//            const int num_intervals = n_intervals[jDim];
//            const int deg = degree[comp][jDim];
//            BasisValues1d bernstein_values(n_derivatives, deg+1, n_points[jDim]);
//
//            const auto &pt_coords = points.get_data_direction(jDim);
//
//            // fill values and derivatives of the Bernstein's polynomials at
//            // quad points in [0,1]
//            for (int order = 0; order < n_derivatives; ++order)
//                bernstein_values.get_derivative(order) =
//                    BernsteinBasis::derivative(order, deg, pt_coords);
//
//            const auto &bez_iComp_jDim = bezier_op.get_operator(comp,jDim);
//            const auto &lengths_jDim = lengths.get_data_direction(jDim);
//
//            // compute the one dimensional B-splines at quad point on the reference interval
//            for (int i = 0 ; i < num_intervals ; ++i)
//            {
//                const auto &M = bez_iComp_jDim[i];
//                const Real one_div_size = 1.0 / lengths_jDim[i];
//                BasisValues1d &basis = splines1d.entry(jDim,i);
//
//                for (int order = 0; order < n_derivatives; ++order)
//                {
//                    const Real scaling_coef = std::pow(one_div_size, order);
//                    basis.get_derivative(order) = scaling_coef * prec_prod(M, bernstein_values.get_derivative(order));
//                } //end loop order
//
//            } // end loop interval
//        } // end loop jDim
//    } // end loop comp
}
#endif


template<int dim_, int range_ , int rank_>
void
NURBSUniformQuadCache<dim_, range_, rank_>::
init_element_cache(ElementAccessor &elem)
{
    base_t::init_element_cache(elem);

    // init the element values for the cache of the BSplineElementAccessor
    bspline_uniform_quad_cache_.init_element_cache(elem.bspline_element_accessor_);
    const auto &quad = bspline_uniform_quad_cache_.get_quad();

    //-----------------------------------------------------------------------------
    auto &cache = elem.local_cache_;
    if (cache == nullptr)
    {
        using Cache = typename ElementAccessor::LocalCache;
        cache = shared_ptr<Cache>(new Cache);
    }


    auto n_basis = space_->get_num_all_element_basis();
    auto &elem_cache = cache->elem_values_;
    elem_cache.resize(flags_, quad, n_basis);

    auto &face_cache = cache->face_values_;
    for (auto f: base_t::faces)
        face_cache[f].resize(face_flags_, quad, n_basis, f);
    //-----------------------------------------------------------------------------
}



template<int dim_, int range_ , int rank_>
void
NURBSUniformQuadCache<dim_, range_, rank_>::
init_element_cache(ElementIterator &elem)
{
    init_element_cache(elem.get_accessor());
}

#if 0
template <int dim, int range, int rank>
void
BSplineUniformQuadCache<dim, range, rank>::
copy_to_inactive_components_values(const vector<Index> &inactive_comp,
                                   const std::array<Index, n_components> &active_map,
                                   ValueTable<Value> &D_phi) const
{
    const Size num_points = quad_.get_num_points_direction().flat_size();
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
            for (int qp = 0; qp < num_points; ++qp)
                inact_D_phi[qp](comp) = act_D_phi[qp](act_comp);
        }
    }
}



template <int dim, int range, int rank>
template <int order>
void
BSplineUniformQuadCache<dim, range, rank>::
copy_to_inactive_components(const vector<Index> &inactive_comp,
                            const std::array<Index, n_components> &active_map,
                            ValueTable<Derivative<order>> &D_phi) const
{
    const Size num_points = quad_.get_num_points_direction().flat_size();
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
            for (int qp = 0; qp < num_points; ++qp)
                for (int der = 0; der < n_ders; ++der)
                    inact_D_phi[qp](der)(comp) = act_D_phi[qp](der)(act_comp);
        }
    }
}




template <int dim, int range, int rank>
void
BSplineUniformQuadCache<dim, range, rank>::
evaluate_bspline_values(
    const ComponentContainer<TensorProductFunctionEvaluator<dim>> &elem_values,
    ValueTable<Value> &D_phi) const
{
    const auto n_points_direction = quad_.get_num_points_direction();
    const Size num_points = n_points_direction.flat_size();
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
            for (int point_id = 0; point_id < num_points; ++point_id)
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
BSplineUniformQuadCache<dim, range, rank>::
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

    const auto n_points_direction = quad_.get_num_points_direction();
    const Size num_points = n_points_direction.flat_size();
    Assert(D_phi.get_num_points() == num_points,
           ExcDimensionMismatch(D_phi.get_num_points(),num_points));


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
                auto const &der_tensor_id = univariate_order[der_id];
                for (int point_id = 0; point_id < num_points; ++point_id)
                {
                    auto const &pts  = values.points_flat_to_tensor(point_id);
                    auto &der = D_phi_i[point_id];
                    der(copy_indices[der_id][0])(comp) = values.evaluate(der_tensor_id, func, pts);
                    for (int k=1; k<copy_indices[der_id].size(); ++k)
                        der(copy_indices[der_id][k])(comp) = der(copy_indices[der_id][0])(comp);
                }
            }

        }
    } // end comp loop

    copy_to_inactive_components<order>(elem_values.get_inactive_components_id(),
                                       elem_values.get_comp_map(), D_phi);
}
#endif



template<int dim_, int range_ , int rank_>
void
NURBSUniformQuadCache<dim_, range_, rank_>::
fill_element_cache(ElementIterator &elem)
{
    fill_element_cache(elem.get_accessor());
}



template<int dim_, int range_ , int rank_>
void
NURBSUniformQuadCache<dim_, range_, rank_>::
fill_element_cache(ElementAccessor &elem)
{
    base_t::fill_element_cache(elem);
    bspline_uniform_quad_cache_.fill_element_cache(elem.bspline_element_accessor_);

    const auto &bspline_cache = elem.bspline_element_accessor_.get_elem_cache();
    auto &nurbs_cache = elem.get_elem_cache();


    const auto &elem_weights = elem.get_local_weights();  // NURBS weights local to the element
    const auto elem_basis_offset = elem.get_basis_offset();

    //--------------------------------------------------------------------------
    auto &flags_handler = nurbs_cache.flags_handler_;
    if (flags_handler.fill_values())
    {
        this->evaluate_nurbs_values(bspline_cache,
                                    elem_weights, elem_basis_offset,nurbs_cache.phi_);
        flags_handler.set_values_filled(true);
    }
    if (flags_handler.fill_gradients())
    {
        this->evaluate_nurbs_gradients(bspline_cache,
                                       elem_weights, elem_basis_offset, nurbs_cache.D1phi_);
        flags_handler.set_gradients_filled(true);
    }

    if (flags_handler.fill_hessians())
    {
        this->evaluate_nurbs_hessians(bspline_cache,
                                      elem_weights, elem_basis_offset, nurbs_cache.D2phi_);
        flags_handler.set_hessians_filled(true);
    }

    if (flags_handler.fill_divergences())
    {
        //TODO(pauletti, Sep 7, 2014): create a specialize exception
        Assert(flags_handler.gradients_filled(),
               ExcMessage("Divergence requires gradient to be filled."));

        auto D1  = nurbs_cache.D1phi_.begin();
        auto div = nurbs_cache.div_phi_.begin();
        auto end = nurbs_cache.D1phi_.end();
        for (; D1 != end; ++D1, ++div)
            *div = trace(*D1);

        flags_handler.set_divergences_filled(true);
    }

    //--------------------------------------------------------------------------

    nurbs_cache.set_filled(true);
}

template <int dim_, int range_, int rank_ >
void
NURBSUniformQuadCache<dim_, range_, rank_>::
evaluate_nurbs_values(
    const ElementCache &bspline_cache,
    const vector<Real> &weights,
    const ComponentContainer<int> &elem_basis_offset,
    ValueTable<Value> &D0_phi_hat) const
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
     * and
     *
     * P_i = w_i * N_i
     *
     *
     *
     *     _n_
     *     \
     * Q = /__  P_i
     *    i = 0
     *
     */

    Assert(bspline_cache.is_initialized(),ExcNotInitialized());
    const auto &bspline_values = bspline_cache.get_values();

    Assert(D0_phi_hat.get_num_functions() == weights.size(),
           ExcDimensionMismatch(D0_phi_hat.get_num_functions(), weights.size()));

    Assert(D0_phi_hat.get_num_functions() == bspline_values.get_num_functions(),
           ExcDimensionMismatch(D0_phi_hat.get_num_functions(), bspline_values.get_num_functions()));
    Assert(D0_phi_hat.get_num_points() == bspline_values.get_num_points(),
           ExcDimensionMismatch(D0_phi_hat.get_num_points(), bspline_values.get_num_points()));


    const auto &w_table = this->space_->weights_;
    const int num_points = D0_phi_hat.get_num_points();

    const auto spline_space = this->space_->get_spline_space();
    const auto num_basis_element =  spline_space->get_num_all_element_basis();

    for (int iComp : w_table.get_active_components_id())
    {
        const int num_basis_comp = num_basis_element.comp_dimension[iComp];

        vector<vector<Real>> P(num_basis_comp, vector<Real>(num_points));
        vector< Real > Q(num_points);

        for (int i = 0; i < num_basis_comp; ++i)
        {
            const int basis_flat_id = elem_basis_offset[iComp] + i;

            const auto &N_i = bspline_values.get_function_view(basis_flat_id);
            const Real w_i = weights[basis_flat_id];

            auto &P_i = P[i];

            for (int iPt = 0; iPt < num_points; iPt++)
            {
                P_i[iPt] = w_i * N_i[iPt](iComp);
                Q[iPt] += P_i[iPt];
            }
        }

        vector< Real >  invQ(num_points);
        for (int iPt = 0; iPt < num_points; ++iPt)
            invQ[iPt] = 1.0 / Q[iPt];

        for (int i = 0; i < num_basis_comp; i++)
        {
            const int basis_flat_id = elem_basis_offset[iComp] + i;
            const auto &P_i = P[i];

            for (int iPt = 0; iPt < num_points; ++iPt)
            {
                auto &R = D0_phi_hat.get_function_view(basis_flat_id)[iPt];
                R(iComp) = invQ[iPt] * P_i[iPt];
            }
        }
    }

//    const auto inactive_components_id = w_table.get_inactive_components_id();
//    const auto comp_map = w_table.get_comp_map();
    bspline_uniform_quad_cache_.copy_to_inactive_components_values(
        w_table.get_inactive_components_id(),
        w_table.get_comp_map(),
        D0_phi_hat);


//    Assert(false,ExcNotImplemented());
//    AssertThrow(false,ExcNotImplemented());

#if 0
    for (int comp : w_table.get_inactive_components_id())
    {
        const auto n_basis = this->bspline_element_accessor_.get_num_basis(comp);
        const Size offset = this->comp_offset_[comp];
        const Size act_offset = this->comp_offset_[w_table.active(comp)];

        for (Size basis_i = 0; basis_i < n_basis;  ++basis_i)
        {
            const auto values_phi_hat_copy_from = D0_phi_hat.get_function_view(act_offset+basis_i);
            auto values_phi_hat_copy_to = D0_phi_hat.get_function_view(offset+basis_i);

            for (int qp = 0; qp < num_points; ++qp)
                values_phi_hat_copy_to[qp](comp) = values_phi_hat_copy_from[qp](w_table.active(comp));
        }
    }
#endif
}



template <int dim_, int range_, int rank_ >
void
NURBSUniformQuadCache<dim_, range_, rank_>::
evaluate_nurbs_gradients(
    const ElementCache &bspline_cache,
    const vector<Real> &weights,
    const ComponentContainer<int> &elem_basis_offset,
    ValueTable<Derivative<1> > &D1_phi_hat) const
{
    /*
     * This function evaluates the derivative of the n+1 NURBS basis function R_0,...,R_n
     * from the set of BSpline basis function N_0,...,N_n
     * where the i-th NURBS basis function is defined as
     *
     *         P_i
     * R_i = -------
     *          Q
     *
     *
     *          dP_i       P_i * dQ
     * dR_i = -------  -  ------------
     *           Q            Q*Q
     *
     * and
     *
     * P_i = w_i * N_i
     *
     *
     * dP_i = w_i * dN_i
     *
     *
     *     _n_
     *     \
     * Q = /__  P_i
     *    i = 0
     *
     *
     *      _n_
     *      \
     * dQ = /__  dP_i
     *     i = 0
     */

    Assert(bspline_cache.is_initialized(),ExcNotInitialized());
    const auto &bspline_values = bspline_cache.get_values();
    const auto &bspline_gradients = bspline_cache.get_gradients();

    Assert(D1_phi_hat.get_num_functions() == weights.size(),
           ExcDimensionMismatch(D1_phi_hat.get_num_functions(), weights.size()));

    Assert(D1_phi_hat.get_num_functions() == bspline_values.get_num_functions(),
           ExcDimensionMismatch(D1_phi_hat.get_num_functions(), bspline_values.get_num_functions()));
    Assert(D1_phi_hat.get_num_points() == bspline_values.get_num_points(),
           ExcDimensionMismatch(D1_phi_hat.get_num_points(), bspline_values.get_num_points()));


    const auto &w_table = this->space_->weights_;
    const int num_points = D1_phi_hat.get_num_points();

    const auto spline_space = this->space_->get_spline_space();
    const auto num_basis_element =  spline_space->get_num_all_element_basis();


    using Grad = special_array<Real,dim>;
    for (int iComp : w_table.get_active_components_id())
    {
        const int num_basis_comp = num_basis_element.comp_dimension[iComp];

        vector<vector<Real> >  P(num_basis_comp, vector<Real>(num_points));
        vector<vector<Grad> > dP(num_basis_comp, vector<Grad>(num_points));

        vector<Real> Q(num_points);
        vector<Grad> dQ(num_points);

        for (int i = 0; i < num_basis_comp; ++i)
        {
            const int basis_flat_id = elem_basis_offset[iComp] + i;
            const auto  &N_i =    bspline_values.get_function_view(basis_flat_id);
            const auto &dN_i = bspline_gradients.get_function_view(basis_flat_id);
            const Real w_i = weights[basis_flat_id];

            auto &P_i =  P[i];
            auto &dP_i = dP[i];

            for (int iPt = 0; iPt < num_points; iPt++)
            {
                P_i[iPt] = w_i * N_i[iPt](iComp);
                Q[iPt] += P_i[iPt];

                auto &dP_i_iPt = dP_i[iPt];
                auto &dQ_iPt = dQ[iPt];
                for (int entry_flat_id = 0; entry_flat_id < dim; ++entry_flat_id)
                {
                    dP_i_iPt[entry_flat_id] = w_i * dN_i[iPt](entry_flat_id)(iComp);
                    dQ_iPt[entry_flat_id] += dP_i_iPt[entry_flat_id];
                }

            }
        }

        vector<Real>   invQ(num_points);
        vector<Real> invQ_2(num_points);
        vector<Grad>  dinvQ(num_points);
        for (int iPt = 0; iPt < num_points; ++iPt)
        {
            const Real invQ_tmp = 1.0 / Q[iPt];
            invQ  [iPt] = invQ_tmp;
            invQ_2[iPt] = invQ_tmp * invQ_tmp;

            for (int j = 0; j < dim; ++j)
                dinvQ[iPt][j] = - invQ_2[iPt] * dQ[iPt][j];
        }

        for (int i = 0; i < num_basis_comp; i++)
        {
            const int basis_flat_id = elem_basis_offset[iComp] + i;
            const auto &P_i =  P[i];
            const auto &dP_i = dP[i];

            for (int iPt = 0; iPt < num_points; ++iPt)
            {
                auto &dR = D1_phi_hat.get_function_view(basis_flat_id)[iPt];

                const Real invQ_tmp = invQ[iPt];
                const Real    P_tmp = P_i[iPt];

                const auto &dinvQ_tmp = dinvQ[iPt];
                const auto &dP_tmp = dP_i[iPt];

                for (int entry_flat_id = 0; entry_flat_id < dim; ++entry_flat_id)
                {
                    dR(entry_flat_id)(iComp) = invQ_tmp * dP_tmp[entry_flat_id] + dinvQ_tmp[entry_flat_id] * P_tmp;
                }
            }
        }
    }

    bspline_uniform_quad_cache_.template copy_to_inactive_components<1>(
        w_table.get_inactive_components_id(),
        w_table.get_comp_map(),
        D1_phi_hat);
#if 0
    for (int comp : w_table.get_inactive_components_id())
    {
        const auto n_basis = this->bspline_element_accessor_.get_num_basis(comp);
        const Size offset = this->comp_offset_[comp];
        const Size act_offset = this->comp_offset_[w_table.active(comp)];

        for (Size basis_i = 0; basis_i < n_basis;  ++basis_i)
        {
            const auto values_phi_hat_copy_from = D1_phi_hat.get_function_view(act_offset+basis_i);
            auto values_phi_hat_copy_to = D1_phi_hat.get_function_view(offset+basis_i);

            for (int qp = 0; qp < num_points; ++qp)
                values_phi_hat_copy_to[qp](comp) = values_phi_hat_copy_from[qp](w_table.active(comp));
        }
    }
#endif
}

template <int dim_, int range_, int rank_ >
void
NURBSUniformQuadCache<dim_, range_, rank_>::
evaluate_nurbs_hessians(
    const ElementCache &bspline_cache,
    const vector<Real> &weights,
    const ComponentContainer<int> &elem_basis_offset,
    ValueTable<Derivative<2> > &D2_phi_hat) const
{
    /*
     * This function evaluates the derivative of the n+1 NURBS basis function R_0,...,R_n
     * from the set of BSpline basis function N_0,...,N_n
     * where the i-th NURBS basis function is defined as
     *
     *         P_i
     * R_i = -------
     *          Q
     *
     *
     *          dP_i       P_i * dQ
     * dR_i = -------  -  ------------
     *           Q            Q*Q
     *
     *
     *           d2P_i     ( 2 * dP_i * dQ + P_i * d2Q )      2 * P_i * dQ * dQ
     * d2R_i = -------- - ------------------------------- + ---------------------
     *             Q                    Q*Q                        Q*Q*Q
     *
     * and
     *
     * P_i = w_i * N_i
     *
     *
     * dP_i = w_i * dN_i
     *
     *
     * d2P_i = w_i * d2N_i
     *
     *
     *     _n_
     *     \
     * Q = /__  P_i
     *    i = 0
     *
     *
     *      _n_
     *      \
     * dQ = /__  dP_i
     *     i = 0
     *
     *
     *       _n_
     *       \
     * d2Q = /__  d2P_i
     *      i = 0
     */

    Assert(bspline_cache.is_initialized(),ExcNotInitialized());
    const auto &bspline_values = bspline_cache.get_values();
    const auto &bspline_gradients = bspline_cache.get_gradients();
    const auto &bspline_hessians = bspline_cache.get_hessians();

    Assert(D2_phi_hat.get_num_functions() == weights.size(),
           ExcDimensionMismatch(D2_phi_hat.get_num_functions(), weights.size()));

    Assert(D2_phi_hat.get_num_functions() == bspline_values.get_num_functions(),
           ExcDimensionMismatch(D2_phi_hat.get_num_functions(), bspline_values.get_num_functions()));
    Assert(D2_phi_hat.get_num_points() == bspline_values.get_num_points(),
           ExcDimensionMismatch(D2_phi_hat.get_num_points(), bspline_values.get_num_points()));


    const auto &w_table = this->space_->weights_;
    const int num_points = D2_phi_hat.get_num_points();

    const auto spline_space = this->space_->get_spline_space();

    const auto num_basis_element =  spline_space->get_num_all_element_basis();


    using Grad = special_array<Real,dim> ;
    using Hess = special_array<special_array<Real,dim>,dim>;
    for (int iComp : w_table.get_active_components_id())
    {
        const int num_basis_comp = num_basis_element.comp_dimension[iComp];

        vector<vector<Real>>   P(num_basis_comp, vector<Real>(num_points));
        vector<vector<Grad>>  dP(num_basis_comp, vector<Grad>(num_points));
        vector<vector<Hess>> d2P(num_basis_comp, vector<Hess>(num_points));

        vector<Real>  Q(num_points);
        vector<Grad> dQ(num_points);
        vector<Hess> d2Q(num_points);

        for (int i = 0; i < num_basis_comp; ++i)
        {
            const int basis_flat_id = elem_basis_offset[iComp] + i;
            const auto   &N_i =    bspline_values.get_function_view(basis_flat_id);
            const auto  &dN_i = bspline_gradients.get_function_view(basis_flat_id);
            const auto &d2N_i =  bspline_hessians.get_function_view(basis_flat_id);

            const Real w_i = weights[basis_flat_id];
            auto   &P_i =   P[i];
            auto  &dP_i =  dP[i];
            auto &d2P_i = d2P[i];

            for (int iPt = 0; iPt < num_points; ++iPt)
            {
                P_i[iPt] = w_i * N_i[iPt](iComp);
                Q[iPt] += P_i[iPt];

                int hess_entry_flat_id = 0;
                for (int j = 0; j < dim; ++j)
                {
                    dP_i[iPt][j] = w_i * dN_i[iPt](j)(iComp);
                    dQ[iPt][j] += dP_i[iPt][j];

                    for (int k = 0; k < dim; ++k, ++hess_entry_flat_id)
                    {
                        d2P_i[iPt][j][k] = w_i * d2N_i[iPt](hess_entry_flat_id)(iComp);
                        d2Q[iPt][j][k] += d2P_i[iPt][j][k];
                    }
                }

            }
        }

        vector<Real> invQ(num_points);
        vector<Real> invQ_2(num_points);
        vector<Real> two_invQ_3(num_points);
        vector<Grad> dinvQ(num_points);
        vector<Hess> d2invQ(num_points);
        for (int iPt = 0; iPt < num_points; ++iPt)
        {
            const Real invQ_tmp = 1.0 / Q[iPt];
            invQ  [iPt] = invQ_tmp;
            invQ_2[iPt] = invQ_tmp * invQ_tmp;
            two_invQ_3[iPt] = 2.0 * invQ_tmp * invQ_tmp * invQ_tmp;

            for (int j = 0; j < dim; ++j)
            {
                dinvQ[iPt][j] = - invQ_2[iPt] * dQ[iPt][j];
                for (int k = 0; k < dim; ++k)
                {
                    d2invQ[iPt][j][k] = - invQ_2[iPt] * d2Q[iPt][j][k] + two_invQ_3[iPt] * dQ[iPt][j] * dQ[iPt][k];
                }
            }
        }

        for (int i = 0; i < num_basis_comp; i++)
        {
            const int basis_flat_id = elem_basis_offset[iComp] + i;

            const auto   &P_i =   P[i];
            const auto  &dP_i =  dP[i];
            const auto &d2P_i = d2P[i];

            for (int iPt = 0; iPt < num_points; iPt++)
            {
                auto &d2R = D2_phi_hat.get_function_view(basis_flat_id)[iPt];
                const Real invQ_tmp = invQ[iPt];
                const Real    P_tmp = P_i[iPt];

                const auto &dinvQ_tmp = dinvQ[iPt];
                const auto    &dP_tmp = dP_i[iPt];

                const auto &d2invQ_tmp = d2invQ[iPt];
                const auto    &d2P_tmp = d2P_i[iPt];

                int hess_entry_flat_id = 0;
                for (int j = 0; j < dim; ++j)
                {
                    for (int k = 0; k < dim; ++k, ++hess_entry_flat_id)
                    {
                        d2R(hess_entry_flat_id)(iComp) = d2invQ_tmp[j][k] *   P_tmp +
                                                         dinvQ_tmp[j]     *  dP_tmp[k] +
                                                         dinvQ_tmp[k]     *  dP_tmp[j] +
                                                         invQ_tmp         * d2P_tmp[j][k];
                    }
                }
            }
        }
        bspline_uniform_quad_cache_.template copy_to_inactive_components<2>(
            w_table.get_inactive_components_id(),
            w_table.get_comp_map(),
            D2_phi_hat);
#if 0
        for (int comp : w_table.get_inactive_components_id())
        {
            const auto n_basis = this->bspline_element_accessor_.get_num_basis(comp);
            const Size offset = this->comp_offset_[comp];
            const Size act_offset = this->comp_offset_[w_table.active(comp)];

            for (Size basis_i = 0; basis_i < n_basis;  ++basis_i)
            {
                const auto values_phi_hat_copy_from = D2_phi_hat.get_function_view(act_offset+basis_i);
                auto values_phi_hat_copy_to = D2_phi_hat.get_function_view(offset+basis_i);

                for (int qp = 0; qp < num_points; ++qp)
                    values_phi_hat_copy_to[qp](comp) = values_phi_hat_copy_from[qp](w_table.active(comp));
            }
        }
#endif
    }
}

template<int dim_, int range_ , int rank_>
void
NURBSUniformQuadCache<dim_, range_, rank_>::
print_info(LogStream &out) const
{
    Assert(false,ExcNotImplemented());
    AssertThrow(false,ExcNotImplemented());
#if 0

    out.begin_item("Grid Cache:");
    base_t::print_info(out);
    out.end_item();


    out.begin_item("One dimensional splines cache:");
    splines1d_.print_info(out);

    out.end_item();
#endif
}

IGA_NAMESPACE_CLOSE

#include <igatools/basis_functions/nurbs_uniform_quad_cache.inst>
