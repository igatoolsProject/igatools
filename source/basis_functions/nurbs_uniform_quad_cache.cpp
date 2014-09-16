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

    bspline_uniform_quad_cache_.copy_to_inactive_components_values(
        w_table.get_inactive_components_id(),
        w_table.get_comp_map(),
        D0_phi_hat);
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
    }
    bspline_uniform_quad_cache_.template copy_to_inactive_components<2>(
        w_table.get_inactive_components_id(),
        w_table.get_comp_map(),
        D2_phi_hat);
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
