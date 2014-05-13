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


#include <igatools/basis_functions/nurbs_element_accessor.h>
#include <igatools/basis_functions/nurbs_space.h>

using std::endl;

using std::array;
using std::vector;
using std::accumulate;

using std::make_shared;
using std::dynamic_pointer_cast;
IGA_NAMESPACE_OPEN


template< int dim, int range, int rank >
NURBSElementAccessor< dim, range, rank >::
NURBSElementAccessor(const std::shared_ptr<ContainerType> space,
                     const int elem_index)
    :
    SpaceElementAccessor<
    NURBSElementAccessor<dim,range,rank>,NURBSSpace<dim,range,rank>,dim,0,range,rank>(space,elem_index),
    bspline_element_accessor_(space->get_spline_space(), elem_index)
{}



template <int dim, int range, int rank>
void
NURBSElementAccessor<dim, range, rank>::
init_values(const ValueFlags fill_flag,
            const Quadrature<dim> &quad)
{
    Assert(contains(fill_flag, ValueFlags::none),
           ExcMessage("Nothing to reset"));

    int max_der_order = -1;
    ValueFlags fill_flag_bspline = fill_flag;
    if (contains(fill_flag, ValueFlags::value))
    {
        max_der_order=std::max(max_der_order,0);
        fill_flag_bspline |= ValueFlags::value;
    }

    if (contains(fill_flag, ValueFlags::face_value))
    {
        max_der_order=std::max(max_der_order,0);
        fill_flag_bspline |= ValueFlags::face_value;
    }

    if (contains(fill_flag, ValueFlags::gradient))
    {
        max_der_order=std::max(max_der_order,1);
        fill_flag_bspline |= ValueFlags::value |
                             ValueFlags::gradient;
    }

    if (contains(fill_flag, ValueFlags::face_gradient))
    {
        max_der_order=std::max(max_der_order,1);
        fill_flag_bspline |= ValueFlags::face_value |
                             ValueFlags::face_gradient;
    }

    if (contains(fill_flag, ValueFlags::hessian))
    {
        max_der_order=std::max(max_der_order,2);
        fill_flag_bspline |= ValueFlags::value |
                             ValueFlags::gradient |
                             ValueFlags::hessian;
    }

    if (contains(fill_flag, ValueFlags::face_hessian))
    {
        max_der_order=std::max(max_der_order,2);
        fill_flag_bspline |= ValueFlags::face_value |
                             ValueFlags::face_gradient |
                             ValueFlags::face_hessian;
    }

    Assert(max_der_order>=0, ExcMessage("Not a right ValueFlag"));

    // init the element values for the cache of the BSplineElementAccessor
    bspline_element_accessor_.init_values(fill_flag_bspline,quad);


    BasisElemValueFlagsHandler elem_flags_handler(fill_flag);
    BasisFaceValueFlagsHandler face_flags_handler(fill_flag);

    // reset the element values for the cache of the NURBSElementAccessor
    this->elem_values_.reset(elem_flags_handler, this->n_basis_direction_, quad);
//    this->elem_values_.reset(*(this->space_),elem_flags_handler,quad);

    Index face_id = 0;
    for (auto& face_value : this->face_values_)
        face_value.reset(face_id++, face_flags_handler, this->n_basis_direction_, quad);
//      face_value.reset(face_id++, *(this->space_), face_flags_handler, quad);
}



template <int dim, int range, int rank>
void
NURBSElementAccessor<dim, range, rank>::
init_face_values(const Index face_id,
                 const ValueFlags fill_flag,
                 const Quadrature<dim-1> &quad)
{
    AssertThrow(false,ExcNotImplemented());
}




template <int dim, int range, int rank  >
void
NURBSElementAccessor<dim, range, rank>::
evaluate_nurbs_values(
    const typename BSplineElementAccessor<dim,range,rank>::ValuesCache &bspline_cache,
    ValueTable<Value> &D0_phi_hat) const
{
    Assert(bspline_cache.is_initialized(),ExcNotInitialized());
    Assert(D0_phi_hat.get_num_functions() == this->get_num_basis(),
           ExcDimensionMismatch(D0_phi_hat.get_num_functions(), this->get_num_basis()));

    const int num_points = D0_phi_hat.get_num_points();

    {
        // here we treat the pure NURBS case
        using space_t = ContainerType;

        typedef Real ValueRange1_t;

        const vector< Real > &weights = this->get_weights();

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

        //----------------------------------------------------------------------------------------------
        const auto &bspline_values = bspline_cache.get_values();
        //----------------------------------------------------------------------------------------------

        if ((this->space_)->is_range_homogeneous() == false)
        {
            //------------------------------------------------------------------------------------------
            int dof_offset = 0;
            for (int iComp = 0; iComp < space_t::n_components; ++iComp)
            {
                const int num_basis_comp = this->get_num_basis(iComp);

                vector< vector<ValueRange1_t> > P(num_basis_comp, vector<Real>(num_points));

                vector< ValueRange1_t > Q(num_points);
                for (int i = 0; i < num_basis_comp; ++i)
                {
                    const int basis_flat_id = dof_offset + i;

                    const auto &N_i = bspline_values.get_function_view(basis_flat_id);
                    const Real w_i = weights[basis_flat_id];

                    auto &P_i = P[i];

                    for (int iPt = 0; iPt < num_points; iPt++)
                    {
                        P_i[iPt] = w_i * N_i[iPt](iComp);

                        Q[iPt] += P_i[iPt];
                    }
                }

                vector< ValueRange1_t >  invQ(num_points);
                for (int iPt = 0; iPt < num_points; ++iPt)
                    invQ[iPt] = 1.0 / Q[iPt];

                for (int i = 0; i < num_basis_comp; i++)
                {
                    const int basis_flat_id = dof_offset + i;

                    const auto &P_i = P[i];

                    for (int iPt = 0; iPt < num_points; ++iPt)
                    {
                        auto &R = D0_phi_hat.get_function_view(basis_flat_id)[iPt];

                        R(iComp) = invQ[iPt] * P_i[iPt];
                    }
                }
                dof_offset += num_basis_comp;

            } // end iComp loop
            //------------------------------------------------------------------------------------------
        }
        else // space_->homogeneous_range_ == true
        {
            //------------------------------------------------------------------------------------------
            const int num_basis_comp = this->get_num_basis(0);

            for (int iComp = 0; iComp < space_t::n_components; ++iComp)
            {
                Assert(this->get_num_basis(iComp) == num_basis_comp,
                       ExcDimensionMismatch(this->get_num_basis(iComp), num_basis_comp));
            }

            vector< vector<ValueRange1_t> > P(num_basis_comp, vector<Real>(num_points));
            vector< ValueRange1_t > Q(num_points);

            for (int i = 0; i < num_basis_comp; ++i)
            {
                const auto &N_i = bspline_values.get_function_view(i);
                const Real w_i = weights[i];

                auto &P_i = P[i];

                for (int iPt = 0; iPt < num_points; iPt++)
                {
                    P_i[iPt] = w_i * N_i[iPt](0);

                    Q[iPt] += P_i[iPt];
                }
            }

            vector< ValueRange1_t > invQ(num_points);
            for (int iPt = 0; iPt < num_points; ++iPt)
                invQ[iPt] = 1.0 / Q[iPt];


            for (int i = 0; i < num_basis_comp; i++)
            {
                const auto &P_i = P[i];

                for (int iPt = 0; iPt < num_points; ++iPt)
                {
                    const ValueRange1_t tmp_R = invQ[iPt] * P_i[iPt];

                    for (int iComp = 0; iComp < space_t::n_components; ++iComp)
                    {
                        const int basis_flat_id = i + iComp * num_basis_comp;

                        auto &R = D0_phi_hat.get_function_view(basis_flat_id)[iPt];
                        R(iComp) = tmp_R;
                    }
                }
            }
        }
    }
}


template <int dim, int range, int rank  >
void
NURBSElementAccessor< dim, range, rank >::
evaluate_nurbs_gradients(
    const typename BSplineElementAccessor<dim,range,rank>::ValuesCache &bspline_cache,
    ValueTable< Derivatives< dim, range, rank, 1 > > &D1_phi_hat) const
{
    Assert(bspline_cache.is_initialized(),ExcNotInitialized());
    Assert(D1_phi_hat.get_num_functions() == this->get_num_basis(),
           ExcDimensionMismatch(D1_phi_hat.get_num_functions(), this->get_num_basis()));

    const int num_points = D1_phi_hat.get_num_points();

    {
        // here we treat the pure NURBS case
        using space_t = ContainerType;

        typedef Real ValueRange1_t;
        typedef array<Real,dim> GradientRange1_t;


        const vector< Real > &weights = this->get_weights();

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

        //----------------------------------------------------------------------------------------------
        const auto &bspline_values = bspline_cache.get_values();
        const auto &bspline_gradients = bspline_cache.get_gradients();
        //----------------------------------------------------------------------------------------------


        if ((this->space_)->is_range_homogeneous() == false)
        {
            //------------------------------------------------------------------------------------------
            int dof_offset = 0;
            for (int iComp = 0; iComp < space_t::n_components; ++iComp)
            {
                const int num_basis_comp = this->get_num_basis(iComp);

                vector< vector< ValueRange1_t > > P(num_basis_comp, vector< ValueRange1_t >(num_points));
                vector< vector< GradientRange1_t > > dP(num_basis_comp, vector< GradientRange1_t >(num_points));

                vector<    ValueRange1_t >  Q(num_points);
                vector< GradientRange1_t > dQ(num_points);

                for (int i = 0; i < num_basis_comp; ++i)
                {
                    const int basis_flat_id = dof_offset + i;

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

                vector<    ValueRange1_t >   invQ(num_points);
                vector<    ValueRange1_t > invQ_2(num_points);
                vector< GradientRange1_t >  dinvQ(num_points);
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
                    const int basis_flat_id = dof_offset + i;

                    const auto &P_i =  P[i];
                    const auto &dP_i = dP[i];

                    for (int iPt = 0; iPt < num_points; iPt++)
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

                dof_offset += num_basis_comp;

            } // end iComp loop
            //------------------------------------------------------------------------------------------
        } // space_->homogeneous_range_ == false
        else // space_->homogeneous_range_ == true
        {
            //------------------------------------------------------------------------------------------
            const int num_basis_comp = this->get_num_basis(0);

            for (int iComp = 0; iComp < space_t::n_components; ++iComp)
            {
                Assert(this->get_num_basis(iComp) == num_basis_comp,
                       ExcDimensionMismatch(this->get_num_basis(iComp), num_basis_comp));
            }

            vector< vector< ValueRange1_t > > P(num_basis_comp, vector< ValueRange1_t >(num_points));
            vector< vector< GradientRange1_t > > dP(num_basis_comp, vector< GradientRange1_t >(num_points));

            vector<    ValueRange1_t >  Q(num_points);
            vector< GradientRange1_t > dQ(num_points);

            for (int i = 0; i < num_basis_comp; ++i)
            {
                const auto  &N_i =    bspline_values.get_function_view(i);
                const auto &dN_i = bspline_gradients.get_function_view(i);

                const Real w_i = weights[i];

                auto &P_i =  P[i];
                auto &dP_i = dP[i];

                for (int iPt = 0; iPt < num_points; iPt++)
                {
                    P_i[iPt] = w_i * N_i[iPt](0);

                    Q[iPt] += P_i[iPt];

                    auto &dP_i_iPt = dP_i[iPt];
                    auto &dQ_iPt = dQ[iPt];
                    for (int entry_flat_id = 0; entry_flat_id < dim; ++entry_flat_id)
                    {
                        dP_i_iPt[entry_flat_id] = w_i * dN_i[iPt](entry_flat_id)(0);

                        dQ_iPt[entry_flat_id] += dP_i_iPt[entry_flat_id];
                    }

                }
            }

            vector<    ValueRange1_t >   invQ(num_points);
            vector<    ValueRange1_t > invQ_2(num_points);
            vector< GradientRange1_t >  dinvQ(num_points);
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
                const auto &P_i =  P[i];
                const auto &dP_i = dP[i];

                for (int iPt = 0; iPt < num_points; iPt++)
                {
                    const Real invQ_tmp = invQ[iPt];
                    const Real    P_tmp = P_i[iPt];

                    const auto &dinvQ_tmp = dinvQ[iPt];
                    const auto &dP_tmp = dP_i[iPt];

                    for (int entry_flat_id = 0; entry_flat_id < dim; ++entry_flat_id)
                    {
                        const auto dR_tmp = invQ_tmp * dP_tmp[entry_flat_id] + dinvQ_tmp[entry_flat_id] * P_tmp;

                        for (int iComp = 0; iComp < space_t::n_components; ++iComp)
                        {
                            auto &dR = D1_phi_hat.get_function_view(num_basis_comp * iComp + i)[iPt];
                            dR(entry_flat_id)(iComp) = dR_tmp;
                        } // end iComp loop
                    }
                }
            }

            //------------------------------------------------------------------------------------------
        }

    }
}


template <int dim, int range, int rank  >
void
NURBSElementAccessor< dim, range, rank >::
evaluate_nurbs_hessians(
    const typename BSplineElementAccessor<dim,range,rank>::ValuesCache &bspline_cache,
    ValueTable< Derivatives< dim, range, rank, 2 > > &D2_phi_hat) const
{
    Assert(bspline_cache.is_initialized(),ExcNotInitialized());
    Assert(D2_phi_hat.get_num_functions() == this->get_num_basis(),
           ExcDimensionMismatch(D2_phi_hat.get_num_functions(), this->get_num_basis()));

    const int num_points = D2_phi_hat.get_num_points();

    {
        // here we treat the pure NURBS case

        using space_t = ContainerType;

        typedef Real ValueRange1_t;
        typedef array<Real,dim> GradientRange1_t;
        typedef array<array<Real,dim>,dim> HessianRange1_t;


        const vector< Real > &weights = this->get_weights();


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
         //*/

        //----------------------------------------------------------------------------------------------
        const auto &bspline_values = bspline_cache.get_values();
        const auto &bspline_gradients = bspline_cache.get_gradients();
        const auto &bspline_hessians = bspline_cache.get_hessians();
        //----------------------------------------------------------------------------------------------


        if ((this->space_)->is_range_homogeneous() == false)
        {
            //------------------------------------------------------------------------------------------
            int dof_offset = 0;
            for (int iComp = 0; iComp < space_t::n_components; ++iComp)
            {
                const int num_basis_comp = this->get_num_basis(iComp);

                vector< vector< ValueRange1_t > >     P(num_basis_comp, vector< ValueRange1_t >(num_points));
                vector< vector< GradientRange1_t > > dP(num_basis_comp, vector< GradientRange1_t >(num_points));
                vector< vector< HessianRange1_t > > d2P(num_basis_comp, vector< HessianRange1_t >(num_points));

                vector<    ValueRange1_t >  Q(num_points);
                vector< GradientRange1_t > dQ(num_points);
                vector< HessianRange1_t > d2Q(num_points);

                for (int iFn = 0; iFn < num_basis_comp; ++iFn)
                {
                    const int basis_flat_id = dof_offset + iFn;

                    const auto   &N_i =    bspline_values.get_function_view(basis_flat_id);
                    const auto  &dN_i = bspline_gradients.get_function_view(basis_flat_id);
                    const auto &d2N_i =  bspline_hessians.get_function_view(basis_flat_id);

                    const Real w_i = weights[basis_flat_id];

                    auto   &P_i =   P[iFn];
                    auto  &dP_i =  dP[iFn];
                    auto &d2P_i = d2P[iFn];

                    for (int iPt = 0; iPt < num_points; iPt++)
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

                vector<    ValueRange1_t >   invQ(num_points);
                vector<    ValueRange1_t > invQ_2(num_points);
                vector<    ValueRange1_t > two_invQ_3(num_points);
                vector< GradientRange1_t >  dinvQ(num_points);
                vector<  HessianRange1_t > d2invQ(num_points);
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

                for (int iFn = 0; iFn < num_basis_comp; ++iFn)
                {
                    const int basis_flat_id = dof_offset + iFn;

                    const auto   &P_i =   P[iFn];
                    const auto  &dP_i =  dP[iFn];
                    const auto &d2P_i = d2P[iFn];

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

                dof_offset += num_basis_comp;

            } // end iComp loop
            //------------------------------------------------------------------------------------------
        } // space_->homogeneous_range_ == false
        else // space_->homogeneous_range_ == true
        {
            const int num_basis_comp = this->get_num_basis(0);

            for (int iComp = 0; iComp < space_t::n_components; ++iComp)
            {
                Assert(this->get_num_basis(iComp) == num_basis_comp,
                       ExcDimensionMismatch(this->get_num_basis(iComp), num_basis_comp));
            }

            //------------------------------------------------------------------------------------------

            vector< vector< ValueRange1_t > > P(num_basis_comp, vector< ValueRange1_t >(num_points));
            vector< vector< GradientRange1_t > > dP(num_basis_comp, vector< GradientRange1_t >(num_points));
            vector< vector< HessianRange1_t > > d2P(num_basis_comp, vector< HessianRange1_t >(num_points));

            vector<    ValueRange1_t >  Q(num_points);
            vector< GradientRange1_t > dQ(num_points);
            vector< HessianRange1_t > d2Q(num_points);

            for (int iFn = 0; iFn < num_basis_comp; ++iFn)
            {
                const auto   &N_i =    bspline_values.get_function_view(iFn);
                const auto  &dN_i = bspline_gradients.get_function_view(iFn);
                const auto &d2N_i =  bspline_hessians.get_function_view(iFn);

                const Real w_i = weights[iFn];

                auto   &P_i =   P[iFn];
                auto  &dP_i =  dP[iFn];
                auto &d2P_i = d2P[iFn];

                for (int iPt = 0; iPt < num_points; iPt++)
                {
                    P_i[iPt] = w_i * N_i[iPt](0);
                    Q[iPt] += P_i[iPt];

                    int hess_entry_flat_id = 0;
                    for (int j = 0; j < dim; ++j)
                    {
                        dP_i[iPt][j] = w_i * dN_i[iPt](j)(0);
                        dQ[iPt][j] += dP_i[iPt][j];

                        for (int k = 0; k < dim; ++k, ++hess_entry_flat_id)
                        {
                            d2P_i[iPt][j][k] = w_i * d2N_i[iPt](hess_entry_flat_id)(0);
                            d2Q[iPt][j][k] += d2P_i[iPt][j][k];
                        }
                    }
                } // end loop iPt
            } // end loop i

            vector<    ValueRange1_t >   invQ(num_points);
            vector<    ValueRange1_t > invQ_2(num_points);
            vector<    ValueRange1_t > two_invQ_3(num_points);
            vector< GradientRange1_t >  dinvQ(num_points);
            vector<  HessianRange1_t > d2invQ(num_points);
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

            for (int iFn = 0; iFn < num_basis_comp; ++iFn)
            {
                const auto   &P_i =   P[iFn];
                const auto  &dP_i =  dP[iFn];
                const auto &d2P_i = d2P[iFn];

                for (int iPt = 0; iPt < num_points; iPt++)
                {

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
                            const auto d2R_tmp = d2invQ_tmp[j][k] *   P_tmp +
                                                 dinvQ_tmp[j]    *  dP_tmp[k] +
                                                 dinvQ_tmp[k]    *  dP_tmp[j] +
                                                 invQ_tmp        * d2P_tmp[j][k];

                            for (int iComp = 0; iComp < space_t::n_components; ++iComp)
                            {
                                auto &d2R = D2_phi_hat.get_function_view(iFn + iComp*num_basis_comp)[iPt];
                                d2R(hess_entry_flat_id)(iComp) = d2R_tmp;
                            } // end iComp loop

                        }
                    }
                }
            }
            //------------------------------------------------------------------------------------------
        } // space_->homogeneous_range_ == true
    }
}



template <int dim, int range, int rank >
void
NURBSElementAccessor< dim, range, rank >::
fill_values()
{
    Assert(this->elem_values_.is_initialized(),ExcNotInitialized());

    // fills the cache of the BSplineElementAccessor
    bspline_element_accessor_.fill_values();
    const auto &bspline_elem_cache = bspline_element_accessor_.get_values_cache();

    if (this->elem_values_.flags_handler_.fill_values())
    {
        evaluate_nurbs_values(
            bspline_elem_cache,
            this->elem_values_.phi_hat_);

        this->elem_values_.flags_handler_.set_values_filled(true);
    }

    if (this->elem_values_.flags_handler_.fill_gradients())
    {
        evaluate_nurbs_gradients(
            bspline_elem_cache,
            this->elem_values_.D1phi_hat_);

        this->elem_values_.flags_handler_.set_gradients_filled(true);
    }

    if (this->elem_values_.flags_handler_.fill_hessians())
    {
        evaluate_nurbs_hessians(
            bspline_elem_cache,
            this->elem_values_.D2phi_hat_);

        this->elem_values_.flags_handler_.set_hessians_filled(true);
    }

    this->elem_values_.set_filled(true);
}



template <int dim, int range, int rank >
void
NURBSElementAccessor< dim, range, rank >::
fill_face_values(const Index face_id)
{
    Assert(face_id < n_faces && face_id >= 0, ExcIndexRange(face_id,0,n_faces));

    auto &face_value = this->face_values_[face_id];

    Assert(face_value.is_initialized(),ExcNotInitialized());

    // fills the cache of the BSplineElementAccessor
    bspline_element_accessor_.fill_face_values(face_id);
    const auto &bspline_face_cache =
        bspline_element_accessor_.get_values_cache(FaceTopology<dim>(face_id));


    if (face_value.flags_handler_.fill_values())
    {
        evaluate_nurbs_values(
            bspline_face_cache,
            face_value.phi_hat_);

        face_value.flags_handler_.set_values_filled(true);
    }

    if (face_value.flags_handler_.fill_gradients())
    {
        evaluate_nurbs_gradients(
            bspline_face_cache,
            face_value.D1phi_hat_);

        face_value.flags_handler_.set_gradients_filled(true);
    }

    if (face_value.flags_handler_.fill_hessians())
    {
        evaluate_nurbs_hessians(
            bspline_face_cache,
            face_value.D2phi_hat_);

        face_value.flags_handler_.set_hessians_filled(true);
    }
    face_value.set_filled(true);
}




template <int dim, int range, int rank >
vector<Real>
NURBSElementAccessor< dim, range, rank >::
get_weights() const
{
    using space_t = BSplineSpace<dim,range,rank>;

    //---------------------------
    // here we compute the offset of the dofs relative to the components of the space
    array<Size,space_t::n_components+1> dofs_offset_comp;
    dofs_offset_comp[0] = 0;
    for (int comp = 0; comp < space_t::n_components; ++comp)
        dofs_offset_comp[comp+1] = dofs_offset_comp[comp] +
                                   (this->space_)->get_num_basis(comp);
    //---------------------------


    vector<Real> weights_element;

    const auto local_to_global = this->get_local_to_global();

    for (uint global_id : local_to_global)
    {
        Index comp_id = 0; // component id of the global index
        Index  dof_id = 0; // flat index of the global index relative to the component
        for (comp_id = 0; comp_id < space_t::n_components; ++comp_id)
        {
            if (global_id < dofs_offset_comp[comp_id+1])
            {
                dof_id = global_id - dofs_offset_comp[comp_id];
                break;
            }
        }

        weights_element.emplace_back((this->space_)->weights_(comp_id)(dof_id));
    }

    return weights_element;
}




template <int dim, int range, int rank>
bool
NURBSElementAccessor<dim, range, rank>::
operator==(const NURBSElementAccessor<dim,range,rank> &a) const
{
    const bool is_same_grid_element_accessor =
        (this->as_cartesian_grid_element_accessor() == a.as_cartesian_grid_element_accessor());
    const bool is_same_bspline_element_accessor =
        (this->bspline_element_accessor_ == a.bspline_element_accessor_);
    return (is_same_grid_element_accessor && is_same_bspline_element_accessor);
}

template <int dim, int range, int rank>
bool
NURBSElementAccessor<dim, range, rank>::
operator!=(const NURBSElementAccessor<dim,range,rank> &a) const
{
    return !((*this) == a);
}

template <int dim, int range, int rank>
void
NURBSElementAccessor<dim, range, rank>::
operator++()
{
    CartesianGridElementAccessor<dim> &grid_element_accessor = this->as_cartesian_grid_element_accessor();
    ++grid_element_accessor;
    ++bspline_element_accessor_;
}



IGA_NAMESPACE_CLOSE

#include <igatools/basis_functions/nurbs_element_accessor.inst>


