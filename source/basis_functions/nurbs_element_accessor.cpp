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


    this->reset_element_and_faces_cache(fill_flag, quad);
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

        const vector< Real > &weights = this->get_local_weights();

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


        const vector< Real > &weights = this->get_local_weights();

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


        const vector< Real > &weights = this->get_local_weights();


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
template <int deriv_order>
auto
NURBSElementAccessor< dim, range, rank >::
evaluate_basis_derivatives_at_points(const std::vector<Point<dim>> &points) const ->
ValueTable< Conditional< deriv_order==0,Value,Derivative<deriv_order> > >
{
    Assert(deriv_order >= 0 && deriv_order <= 2, ExcIndexRange(deriv_order,0,3));

    const int n_points = points.size();
    Assert(points.size() > 0, ExcEmptyObject());

    const int n_basis = this->get_num_basis();

    ValueTable< Conditional< deriv_order==0,Value,Derivative<deriv_order> > > result(n_basis,n_points);

    ComponentTable<int> n_basis_component;
    for (int comp = 0 ; comp < Space::n_components ; ++comp)
        n_basis_component(comp) = this->get_num_basis(comp);

    const auto W_vector = this->get_local_weights();

    if (deriv_order == 0)
    {
        const auto P_table = bspline_element_accessor_.evaluate_basis_values_at_points(points);
        const auto Q_table = bspline_element_accessor_.evaluate_field_values_at_points(this->get_local_weights(),points);

        LogStream out ;
        out << "P table=" << std::endl;
        P_table.print_info(out);
        out << std::endl;

        out << "Q table=" << std::endl;
        Q_table.print_info(out);
        out << std::endl;

        auto P_it = P_table.cbegin();

        auto R_it = result.begin();

        auto W_it = W_vector.cbegin();
        for (int comp = 0 ; comp < Space::n_components ; ++comp)
        {
            for (int ifn = 0 ; ifn < n_basis_component(comp) ; ifn++)
            {
                const auto &W = *W_it;
                auto Q_it = Q_table.cbegin();
                for (int jpt = 0 ; jpt < n_points ; jpt++)
                {
                    const auto &P = (*P_it)[comp];
                    const auto &Q = (*Q_it)[comp];
                    auto &R = (*R_it)[comp];

                    R = W * P / Q;

                    ++P_it;
                    ++Q_it;

                    ++R_it;
                } // end loop jpt
                ++W_it;
            } // end loop ifn
        } // end loop comp

    } // end if (deriv_order == 0)
    else if (deriv_order == 1)
    {
        const auto P_table = bspline_element_accessor_.evaluate_basis_values_at_points(points);
        const auto Q_table = bspline_element_accessor_.evaluate_field_values_at_points(this->get_local_weights(),points);

        const auto DP_table = bspline_element_accessor_.evaluate_basis_gradients_at_points(points);
        const auto DQ_table = bspline_element_accessor_.evaluate_field_gradients_at_points(this->get_local_weights(),points);

        auto P_it = P_table.cbegin();
        auto DP_it = DP_table.cbegin();

        auto DR_it = result.begin();


        auto W_it = W_vector.cbegin();
        for (int comp = 0 ; comp < Space::n_components ; ++comp)
        {
            for (int ifn = 0 ; ifn < n_basis_component(comp) ; ifn++)
            {
                const auto &W = *W_it;

                auto Q_it = Q_table.cbegin();
                auto DQ_it = DQ_table.cbegin();
                for (int jpt = 0 ; jpt < n_points ; jpt++)
                {
                    const auto &P = (*P_it)[comp];
                    const auto &Q = (*Q_it)[comp];

                    const auto &DP = (*DP_it);
                    const auto &DQ = (*DQ_it);

                    auto &DR = (*DR_it);

                    for (int i = 0 ; i < dim ; ++i)
                        DR(i)(comp) = (W / (Q*Q)) * (DP(i)(comp) * Q - P * DQ(i)(comp));

                    ++P_it;
                    ++Q_it;

                    ++DP_it;
                    ++DQ_it;

                    ++DR_it;
                } // end loop jpt
                ++W_it;
            } // end loop ifn
        } // end loop comp

    } // end else if (deriv_order == 1)
    else if (deriv_order == 2)
    {
        const auto P_table = bspline_element_accessor_.evaluate_basis_values_at_points(points);
        const auto Q_table = bspline_element_accessor_.evaluate_field_values_at_points(this->get_local_weights(),points);

        const auto DP_table = bspline_element_accessor_.evaluate_basis_gradients_at_points(points);
        const auto DQ_table = bspline_element_accessor_.evaluate_field_gradients_at_points(this->get_local_weights(),points);

        const auto D2P_table = bspline_element_accessor_.evaluate_basis_hessians_at_points(points);
        const auto D2Q_table = bspline_element_accessor_.evaluate_field_hessians_at_points(this->get_local_weights(),points);

        auto P_it = P_table.cbegin();
        auto DP_it = DP_table.cbegin();
        auto D2P_it = D2P_table.cbegin();

        auto D2R_it = result.begin();


        auto W_it = W_vector.cbegin();
        for (int comp = 0 ; comp < Space::n_components ; ++comp)
        {
            for (int ifn = 0 ; ifn < n_basis_component(comp) ; ifn++)
            {
                const auto &W = *W_it;

                auto Q_it = Q_table.cbegin();
                auto DQ_it = DQ_table.cbegin();
                auto D2Q_it = D2Q_table.cbegin();
                for (int jpt = 0 ; jpt < n_points ; jpt++)
                {
                    const auto &P = (*P_it)[comp];
                    const auto &Q = (*Q_it)[comp];

                    const auto &DP = (*DP_it);
                    const auto &DQ = (*DQ_it);

                    const auto &D2P = (*D2P_it);
                    const auto &D2Q = (*D2Q_it);

                    auto &D2R = (*D2R_it);

                    int der_entry_id = 0;
                    for (int i = 0 ; i < dim ; ++i)
                        for (int j = 0 ; j < dim ; ++j, ++der_entry_id)
                            D2R(der_entry_id)(comp) = (W/Q) *(D2P(der_entry_id)(comp)
                                                              - (P * D2Q(der_entry_id)(comp) +
                                                                 DP(i)(comp) * DQ(j)(comp) +
                                                                 DP(j)(comp) * DQ(i)(comp)) / Q +
                                                              DQ(i)(comp) * DQ(j)(comp) * (2.0 * P) / (Q*Q));
                    ++P_it;
                    ++Q_it;

                    ++DP_it;
                    ++DQ_it;

                    ++D2P_it;
                    ++D2Q_it;

                    ++D2R_it;
                } // end loop jpt
                ++W_it;
            } // end loop ifn
        } // end loop comp

    } // end else if (deriv_order == 2)

    return result;
}

template <int dim, int range, int rank >
void
NURBSElementAccessor< dim, range, rank >::
fill_values(const TopologyId<dim> &topology_id)
{
    Assert(topology_id.is_element() || topology_id.is_face(),
           ExcMessage("Only element or face topology is allowed."));

    auto &cache = this->get_values_cache(topology_id);
    Assert(cache.is_initialized(), ExcNotInitialized());


    // fills the cache of the BSplineElementAccessor
    bspline_element_accessor_.fill_values(topology_id);
    const auto &bspline_cache = bspline_element_accessor_.get_values_cache(topology_id);

    if (cache.flags_handler_.fill_values())
    {
        evaluate_nurbs_values(bspline_cache, cache.phi_);

        cache.flags_handler_.set_values_filled(true);
    }

    if (cache.flags_handler_.fill_gradients())
    {
        evaluate_nurbs_gradients(bspline_cache, cache.D1phi_);

        cache.flags_handler_.set_gradients_filled(true);
    }

    if (cache.flags_handler_.fill_hessians())
    {
        evaluate_nurbs_hessians(bspline_cache, cache.D2phi_);

        cache.flags_handler_.set_hessians_filled(true);
    }

    cache.set_filled(true);
}



template <int dim, int range, int rank >
vector<Real>
NURBSElementAccessor< dim, range, rank >::
get_local_weights() const
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
//*/

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


