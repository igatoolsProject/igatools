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

#ifndef ELLIPTIC_OPERATORS_SF_INTEGRATION_H_
#define ELLIPTIC_OPERATORS_SF_INTEGRATION_H_

#include <igatools/base/config.h>
#include <igatools/operators/elliptic_operators.h>
#include <igatools/utils/multi_array_utils.h>

#include <vector>

IGA_NAMESPACE_OPEN




/**
 * @brief Class containing the methods for the evaluation of some elliptic operators
 * (mass matrix, stiffness matrix, etc.) on a Bezier element
 * using the <em>sum-factorization quadrature approach</em>.
 *
 * All public methods take as input two PhysicalSpaceElementAccessor object
 * (one for the test space and the other for the trial space) and has an output argument that is
 * the local matrix relative to the elliptic operator that is evaluated.
 *
 *
 * @author M. Martinelli
 * @date 16 Apr 2014
 */

template <class PhysSpaceTest,class PhysSpaceTrial>
class EllipticOperatorsSFIntegration :
    public EllipticOperators<PhysSpaceTest,PhysSpaceTrial>
{
public:
    using base_t = EllipticOperators<PhysSpaceTest,PhysSpaceTrial>;

    using base_t::dim;
    using base_t::space_dim;

    using Clock = typename base_t::Clock;
    using TimePoint = typename base_t::TimePoint;
    using Duration = typename base_t::Duration;


    /** Type for the element accessor of the <em>test</em> physical space. */
    using ElemTest = typename base_t::ElemTest;

    /** Type for the element accessor of the <em>trial</em> physical space. */
    using ElemTrial = typename base_t::ElemTest;

    /** The constructors and destructor are inherithed from the base class. */
    using base_t::base_t;


    /** @name Assignment operators */
    ///@{
    /** Copy assignment operator. */
    EllipticOperatorsSFIntegration<PhysSpaceTest,PhysSpaceTrial> &
    operator=(const EllipticOperatorsSFIntegration<PhysSpaceTest,PhysSpaceTrial> &in) = default;


    /** Move assignment operator. */
    EllipticOperatorsSFIntegration<PhysSpaceTest,PhysSpaceTrial> &
    operator=(EllipticOperatorsSFIntegration<PhysSpaceTest,PhysSpaceTrial> &&in) = default;
    ///@}


    /**
     * This function evaluates the local (i.e. element-based) matrix \f$ A_e \f$
     * for witch its entries are
     * \f[
          (A_e)_{ij} = \int_{\Omega_e} \phi^{e,\text{test}}_i(x) C(x) \phi^{e,\text{trial}}_j(x) \; d \Omega.
       \f]
     * The matrix \f$ A_e \f$ is commonly referred as <em>local mass-matrix</em>.
     */
    virtual void eval_operator_u_v(
        const ElemTest &elem_test,
        const ElemTrial &elem_trial,
        const ValueVector<Real> &c,
        DenseMatrix &operator_u_v) const override final;

    /**
     * This function evaluates the local (i.e. element-based) matrix \f$ A_e \f$
     * for witch its entries are
     * \f[
          (A_e)_{ij} = \int_{\Omega_e} \sum_{r=1}^{sp\_dim} \sum_{s=1}^{sp\_dim}
          \bigl( \nabla \phi^{e,\text{test}}_i \bigr)_r
          \, C_{sr}(x) \,
          \bigl( \nabla \phi^{e,\text{trial}}_j \bigr)_s \; d \Omega.
       \f]
     * The matrix \f$ A_e \f$ is commonly referred as <em>local stiffness-matrix</em>.
     */
    virtual void eval_operator_gradu_gradv(
        const ElemTest &elem_test,
        const ElemTrial &elem_trial,
        const vector<TMatrix<space_dim,space_dim>> &coeffs,
        DenseMatrix &operator_gradu_gradv) const override final;


    /**
     * This function evaluates the local (i.e. element-based) vector \f$ f_e \f$
     * for witch its entries are
     * \f[
          (f_e)_{i} = \int_{\Omega_e} \phi^{e,\text{test}}_i
          f(x)  \; d \Omega.
       \f]
     */
    virtual void eval_operator_rhs_v(
        const ElemTest &elem_test,
        const ValueVector<typename PhysSpaceTrial::Value> &f,
        DenseVector &operator_rhs_v) const override final;


    /**
     * This function evaluates the local (i.e. element-based) matrix \f$ A_e \f$
     * for witch its entries are
     * \f[
          (A_e)_{ij} = \int_{\Omega_e} \sum_{s=1}^{sp\_dim}
          \phi^{e,\text{test}}_i
          \, \beta_{s}(x) \,
          \bigl( \nabla \phi^{e,\text{trial}}_j \bigr)_s \; d \Omega
          = \int_{\Omega_e}
          \phi^{e,\text{test}}_i
          \, \vec{\beta}(x) \, \cdot \,
          \nabla \phi^{e,\text{trial}}_j \; d \Omega .
       \f]
     */
    virtual void eval_operator_gradu_v(
        const ElemTest &elem_test,
        const ElemTrial &elem_trial,
        const ValueVector<typename PhysSpaceTrial::Gradient> &beta,
        DenseMatrix &operator_gradu_v) const override final;

private:

    /**
     * Returns the quadrature weights multiplied by the one-dimensional basis
     * for test and trial space, as needed by the integration using the
     * sum factorization technique.
     * \f[ J[i]_{\theta_i,\alpha_i,\beta_i} =
       w_{i,\theta_i}
       \phi^{\text{test}}_{\beta_i}(x_{\theta_i})
       \phi^{\text{trial}}_{\alpha_i}(x_{\theta_i})
       \f]
     * where \f$ w_{i,\theta_i} \f$ is the quadrature weight
     * relative to the the point \f$x_{\theta_i}\f$
     * along the \f$i\f$-th direction.
     */
    std::array<DynamicMultiArray<Real,3>,dim>
    evaluate_w_phi1Dtrial_phi1Dtest(
        const std::array<ValueTable<Real>,dim> &phi_1D_test,
        const std::array<ValueTable<Real>,dim> &phi_1D_trial,
        const TensorProductArray<dim> &quad_weights,
        const std::array<Real,dim> &length_element_edge) const;
};





template <int dim, int r=dim>
class SumFactorizationIntegrator
{
public:
    void
    operator()(
        const bool is_symmetric,
        const TensorSize<dim> &t_size_theta,
        const TensorSize<dim> &t_size_alpha,
        const TensorSize<dim> &t_size_beta,
        const std::array<DynamicMultiArray<Real,3>,dim> &J,
        const DynamicMultiArray<Real,3> &Cpre,
        DenseMatrix &local_operator) const
    {
        const int k = dim-r+1;

        // (alpha_1,...alpha_{k-1})
        TensorSize<k-1> t_size_alpha_1_km1;
        // (beta_1,...beta_{k-1})
        TensorSize<k-1> t_size_beta_1_km1;
        for (int i = 0 ; i < k-1 ; ++i)
        {
            t_size_alpha_1_km1[i] = t_size_alpha[i];
            t_size_beta_1_km1 [i] = t_size_beta[i];
        }

        // (alpha_1,...alpha_k)
        TensorSize<k> t_size_alpha_1_k;
        // (beta_1,...beta_k)
        TensorSize<k> t_size_beta_1_k;
        for (int i = 0 ; i < k ; ++i)
        {
            t_size_alpha_1_k[i] = t_size_alpha[i];
            t_size_beta_1_k [i] = t_size_beta[i];
        }

        // (theta_{k+1},...theta_{dim})
        TensorSize<dim-k> t_size_theta_kp1_d;
        for (int i = 0 ; i < dim-k ; ++i)
            t_size_theta_kp1_d[i] = t_size_theta[i+k];


        // (theta_k,...theta_{dim})
        TensorSize<dim-k+1> t_size_theta_k_d;
        for (int i = 0 ; i <= dim-k ; ++i)
            t_size_theta_k_d[i] = t_size_theta[i+k-1];


        const Size f_size_theta_kp1_d = (dim-k>0)?t_size_theta_kp1_d.flat_size():1;
        const Size f_size_alpha_1_km1 = t_size_alpha_1_km1.flat_size();
        const Size f_size_beta_1_km1  = t_size_beta_1_km1.flat_size();


        const auto &Jk = J[k-1];
        TensorSize<3> t_size_Jk = Jk.tensor_size();
        Assert(t_size_Jk[0] == t_size_theta[k-1],ExcDimensionMismatch(t_size_Jk[0],t_size_theta[k-1]));
        Assert(t_size_Jk[1] == t_size_alpha[k-1],ExcDimensionMismatch(t_size_Jk[1],t_size_alpha[k-1]));
        Assert(t_size_Jk[2] == t_size_beta [k-1],ExcDimensionMismatch(t_size_Jk[2],t_size_beta [k-1]));
        TensorIndex<3> t_wgt_Jk = MultiArrayUtils<3>::compute_weight(t_size_Jk);


        const TensorSize<3> t_size_Cpre = Cpre.tensor_size();
        Assert(t_size_Cpre[0] == t_size_theta_k_d.flat_size(),
               ExcDimensionMismatch(t_size_Cpre[0],t_size_theta_k_d.flat_size()));
        Assert(t_size_Cpre[1] == t_size_alpha_1_km1.flat_size(),
               ExcDimensionMismatch(t_size_Cpre[1],t_size_alpha_1_km1.flat_size()));
        Assert(t_size_Cpre[2] == t_size_beta_1_km1.flat_size(),
               ExcDimensionMismatch(t_size_Cpre[2],t_size_beta_1_km1.flat_size()));
        TensorIndex<3> t_wgt_Cpre = MultiArrayUtils<3>::compute_weight(t_size_Cpre);


        TensorSize<3> t_size_Cpost;
        t_size_Cpost[0] = f_size_theta_kp1_d;
        t_size_Cpost[1] = t_size_alpha_1_k.flat_size();
        t_size_Cpost[2] = t_size_beta_1_k.flat_size();
        TensorIndex<3> t_wgt_Cpost = MultiArrayUtils<3>::compute_weight(t_size_Cpost);
        DynamicMultiArray<Real,3> Cpost(t_size_Cpost);


        TensorIndex<3> tid_Jk;
        TensorIndex<3> tid_Cpre;
        TensorIndex<3> tid_Cpost;

        if (!is_symmetric)
        {

            tid_Cpost[2] = 0;
            tid_Jk[0] = 0;
            for (Index flat_beta_k_1 = 0 ; flat_beta_k_1 < f_size_beta_1_km1 ; ++flat_beta_k_1)
            {
                tid_Cpre[2] = flat_beta_k_1;

                for (int beta_k = 0 ; beta_k < t_size_beta[k-1] ; ++beta_k)
                {
                    tid_Jk[2] = beta_k;
                    tid_Cpost[1] = 0 ;

                    for (Index flat_alpha_k_1 = 0 ; flat_alpha_k_1 < f_size_alpha_1_km1 ; ++flat_alpha_k_1)
                    {
                        tid_Cpre[1] = flat_alpha_k_1;

                        for (int alpha_k = 0 ; alpha_k < t_size_alpha[k-1] ; ++alpha_k)
                        {
                            tid_Jk[1] = alpha_k;
                            tid_Cpre[0] = 0 ;
                            tid_Cpost[0] = 0 ;
                            for (Index fid_theta_kp1_d = 0 ; fid_theta_kp1_d < f_size_theta_kp1_d ; ++fid_theta_kp1_d)
                            {
                                const Index f_id_Cpost = MultiArrayUtils<3>::tensor_to_flat_index(tid_Cpost,t_wgt_Cpost);

                                const Index f_id_Cpre = MultiArrayUtils<3>::tensor_to_flat_index(tid_Cpre,t_wgt_Cpre);

                                const Index f_id_Jk = MultiArrayUtils<3>::tensor_to_flat_index(tid_Jk,t_wgt_Jk);

                                Cpost[f_id_Cpost] = std::inner_product(
                                                        &Jk[f_id_Jk],
                                                        &Jk[f_id_Jk]+t_size_theta[k-1],
                                                        &Cpre[f_id_Cpre],
                                                        0.0);

                                tid_Cpost[0]++;
                                tid_Cpre [0] += t_size_theta[k-1];
                            } //end loop flat_theta_k_1

                            tid_Cpost[1]++;

                        } // end loop alpha_k
                    } // end loop flat_alpha_k_1

                    tid_Cpost[2]++;
                } // end loop beta_k

            } // end loop flat_beta_k_1


        } // end if(!is_symmetric)
        else
        {
            using MAUtils_k = MultiArrayUtils<k>;
            using MAUtils_km1 = MultiArrayUtils<k-1>;

            const Size f_size_alpha_1_k = t_size_alpha_1_k.flat_size();
            const TensorIndex<k> wgt_alpha_1_k = MAUtils_k::compute_weight(t_size_alpha_1_k);
            const TensorIndex<k-1> wgt_alpha_1_km1 = MAUtils_km1::compute_weight(t_size_alpha_1_km1);

            const Size f_size_beta_1_k = t_size_beta_1_k.flat_size();
            const TensorIndex<k> wgt_beta_1_k = MAUtils_k::compute_weight(t_size_beta_1_k);
            const TensorIndex<k-1> wgt_beta_1_km1 = MAUtils_km1::compute_weight(t_size_beta_1_km1);

            TensorIndex<k-1> tid_alpha_1_km1;
            TensorIndex<k-1> tid_beta_1_km1;
            for (Index fid_beta_1_k = 0 ; fid_beta_1_k < f_size_beta_1_k ; ++fid_beta_1_k)
            {
                const TensorIndex<k> tid_beta_1_k =
                    MAUtils_k::flat_to_tensor_index(fid_beta_1_k,wgt_beta_1_k);

                for (int i = 0 ; i < k-1 ; ++i)
                    tid_beta_1_km1(i) = tid_beta_1_k(i);

                const Index beta_k = tid_beta_1_k(k-1);

                const Index fid_beta_1_km1 =
                    (k>1)?MAUtils_km1::tensor_to_flat_index(tid_beta_1_km1,wgt_beta_1_km1):0;


                for (Index fid_alpha_1_k = fid_beta_1_k ; fid_alpha_1_k < f_size_alpha_1_k ; ++fid_alpha_1_k)
                {
                    const TensorIndex<k> tid_alpha_1_k =
                        MAUtils_k::flat_to_tensor_index(fid_alpha_1_k,wgt_alpha_1_k);
                    for (int i = 0 ; i < k-1 ; ++i)
                        tid_alpha_1_km1(i) = tid_alpha_1_k(i);

                    const Index alpha_k = tid_alpha_1_k(k-1);

                    const Index fid_alpha_1_km1 =
                        (k>1)?MAUtils_km1::tensor_to_flat_index(tid_alpha_1_km1,wgt_alpha_1_km1):0;

                    tid_Cpre[1] = std::max(fid_alpha_1_km1,fid_beta_1_km1);
                    tid_Cpre[2] = std::min(fid_alpha_1_km1,fid_beta_1_km1);

                    tid_Jk[1] = std::max(alpha_k,beta_k);
                    tid_Jk[2] = std::min(alpha_k,beta_k);

                    tid_Cpost[2] = fid_beta_1_k ;
                    tid_Cpost[1] = fid_alpha_1_k;

                    tid_Cpre [0] = 0;
                    tid_Cpost[0] = 0;
                    for (Index fid_theta_kp1_d = 0 ; fid_theta_kp1_d < f_size_theta_kp1_d ; ++fid_theta_kp1_d)
                    {
                        const Index f_id_Cpost = MultiArrayUtils<3>::tensor_to_flat_index(tid_Cpost,t_wgt_Cpost);

                        const Index f_id_Cpre = MultiArrayUtils<3>::tensor_to_flat_index(tid_Cpre,t_wgt_Cpre);

                        const Index f_id_Jk = MultiArrayUtils<3>::tensor_to_flat_index(tid_Jk,t_wgt_Jk);

                        Cpost[f_id_Cpost] = std::inner_product(
                                                &Jk[f_id_Jk],
                                                &Jk[f_id_Jk]+t_size_theta[k-1],
                                                &Cpre[f_id_Cpre],
                                                0.0);


                        tid_Cpost[0]++;
                        tid_Cpre [0] += t_size_theta[k-1];
                    } //end loop flat_theta_k_1

                }// end loop fid_alpha_1_k

            } // end loop fid_beta_1_k
            //*/

        } // end if (symmetric)

        SumFactorizationIntegrator<dim,r-1> integrate_sf;
        integrate_sf(
            is_symmetric,
            t_size_theta,
            t_size_alpha,
            t_size_beta,
            J,
            Cpost,
            local_operator);
    }
};


template <int dim>
class SumFactorizationIntegrator<dim,1>
{
public:
    void operator()(
        const bool is_symmetric,
        const TensorSize<dim> &t_size_theta,
        const TensorSize<dim> &t_size_alpha,
        const TensorSize<dim> &t_size_beta,
        const std::array<DynamicMultiArray<Real,3>,dim> &J,
        const DynamicMultiArray<Real,3> &Cpre,
        DenseMatrix &local_operator) const
    {
        const int k = dim;

        // (alpha_1,...alpha_{k-1})
        TensorSize<k-1> t_size_alpha_1_km1;
        // (beta_1,...beta_{k-1})
        TensorSize<k-1> t_size_beta_1_km1;
        for (int i = 0 ; i < k-1 ; ++i)
        {
            t_size_alpha_1_km1[i] = t_size_alpha[i];
            t_size_beta_1_km1 [i] = t_size_beta[i];
        }

        // (alpha_1,...alpha_k)
        TensorSize<k> t_size_alpha_1_k;
        // (beta_1,...beta_k)
        TensorSize<k> t_size_beta_1_k;
        for (int i = 0 ; i < k ; ++i)
        {
            t_size_alpha_1_k[i] = t_size_alpha[i];
            t_size_beta_1_k [i] = t_size_beta[i];
        }

        // (theta_{k+1},...theta_{dim})
        TensorSize<dim-k> t_size_theta_kp1_d;
        for (int i = 0 ; i < dim-k ; ++i)
            t_size_theta_kp1_d[i] = t_size_theta[i+k];


        // (theta_k,...theta_{dim})
        TensorSize<dim-k+1> t_size_theta_k_d;
        for (int i = 0 ; i <= dim-k ; ++i)
            t_size_theta_k_d[i] = t_size_theta[i+k-1];


        const Size f_size_alpha_1_km1 = t_size_alpha_1_km1.flat_size();
        const Size f_size_beta_1_km1  = t_size_beta_1_km1.flat_size();


        const auto &Jk = J[k-1];
        TensorSize<3> t_size_Jk = Jk.tensor_size();
        Assert(t_size_Jk[0] == t_size_theta[k-1],ExcDimensionMismatch(t_size_Jk[0],t_size_theta[k-1]));
        Assert(t_size_Jk[1] == t_size_alpha[k-1],ExcDimensionMismatch(t_size_Jk[1],t_size_alpha[k-1]));
        Assert(t_size_Jk[2] == t_size_beta [k-1],ExcDimensionMismatch(t_size_Jk[2],t_size_beta [k-1]));
        TensorIndex<3> t_wgt_Jk = MultiArrayUtils<3>::compute_weight(t_size_Jk);


        const TensorSize<3> t_size_Cpre = Cpre.tensor_size();
        Assert(t_size_Cpre[0] == t_size_theta_k_d.flat_size(),
               ExcDimensionMismatch(t_size_Cpre[0],t_size_theta_k_d.flat_size()));
        Assert(t_size_Cpre[1] == t_size_alpha_1_km1.flat_size(),
               ExcDimensionMismatch(t_size_Cpre[1],t_size_alpha_1_km1.flat_size()));
        Assert(t_size_Cpre[2] == t_size_beta_1_km1.flat_size(),
               ExcDimensionMismatch(t_size_Cpre[2],t_size_beta_1_km1.flat_size()));
        TensorIndex<3> t_wgt_Cpre = MultiArrayUtils<3>::compute_weight(t_size_Cpre);


        TensorIndex<3> tid_Jk;
        TensorIndex<3> tid_Cpre;
        TensorIndex<3> tid_Cpost;


        const Size f_size_alpha_1_k = t_size_alpha_1_k.flat_size();
        const Size f_size_beta_1_k = t_size_beta_1_k.flat_size();

        Assert(local_operator.get_num_rows() == f_size_beta_1_k,
               ExcDimensionMismatch(local_operator.get_num_rows(),f_size_beta_1_k));
        Assert(local_operator.get_num_cols() == f_size_alpha_1_k,
               ExcDimensionMismatch(local_operator.get_num_cols(),f_size_alpha_1_k));

        if (!is_symmetric)
        {

            tid_Jk[0] = 0;
            tid_Cpre[0] = 0 ;
            for (Index flat_beta_k_1 = 0 ; flat_beta_k_1 < f_size_beta_1_km1 ; ++flat_beta_k_1)
            {
                tid_Cpre[2] = flat_beta_k_1;

                for (int beta_k = 0 ; beta_k < t_size_beta[k-1] ; ++beta_k)
                {
                    tid_Jk[2] = beta_k;

                    const Index f_id_test = flat_beta_k_1 + beta_k;

                    for (Index flat_alpha_k_1 = 0 ; flat_alpha_k_1 < f_size_alpha_1_km1 ; ++flat_alpha_k_1)
                    {
                        tid_Cpre[1] = flat_alpha_k_1;

                        for (int alpha_k = 0 ; alpha_k < t_size_alpha[k-1] ; ++alpha_k)
                        {
                            tid_Jk[1] = alpha_k;

                            const Index f_id_trial = flat_alpha_k_1 + alpha_k;

                            const Index f_id_Cpre = MultiArrayUtils<3>::tensor_to_flat_index(tid_Cpre,t_wgt_Cpre);

                            const Index f_id_Jk = MultiArrayUtils<3>::tensor_to_flat_index(tid_Jk,t_wgt_Jk);

                            local_operator(f_id_test,f_id_trial) =
                                std::inner_product(
                                    &Jk[f_id_Jk],
                                    &Jk[f_id_Jk]+t_size_theta[k-1],
                                    &Cpre[f_id_Cpre],
                                    0.0);
                        } // end loop alpha_k
                    } // end loop flat_alpha_k_1

                } // end loop beta_k

            } // end loop flat_beta_k_1


        } // end if(!is_symmetric)
        else
        {
            using MAUtils_k = MultiArrayUtils<k>;
            using MAUtils_km1 = MultiArrayUtils<k-1>;

            const TensorIndex<k> wgt_alpha_1_k = MAUtils_k::compute_weight(t_size_alpha_1_k);
            const TensorIndex<k-1> wgt_alpha_1_km1 = MAUtils_km1::compute_weight(t_size_alpha_1_km1);

            const TensorIndex<k> wgt_beta_1_k = MAUtils_k::compute_weight(t_size_beta_1_k);
            const TensorIndex<k-1> wgt_beta_1_km1 = MAUtils_km1::compute_weight(t_size_beta_1_km1);

            TensorIndex<k-1> tid_alpha_1_km1;
            TensorIndex<k-1> tid_beta_1_km1;
            tid_Cpre [0] = 0;
            for (Index fid_beta_1_k = 0 ; fid_beta_1_k < f_size_beta_1_k ; ++fid_beta_1_k)
            {
                const TensorIndex<k> tid_beta_1_k =
                    MAUtils_k::flat_to_tensor_index(fid_beta_1_k,wgt_beta_1_k);

                for (int i = 0 ; i < k-1 ; ++i)
                    tid_beta_1_km1(i) = tid_beta_1_k(i);

                const Index beta_k = tid_beta_1_k(k-1);

                const Index fid_beta_1_km1 =
                    (k>1)?MAUtils_km1::tensor_to_flat_index(tid_beta_1_km1,wgt_beta_1_km1):0;


                for (Index fid_alpha_1_k = fid_beta_1_k ; fid_alpha_1_k < f_size_alpha_1_k ; ++fid_alpha_1_k)
                {
                    const TensorIndex<k> tid_alpha_1_k =
                        MAUtils_k::flat_to_tensor_index(fid_alpha_1_k,wgt_alpha_1_k);
                    for (int i = 0 ; i < k-1 ; ++i)
                        tid_alpha_1_km1(i) = tid_alpha_1_k(i);

                    const Index alpha_k = tid_alpha_1_k(k-1);

                    const Index fid_alpha_1_km1 =
                        (k>1)?MAUtils_km1::tensor_to_flat_index(tid_alpha_1_km1,wgt_alpha_1_km1):0;

                    tid_Cpre[1] = std::max(fid_alpha_1_km1,fid_beta_1_km1);
                    tid_Cpre[2] = std::min(fid_alpha_1_km1,fid_beta_1_km1);

                    tid_Jk[1] = std::max(alpha_k,beta_k);
                    tid_Jk[2] = std::min(alpha_k,beta_k);

                    const Index f_id_Cpre = MultiArrayUtils<3>::tensor_to_flat_index(tid_Cpre,t_wgt_Cpre);

                    const Index f_id_Jk = MultiArrayUtils<3>::tensor_to_flat_index(tid_Jk,t_wgt_Jk);

                    local_operator(fid_beta_1_k,fid_alpha_1_k) =
                        std::inner_product(
                            &Jk[f_id_Jk],
                            &Jk[f_id_Jk]+t_size_theta[k-1],
                            &Cpre[f_id_Cpre],
                            0.0);

                }// end loop fid_alpha_1_k

            } // end loop fid_beta_1_k
            //*/

            // here we copy the upper triangular part of the matrix on the lower triangular part
            for (int test_id = 0 ; test_id < f_size_beta_1_k ; ++test_id)
                for (int trial_id = 0; trial_id < test_id ; ++trial_id)
                    local_operator(test_id,trial_id) = local_operator(trial_id,test_id);

        } // end if (symmetric)
    }
};

#define SPECIALIZED
#ifdef SPECIALIZED
template <>
class SumFactorizationIntegrator<1,1>
{
public:
    void operator()(
        const bool is_symmetric,
        const TensorSize<1> &t_size_theta,
        const TensorSize<1> &t_size_alpha,
        const TensorSize<1> &t_size_beta,
        const std::array<DynamicMultiArray<Real,3>,1> &J,
        const DynamicMultiArray<Real,3> &C,
        DenseMatrix &local_operator) const
    {
        Assert(t_size_alpha.flat_size() == local_operator.get_num_cols(),
               ExcDimensionMismatch(t_size_alpha.flat_size(),local_operator.get_num_cols()));
        Assert(t_size_beta.flat_size() == local_operator.get_num_rows(),
               ExcDimensionMismatch(t_size_beta.flat_size(),local_operator.get_num_rows()));


        TensorIndex<3> t_id_J;
        TensorIndex<3> t_id_C;

        if (!is_symmetric)
        {
            //--------------------------------------------------------------
            t_id_C[1] = 0;
            t_id_C[2] = 0;

            for (Index beta_0 = 0 ; beta_0 < t_size_beta[0] ; ++beta_0)
            {
                t_id_J[2] = beta_0;
                for (Index alpha_0 = 0 ; alpha_0 < t_size_alpha[0] ; ++alpha_0)
                {
                    t_id_J[1] = alpha_0;

                    Real sum = 0.0;
                    for (Index theta_0 = 0; theta_0 < t_size_theta[0] ; ++theta_0)
                    {
                        t_id_J[0] = theta_0;
                        t_id_C[0] = theta_0;
                        sum += C(t_id_C) * J[0](t_id_J);
                    }

                    local_operator(beta_0,alpha_0) = sum;
                }
            }
            //--------------------------------------------------------------
        } // end if (!is_symmetric)
        else
        {
            //--------------------------------------------------------------
            t_id_C[1] = 0;
            t_id_C[2] = 0;
            for (Index beta_0 = 0 ; beta_0 < t_size_beta[0] ; ++beta_0)
            {
                t_id_J[2] = beta_0;
                for (Index alpha_0 = beta_0 ; alpha_0 < t_size_alpha[0] ; ++alpha_0)
                {
                    t_id_J[1] = alpha_0;

                    Real sum = 0.0;
                    for (Index theta_0 = 0; theta_0 < t_size_theta[0] ; ++theta_0)
                    {
                        t_id_J[0] = theta_0;
                        t_id_C[0] = theta_0;
                        sum += C(t_id_C) * J[0](t_id_J);
                    }

                    local_operator(beta_0,alpha_0) = sum;
                }
            }

            // here we copy the upper triangular part of the matrix on the lower triangular part
            const Size n_basis_test = local_operator.get_num_rows();
            for (int test_id = 0 ; test_id < n_basis_test ; ++test_id)
                for (int trial_id = 0; trial_id < test_id ; ++trial_id)
                    local_operator(test_id,trial_id) = local_operator(trial_id,test_id);

            //--------------------------------------------------------------

        } // end if (is_symmetric)
    }
};

template <>
class SumFactorizationIntegrator<2,2>
{
public:
    void operator()(
        const bool is_symmetric,
        const TensorSize<2> &t_size_theta,
        const TensorSize<2> &t_size_alpha,
        const TensorSize<2> &t_size_beta,
        const std::array<DynamicMultiArray<Real,3>,2> &J,
        const DynamicMultiArray<Real,3> &C,
        DenseMatrix &local_operator) const
    {
        Assert(t_size_alpha.flat_size() == local_operator.get_num_cols(),
               ExcDimensionMismatch(t_size_alpha.flat_size(),local_operator.get_num_cols()));
        Assert(t_size_beta.flat_size() == local_operator.get_num_rows(),
               ExcDimensionMismatch(t_size_beta.flat_size(),local_operator.get_num_rows()));


        TensorIndex<3> t_id_J;
        TensorIndex<3> t_id_C;


        //--------------------------------------------------------------
        TensorIndex<3> t_id_C1;
        TensorSize<3> t_size_C1;
        t_size_C1[0] = t_size_theta[1];
        t_size_C1[1] = t_size_alpha[0];
        t_size_C1[2] = t_size_beta[0];
        DynamicMultiArray<Real,3> C1(t_size_C1);
        //--------------------------------------------------------------


        //--------------------------------------------------------------
        t_id_C[1] = 0;
        t_id_C[2] = 0;
        for (Index beta_0 = 0 ; beta_0 < t_size_beta[0] ; ++beta_0)
        {
            t_id_J [2] = beta_0;
            t_id_C1[2] = beta_0;
            for (Index alpha_0 = 0 ; alpha_0 < t_size_alpha[0] ; ++alpha_0)
            {
                t_id_J [1] = alpha_0;
                t_id_C1[1] = alpha_0;

                Index theta_0_1 = 0;
                for (Index theta_1 = 0; theta_1 < t_size_theta[1] ; ++theta_1)
                {
                    Real sum = 0.0;
                    for (Index theta_0 = 0; theta_0 < t_size_theta[0] ; ++theta_0,++theta_0_1)
                    {
                        t_id_J[0] = theta_0;
                        t_id_C[0] = theta_0_1;
                        sum += C(t_id_C) * J[0](t_id_J);
                    }

                    t_id_C1[0] = theta_1;
                    C1(t_id_C1) = sum;
                } //end loop theta_1
            } //end loop alpha_0
        } // end loop beta_0
        //--------------------------------------------------------------


        if (!is_symmetric)
        {
            //--------------------------------------------------------------
            for (Index beta_1 = 0 ; beta_1 < t_size_beta[1] ; ++beta_1)
            {
                t_id_J[2] = beta_1;
                for (Index beta_0 = 0 ; beta_0 < t_size_beta[0] ; ++beta_0)
                {
                    Index beta_0_1 = beta_1*t_size_beta[0] + beta_0;

                    t_id_C1[2] = beta_0;

                    for (Index alpha_1 = 0 ; alpha_1 < t_size_alpha[1] ; ++alpha_1)
                    {
                        t_id_J[1] = alpha_1;
                        for (Index alpha_0 = 0 ; alpha_0 < t_size_alpha[0] ; ++alpha_0)
                        {
                            Index alpha_0_1 = alpha_1*t_size_alpha[0] + alpha_0;

                            t_id_C1[1] = alpha_0;

                            Real sum = 0.0;
                            for (Index theta_1 = 0; theta_1 < t_size_theta[1] ; ++theta_1)
                            {
                                t_id_J [0] = theta_1;
                                t_id_C1[0] = theta_1;
                                sum += C1(t_id_C1) * J[1](t_id_J);
                            } // end loop theta_1

                            local_operator(beta_0_1,alpha_0_1) = sum;
                        } //end loop alpha_0
                    } //end loop alpha_1
                } // end loop beta_0
            } // end loop beta_1
            //--------------------------------------------------------------
        }//end if (!is_symmetric)
        else
        {
            //--------------------------------------------------------------
            for (Index beta_1 = 0 ; beta_1 < t_size_beta[1] ; ++beta_1)
            {
                t_id_J[2] = beta_1;
                for (Index beta_0 = 0 ; beta_0 < t_size_beta[0] ; ++beta_0)
                {
                    Index beta_0_1 = beta_1*t_size_beta[0] + beta_0;

                    t_id_C1[2] = beta_0;

                    for (Index alpha_1 = beta_1 ; alpha_1 < t_size_alpha[1] ; ++alpha_1)
                    {
                        t_id_J[1] = alpha_1;
                        for (Index alpha_0 = 0 ; alpha_0 < t_size_alpha[0] ; ++alpha_0)
                        {
                            Index alpha_0_1 = alpha_1*t_size_alpha[0] + alpha_0;

                            t_id_C1[1] = alpha_0;

                            Real sum = 0.0;
                            for (Index theta_1 = 0; theta_1 < t_size_theta[1] ; ++theta_1)
                            {
                                t_id_J [0] = theta_1;
                                t_id_C1[0] = theta_1;
                                sum += C1(t_id_C1) * J[1](t_id_J);
                            } // end loop theta_1

                            local_operator(beta_0_1,alpha_0_1) = sum;
                        } //end loop alpha_0
                    } //end loop alpha_1
                } // end loop beta_0
            } // end loop beta_1


            // here we copy the upper triangular part of the matrix on the lower triangular part
            const Size n_basis_test = local_operator.get_num_rows();
            for (int test_id = 0 ; test_id < n_basis_test ; ++test_id)
                for (int trial_id = 0; trial_id < test_id ; ++trial_id)
                    local_operator(test_id,trial_id) = local_operator(trial_id,test_id);

            //--------------------------------------------------------------
        }//end if (is_symmetric)
    }
};


template <>
class SumFactorizationIntegrator<3,3>
{
public:
    void operator()(
        const bool is_symmetric,
        const TensorSize<3> &t_size_theta,
        const TensorSize<3> &t_size_alpha,
        const TensorSize<3> &t_size_beta,
        const std::array<DynamicMultiArray<Real,3>,3> &J,
        const DynamicMultiArray<Real,3> &C,
        DenseMatrix &local_operator) const
    {
        Assert(t_size_alpha.flat_size() == local_operator.get_num_cols(),
               ExcDimensionMismatch(t_size_alpha.flat_size(),local_operator.get_num_cols()));
        Assert(t_size_beta.flat_size() == local_operator.get_num_rows(),
               ExcDimensionMismatch(t_size_beta.flat_size(),local_operator.get_num_rows()));



        TensorIndex<3> t_id_J;
        TensorIndex<3> t_id_C;

        //--------------------------------------------------------------
        TensorIndex<3> t_id_C1;
        TensorSize<3> t_size_C1;
        t_size_C1[0] = t_size_theta[1] * t_size_theta[2];
        t_size_C1[1] = t_size_alpha[0];
        t_size_C1[2] = t_size_beta[0];
        DynamicMultiArray<Real,3> C1(t_size_C1);
        //--------------------------------------------------------------


        //--------------------------------------------------------------
        TensorIndex<3> t_id_C2;
        TensorSize<3> t_size_C2;
        t_size_C2[0] = t_size_theta[2];
        t_size_C2[1] = t_size_alpha[0] * t_size_alpha[1];
        t_size_C2[2] = t_size_beta [0] * t_size_beta [1];
        DynamicMultiArray<Real,3> C2(t_size_C2);
        //--------------------------------------------------------------


        //--------------------------------------------------------------
        t_id_C[1] = 0;
        t_id_C[2] = 0;
        for (Index beta_0 = 0 ; beta_0 < t_size_beta[0] ; ++beta_0)
        {
            t_id_J [2] = beta_0;
            t_id_C1[2] = beta_0;
            for (Index alpha_0 = 0 ; alpha_0 < t_size_alpha[0] ; ++alpha_0)
            {
                t_id_J [1] = alpha_0;
                t_id_C1[1] = alpha_0;

                Index theta_1_2 = 0;
                Index theta_0_1_2 = 0;
                for (Index theta_2 = 0; theta_2 < t_size_theta[2] ; ++theta_2)
                {
                    for (Index theta_1 = 0; theta_1 < t_size_theta[1] ; ++theta_1,++theta_1_2)
                    {
                        Real sum = 0.0;
                        for (Index theta_0 = 0; theta_0 < t_size_theta[0] ; ++theta_0,++theta_0_1_2)
                        {
                            t_id_J[0] = theta_0;
                            t_id_C[0] = theta_0_1_2;
                            sum += C(t_id_C) * J[0](t_id_J);
                        }

                        t_id_C1[0] = theta_1_2;
                        C1(t_id_C1) = sum;
                    } // end loop theta_1
                } // end loop theta_2
            } // end loop alpha_0
        } // end loop beta_0
        //--------------------------------------------------------------



        //--------------------------------------------------------------
        for (Index beta_1 = 0 ; beta_1 < t_size_beta[1] ; ++beta_1)
        {
            t_id_J[2] = beta_1;
            for (Index beta_0 = 0 ; beta_0 < t_size_beta[0] ; ++beta_0)
            {
                Index beta_0_1 = beta_1*t_size_beta[0] + beta_0;

                t_id_C1[2] = beta_0;
                t_id_C2[2] = beta_0_1;

                for (Index alpha_1 = 0 ; alpha_1 < t_size_alpha[1] ; ++alpha_1)
                {
                    t_id_J[1] = alpha_1;
                    for (Index alpha_0 = 0 ; alpha_0 < t_size_alpha[0] ; ++alpha_0)
                    {
                        Index alpha_0_1 = alpha_1*t_size_alpha[0] + alpha_0;

                        t_id_C1[1] = alpha_0;
                        t_id_C2[1] = alpha_0_1;

                        Index theta_1_2 = 0;
                        for (Index theta_2 = 0; theta_2 < t_size_theta[2] ; ++theta_2)
                        {
                            Real sum = 0.0;
                            for (Index theta_1 = 0; theta_1 < t_size_theta[1] ; ++theta_1,++theta_1_2)
                            {
                                t_id_J [0] = theta_1;
                                t_id_C1[0] = theta_1_2;
                                sum += C1(t_id_C1) * J[1](t_id_J);
                            } // end loop theta_1

                            t_id_C2[0] = theta_2;
                            C2(t_id_C2) = sum;
                        } // end loop theta_2
                    } //end loop alpha_0
                } //end loop alpha_1
            } // end loop beta_0
        } // end loop beta_1
        //--------------------------------------------------------------

        if (!is_symmetric)
        {

            //--------------------------------------------------------------
            for (Index beta_2 = 0 ; beta_2 < t_size_beta[2] ; ++beta_2)
            {
                t_id_J[2] = beta_2;
                for (Index beta_1 = 0 ; beta_1 < t_size_beta[1] ; ++beta_1)
                {
                    for (Index beta_0 = 0 ; beta_0 < t_size_beta[0] ; ++beta_0)
                    {
                        Index beta_0_1 = beta_1 * t_size_beta[0] + beta_0 ;
                        Index beta_0_1_2 = (beta_2*t_size_beta[1] + beta_1) * t_size_beta[0] + beta_0 ;

                        t_id_C2[2] = beta_0_1 ;

                        for (Index alpha_2 = 0 ; alpha_2 < t_size_alpha[2] ; ++alpha_2)
                        {
                            t_id_J[1] = alpha_2;
                            for (Index alpha_1 = 0 ; alpha_1 < t_size_alpha[1] ; ++alpha_1)
                            {
                                for (Index alpha_0 = 0 ; alpha_0 < t_size_alpha[0] ; ++alpha_0)
                                {
                                    Index alpha_0_1 = alpha_1 * t_size_alpha[0] + alpha_0 ;
                                    Index alpha_0_1_2 = (alpha_2*t_size_alpha[1] + alpha_1) * t_size_alpha[0] + alpha_0 ;

                                    t_id_C2[1] = alpha_0_1 ;

                                    Real sum = 0.0;
                                    for (Index theta_2 = 0; theta_2 < t_size_theta[2] ; ++theta_2)
                                    {
                                        t_id_J [0] = theta_2;
                                        t_id_C2[0] = theta_2;
                                        sum += C2(t_id_C2) * J[2](t_id_J);
                                    } // end loop theta_1

                                    local_operator(beta_0_1_2,alpha_0_1_2) = sum;

                                } // end loop alpha_0
                            } // end loop alpha_1
                        } // end loop alpha_2
                    } // end loop beta_0
                } // end loop beta_1
            } // end loop beta_2
            //--------------------------------------------------------------
        } //end if (!is_symmetric)
        else
        {

            //--------------------------------------------------------------
            for (Index beta_2 = 0 ; beta_2 < t_size_beta[2] ; ++beta_2)
            {
                t_id_J[2] = beta_2;

                const Index b_tmp_2 = beta_2*t_size_beta[1];

                Index beta_0_1 = 0;
                for (Index beta_1 = 0 ; beta_1 < t_size_beta[1] ; ++beta_1)
                {
                    const Index b_tmp_1 = (b_tmp_2 + beta_1) * t_size_beta[0];

                    for (Index beta_0 = 0 ; beta_0 < t_size_beta[0] ; ++beta_0, ++beta_0_1)
                    {
                        Index beta_0_1_2 = b_tmp_1 + beta_0 ;

                        t_id_C2[2] = beta_0_1;

                        for (Index alpha_2 = beta_2 ; alpha_2 < t_size_alpha[2] ; ++alpha_2)
                        {
                            t_id_J[1] = alpha_2;

                            const Index a_tmp_2 = alpha_2*t_size_alpha[1];

                            Index alpha_0_1 = 0;
                            for (Index alpha_1 = 0 ; alpha_1 < t_size_alpha[1] ; ++alpha_1)
                            {
                                const Index a_tmp_1 = (a_tmp_2 + alpha_1) * t_size_alpha[0];

                                for (Index alpha_0 = 0 ; alpha_0 < t_size_alpha[0] ; ++alpha_0,++alpha_0_1)
                                {
                                    const Index alpha_0_1_2 = a_tmp_1 + alpha_0 ;

                                    t_id_C2[1] = alpha_0_1;

                                    Real sum = 0.0;
                                    for (Index theta_2 = 0; theta_2 < t_size_theta[2] ; ++theta_2)
                                    {
                                        t_id_J [0] = theta_2;
                                        t_id_C2[0] = theta_2;

                                        sum += C2(t_id_C2) * J[2](t_id_J);
                                    } // end loop theta_1

                                    local_operator(beta_0_1_2,alpha_0_1_2) = sum;
                                } // end loop alpha_0
                            } // end loop alpha_1
                        } // end loop alpha_2
                    } // end loop beta_0
                } // end loop beta_1
            } // end loop beta_2

            // here we copy the upper triangular part of the matrix on the lower triangular part
            const Size n_basis_test = local_operator.get_num_rows();
            for (int test_id = 0 ; test_id < n_basis_test ; ++test_id)
                for (int trial_id = 0; trial_id < test_id ; ++trial_id)
                    local_operator(test_id,trial_id) = local_operator(trial_id,test_id);

            //--------------------------------------------------------------
        } // end if (is_symmetric)
    }
};
#endif






template<class PhysSpaceTest, class PhysSpaceTrial>
inline
auto
EllipticOperatorsSFIntegration<PhysSpaceTest,PhysSpaceTrial>::
evaluate_w_phi1Dtrial_phi1Dtest(
    const std::array<ValueTable<Real>,dim> &phi_1D_test,
    const std::array<ValueTable<Real>,dim> &phi_1D_trial,
    const TensorProductArray<dim> &quad_weights,
    const std::array<Real,dim> &length_element_edge) const -> std::array<DynamicMultiArray<Real,3>,dim>
{
    std::array<DynamicMultiArray<Real,3>,dim> moments;

    for (int dir = 0 ; dir < dim ; ++dir)
    {
        const auto &phi_test  = phi_1D_test [dir];
        const auto &phi_trial = phi_1D_trial[dir];
        const auto &w = quad_weights.get_data_direction(dir);

        const Size n_basis_test  = phi_test.get_num_functions();
        const Size n_basis_trial = phi_trial.get_num_functions();
        const Size n_pts = w.size();

        TensorSize<3> moments1D_tensor_size;
        moments1D_tensor_size[0] = n_pts;
        moments1D_tensor_size[1] = n_basis_trial;
        moments1D_tensor_size[2] = n_basis_test;

        auto &moments1D = moments[dir];
        moments1D.resize(moments1D_tensor_size);

        Assert(phi_test.get_num_points() == n_pts,
        ExcDimensionMismatch(phi_test.get_num_points(),n_pts));
        Assert(phi_trial.get_num_points() == n_pts,
        ExcDimensionMismatch(phi_trial.get_num_points(),n_pts));


        vector<Real> w_times_edge_length(n_pts);

        const Real edge_length = length_element_edge[dir];
        for (int jpt = 0 ; jpt < n_pts ; ++jpt)
            w_times_edge_length[jpt] = w[jpt] * edge_length;


        Index flat_id_I = 0 ;
        for (Index f_id_test = 0 ; f_id_test < n_basis_test ; ++f_id_test)
        {
            const auto phi_1D_test = phi_test.get_function_view(f_id_test);

            for (Index f_id_trial = 0 ; f_id_trial < n_basis_trial ; ++f_id_trial)
            {
                const auto phi_1D_trial = phi_trial.get_function_view(f_id_trial);

                for (int jpt = 0 ; jpt < n_pts ; ++jpt)
                    moments1D[flat_id_I++] =
                    w_times_edge_length[jpt] * phi_1D_test[jpt] * phi_1D_trial[jpt];
            } // end loop mu1
        } // end loop mu2
    } // end loop dir

    return moments;
}


template<class PhysSpaceTest, class PhysSpaceTrial>
inline
void
EllipticOperatorsSFIntegration<PhysSpaceTest,PhysSpaceTrial>::
eval_operator_u_v(
    const ElemTest &elem_test,
    const ElemTrial &elem_trial,
    const ValueVector<Real> &coeffs,
    DenseMatrix &operator_u_v) const
{
    //TODO: only the symmetric case is tested. In the non symmetric case, we need to check that
    // the physical space iterators have the same grid, map, reference space, index, etc.
    Assert(&elem_test == &elem_trial,ExcNotImplemented());



    //----------------------------------------------------
    // Assembly of the local mass matrix using sum-factorization -- begin
#ifdef TIME_PROFILING
    const TimePoint start_assembly_mass_matrix = Clock::now();
#endif //#ifdef TIME_PROFILING




    //--------------------------------------------------------------------------
    bool is_symmetric = this->test_if_same_space(elem_test,elem_trial);
    //--------------------------------------------------------------------------




    //--------------------------------------------------------------------------
#ifdef TIME_PROFILING
    const auto start_initialization = Clock::now();
#endif //#ifdef TIME_PROFILING




    //--------------------------------------------------------------------------
    // getting the number of basis along each coordinate direction for the test and trial space

    const Index comp = 0; // only scalar spaces for the moment

    // test space -- begin
    TensorIndex<dim> degree_test = elem_test.get_physical_space()->get_reference_space()->get_degree()[comp];
    TensorSize<dim> n_basis_elem_test(degree_test + 1);

    Assert(n_basis_elem_test.flat_size()==elem_test.get_num_basis(),
           ExcDimensionMismatch(n_basis_elem_test.flat_size(),elem_test.get_num_basis()));

    const auto weight_basis_test = MultiArrayUtils<dim>::compute_weight(n_basis_elem_test);

    const TensorSize<dim> n_basis_test = n_basis_elem_test;
    // test space -- end


    // trial space -- begin
    TensorIndex<dim> degree_trial = elem_trial.get_physical_space()->get_reference_space()->get_degree()[comp];
    TensorSize<dim> n_basis_elem_trial(degree_trial+1);

    Assert(n_basis_elem_trial.flat_size()==elem_trial.get_num_basis(),
           ExcDimensionMismatch(n_basis_elem_trial.flat_size(),elem_trial.get_num_basis()));

    const auto weight_basis_trial = MultiArrayUtils<dim>::compute_weight(n_basis_elem_trial);

    const TensorSize<dim> n_basis_trial = n_basis_elem_trial;
    // trial space -- end

//    const Size n_basis = n_basis_elem.flat_size();
    //--------------------------------------------------------------------------


    //--------------------------------------------------------------------------
    // getting the 1D values for the test space -- begin
    std::array< ValueTable<Real>,dim> phi_1D_test;
    {
        const auto &ref_elem_accessor = elem_test.get_ref_space_accessor();

        const auto &quad_points = ref_elem_accessor.get_quad_points();

        const auto phi_1D_test_table =
            ref_elem_accessor.evaluate_univariate_derivatives_at_points(0,quad_points);

        phi_1D_test = phi_1D_test_table[0]; // only valid for scalar spaces
    }
    // getting the 1D values for the test space -- end
    //--------------------------------------------------------------------------


    //--------------------------------------------------------------------------
    // getting the 1D values for the trial space -- begin
    std::array< ValueTable<Real>,dim> phi_1D_trial;
    {
        const auto &ref_elem_accessor = elem_trial.get_ref_space_accessor();

        const auto &quad_points = ref_elem_accessor.get_quad_points();

        const auto phi_1D_trial_table =
            ref_elem_accessor.evaluate_univariate_derivatives_at_points(0,quad_points);

        phi_1D_trial = phi_1D_trial_table[0]; // only valid for scalar spaces
    }
    // getting the 1D values for the trial space -- end
    //--------------------------------------------------------------------------


#ifdef TIME_PROFILING
    const auto end_initialization = Clock::now();
    const Duration elapsed_time_initialization = end_initialization - start_initialization;
    std::cout << "Elapsed_seconds initialization = " << elapsed_time_initialization.count() << std::endl;
#endif //#ifdef TIME_PROFILING
    //--------------------------------------------------------------------------



    //----------------------------------------------------
    // Coefficient evaluation phase -- begin
#ifdef TIME_PROFILING
    const TimePoint start_coefficient_evaluation = Clock::now();
#endif //#ifdef TIME_PROFILING


    // checks that the mapping used in the test space and in the trial space is the same
    Assert(elem_test.get_physical_space()->get_push_forward()->get_mapping() ==
           elem_trial.get_physical_space()->get_push_forward()->get_mapping(),
           ExcMessage("Test and trial spaces must have the same mapping (and the same grid)!"));


    // checks that the elements on the grid are the same
    Assert(static_cast<const CartesianGridElementAccessor<dim> &>(elem_test.get_ref_space_accessor()) ==
           static_cast<const CartesianGridElementAccessor<dim> &>(elem_trial.get_ref_space_accessor()),
           ExcMessage("Different elements for test space and trial space."));


    // performs the evaluation of the function coeffs*det(DF) at the quadrature points
    const auto &det_DF = elem_test.get_measures() ;
    Assert(det_DF.size() == coeffs.size(),
           ExcDimensionMismatch(det_DF.size(), coeffs.size()));
    Size n_points = det_DF.size();


    TensorSize<dim> n_points_1D = elem_test.get_ref_space_accessor().get_quad_points().get_num_points_direction();
    Assert(n_points_1D.flat_size() == n_points,
           ExcDimensionMismatch(n_points_1D.flat_size(),n_points));


    DynamicMultiArray<Real,dim> c_times_detDF(n_points_1D);
    for (Index ipt = 0 ; ipt < n_points ; ++ipt)
        c_times_detDF[ipt] = coeffs[ipt] * det_DF[ipt];


#ifdef TIME_PROFILING
    const TimePoint end_coefficient_evaluation = Clock::now();
    const Duration elapsed_time_coefficient_evaluation =
        end_coefficient_evaluation - start_coefficient_evaluation;
    std::cout << "Elapsed seconds coefficient evaluation mass= "
              << elapsed_time_coefficient_evaluation.count() << std::endl;
#endif //#ifdef TIME_PROFILING
    // Coefficient evaluation phase -- end
    //----------------------------------------------------





    //----------------------------------------------------
    // precalculation of the J[i](theta_i,alpha_i,beta_i) terms
    // (i.e. the weigths[theta_i] * phi_trial[alpha_i] * phi_test[beta_i] )
#ifdef TIME_PROFILING
    const auto start_compute_phi1Dtest_phi1Dtrial = Clock::now();
#endif //#ifdef TIME_PROFILING

    const std::array<Real,dim> length_element_edge =
        elem_test.get_ref_space_accessor().get_coordinate_lengths();

    const auto w_phi1Dtrial_phi1Dtest = evaluate_w_phi1Dtrial_phi1Dtest(
                                            phi_1D_test,
                                            phi_1D_trial,
                                            elem_test.get_ref_space_accessor().get_quad_points().get_weights(),
                                            length_element_edge);

#ifdef TIME_PROFILING
    const auto end_compute_phi1Dtest_phi1Dtrial = Clock::now();
    Duration elapsed_time_compute_phi1Dtest_phi1Dtrial =
        end_compute_phi1Dtest_phi1Dtrial- start_compute_phi1Dtest_phi1Dtrial;
    std::cout << "Elapsed seconds w * phi1d_trial * phi1d_test = "
              << elapsed_time_compute_phi1Dtest_phi1Dtrial.count() << std::endl;
#endif //#ifdef TIME_PROFILING
    //----------------------------------------------------




    //----------------------------------------------------
    // Assembly of the local mass matrix using sum-factorization -- begin
#ifdef TIME_PROFILING
    const auto start_sum_factorization = Clock::now();
#endif //#ifdef TIME_PROFILING

    TensorSize<3> tensor_size_C0;
    tensor_size_C0[0] = n_points_1D.flat_size(); // theta size
    tensor_size_C0[1] = 1; // alpha size
    tensor_size_C0[2] = 1; // beta size

    DynamicMultiArray<Real,3> C0(tensor_size_C0);
    const Size n_entries = tensor_size_C0.flat_size();

    Assert(n_entries == c_times_detDF.flat_size(),
           ExcDimensionMismatch(n_entries,c_times_detDF.flat_size()));
    for (Index entry_id = 0 ; entry_id < n_entries ; ++entry_id)
        C0[entry_id] = c_times_detDF[entry_id];


    SumFactorizationIntegrator<dim> integrate_sf;
    integrate_sf(is_symmetric,
                 n_points_1D,
                 n_basis_trial,
                 n_basis_test,
                 w_phi1Dtrial_phi1Dtest,
                 C0,
                 operator_u_v);


#ifdef TIME_PROFILING
    const auto end_sum_factorization = Clock::now();
    Duration elapsed_time_sum_factorization = end_sum_factorization - start_sum_factorization;
    std::cout << "Elapsed seconds sum-factorization = " << elapsed_time_sum_factorization.count() << std::endl;
#endif //#ifdef TIME_PROFILING
    // Assembly of the local mass matrix using sum-factorization -- end
    //----------------------------------------------------


#ifdef TIME_PROFILING
    const Duration elapsed_time_assemble = elapsed_time_sum_factorization +
                                           elapsed_time_compute_phi1Dtest_phi1Dtrial +
                                           elapsed_time_coefficient_evaluation +
                                           elapsed_time_initialization ;
    std::cout << "Elapsed seconds assemblying = " << elapsed_time_assemble.count() << std::endl;


    const TimePoint end_assembly_mass_matrix = Clock::now();

    const_cast<Duration &>(this->elapsed_time_operator_u_v_) += end_assembly_mass_matrix - start_assembly_mass_matrix;
    std::cout << "Elapsed seconds operator u_v sum-factorization= "
              << this->elapsed_time_operator_u_v_.count() << std::endl;
#endif //#ifdef TIME_PROFILING
    // Assembly of the local mass matrix using sum-factorization -- end
    //----------------------------------------------------
}



template<class PhysSpaceTest, class PhysSpaceTrial>
inline
void
EllipticOperatorsSFIntegration<PhysSpaceTest,PhysSpaceTrial>::
eval_operator_gradu_gradv(
    const ElemTest &elem_test,
    const ElemTrial &elem_trial,
    const vector<TMatrix<space_dim,space_dim>> &coeffs,
    DenseMatrix &operator_gradu_gradv) const
{

    //----------------------------------------------------
    // Assembly of the local stiffness matrix using sum-factorization -- begin
#ifdef TIME_PROFILING
    const TimePoint start_assembly_stiffness_matrix = Clock::now();

    //--------------------------------------------------------------------------
    const auto start_initialization = Clock::now();
#endif //#ifdef TIME_PROFILING




    //--------------------------------------------------------------------------
    // getting the number of basis along each coordinate direction for the test and trial space

    const Index comp = 0; // only scalar spaces for the moment

    // test space -- begin
    TensorIndex<dim> degree_test = elem_test.get_physical_space()->get_reference_space()->get_degree()[comp];
    TensorSize<dim> n_basis_elem_test(degree_test + 1);

    Assert(n_basis_elem_test.flat_size()==elem_test.get_num_basis(),
           ExcDimensionMismatch(n_basis_elem_test.flat_size(),elem_test.get_num_basis()));

    const auto weight_basis_test = MultiArrayUtils<dim>::compute_weight(n_basis_elem_test);

    const TensorSize<dim> n_basis_test = n_basis_elem_test;
    // test space -- end


    // trial space -- begin
    TensorIndex<dim> degree_trial = elem_trial.get_physical_space()->get_reference_space()->get_degree()[comp];
    TensorSize<dim> n_basis_elem_trial(degree_trial + 1);

    Assert(n_basis_elem_trial.flat_size()==elem_trial.get_num_basis(),
           ExcDimensionMismatch(n_basis_elem_trial.flat_size(),elem_trial.get_num_basis()));

    const auto weight_basis_trial = MultiArrayUtils<dim>::compute_weight(n_basis_elem_trial);

    const TensorSize<dim> n_basis_trial = n_basis_elem_trial;
    // trial space -- end

//    const Size n_basis = n_basis_elem.flat_size();
    //--------------------------------------------------------------------------


    //--------------------------------------------------------------------------
    // getting the 1D values for the test space -- begin
    std::array< ValueTable<Real>,dim> phi_1D_test;
    std::array< ValueTable<Real>,dim> grad_phi_1D_test;
    {
        const auto &ref_elem_accessor = elem_test.get_ref_space_accessor();

        const auto &quad_points = ref_elem_accessor.get_quad_points();

        const auto phi_1D_test_table =
            ref_elem_accessor.evaluate_univariate_derivatives_at_points(0,quad_points);

        const auto grad_phi_1D_test_table =
            ref_elem_accessor.evaluate_univariate_derivatives_at_points(1,quad_points);

        phi_1D_test = phi_1D_test_table[0]; // only valid for scalar spaces
        grad_phi_1D_test = grad_phi_1D_test_table[0];  // only valid for scalar spaces
    }
    // getting the 1D values for the test space -- end
    //--------------------------------------------------------------------------


    //--------------------------------------------------------------------------
    // getting the 1D values for the trial space -- begin
    std::array< ValueTable<Real>,dim> phi_1D_trial;
    std::array< ValueTable<Real>,dim> grad_phi_1D_trial;
    {
        const auto &ref_elem_accessor = elem_trial.get_ref_space_accessor();

        const auto &quad_points = ref_elem_accessor.get_quad_points();

        const auto phi_1D_trial_table =
            ref_elem_accessor.evaluate_univariate_derivatives_at_points(0,quad_points);

        const auto grad_phi_1D_trial_table =
            ref_elem_accessor.evaluate_univariate_derivatives_at_points(1,quad_points);

        phi_1D_trial = phi_1D_trial_table[0]; // only valid for scalar spaces
        grad_phi_1D_trial = grad_phi_1D_trial_table[0];  // only valid for scalar spaces
    }
    // getting the 1D values for the trial space -- end
    //--------------------------------------------------------------------------

#ifdef TIME_PROFILING
    const auto end_initialization = Clock::now();
    const Duration elapsed_time_initialization = end_initialization - start_initialization;
    std::cout << "Elapsed_seconds initialization = " << elapsed_time_initialization.count() << std::endl;
#endif //#ifdef TIME_PROFILING
    //--------------------------------------------------------------------------






    //----------------------------------------------------
    // Coefficient evaluation phase -- begin
#ifdef TIME_PROFILING
    const TimePoint start_coefficient_evaluation = Clock::now();
#endif //#ifdef TIME_PROFILING


    // checks that the mapping used in the test space and in the trial space is the same
    Assert(elem_test.get_physical_space()->get_push_forward()->get_mapping() ==
           elem_trial.get_physical_space()->get_push_forward()->get_mapping(),
           ExcMessage("Test and trial spaces must have the same mapping (and the same grid)!"));


    // checks that the elements on the grid are the same
    Assert(static_cast<const CartesianGridElementAccessor<dim> &>(elem_test.get_ref_space_accessor()) ==
           static_cast<const CartesianGridElementAccessor<dim> &>(elem_trial.get_ref_space_accessor()),
           ExcMessage("Different elements for test space and trial space."));


    // performs the evaluation of the function DF^{-1} * C * DF^{-T} * det(DF) at the quadrature points
    const auto &det_DF = elem_test.get_measures() ;

    const auto &invDF = elem_test.get_push_forward_accessor().get_inv_gradients();

    Size n_points = coeffs.size();
    Assert(det_DF.size() == n_points,
           ExcDimensionMismatch(det_DF.size(), n_points));
    Assert(invDF.size() == n_points,
           ExcDimensionMismatch(invDF.size(), n_points));


    TensorSize<dim> n_points_1D = elem_test.get_ref_space_accessor().get_quad_points().get_num_points_direction();
    Assert(n_points_1D.flat_size() == n_points,
           ExcDimensionMismatch(n_points_1D.flat_size(),n_points));

//    LogStream out;
    vector<TMatrix<dim,dim>> C_hat(n_points);
    for (Index ipt = 0 ; ipt < n_points ; ++ipt)
    {
        TMatrix<dim,dim> &C_hat_ipt = C_hat[ipt];

        const TMatrix<space_dim,space_dim> C_ipt = coeffs[ipt] * det_DF[ipt];
        const auto &invDF_ipt = invDF[ipt];
//        out << "C at point " << ipt << "= " << C_ipt << endl;

        for (Index r = 0 ; r < dim ; ++r)
            for (Index s = 0 ; s < dim ; ++s)
                for (Index k = 0 ; k < space_dim ; ++k)
                    for (Index l = 0 ; l < space_dim ; ++l)
                        C_hat_ipt[r][s] += invDF_ipt[k][r](0) * C_ipt[k][l] * invDF_ipt[l][s](0);


//        out << "C_hat at point " << ipt << "= " << C_hat_ipt << endl;
    }


#ifdef TIME_PROFILING
    const TimePoint end_coefficient_evaluation = Clock::now();
    const Duration elapsed_time_coefficient_evaluation =
        end_coefficient_evaluation - start_coefficient_evaluation;
    std::cout << "Elapsed seconds coefficient evaluation stiffness = "
              << elapsed_time_coefficient_evaluation.count() << std::endl;
#endif //#ifdef TIME_PROFILING
    // Coefficient evaluation phase -- end
    //----------------------------------------------------






    //----------------------------------------------------
    // Assembly of the local stiffness matrix using sum-factorization -- begin
#ifdef TIME_PROFILING
    Duration elapsed_time_compute_phi1Dtest_phi1Dtrial;

    const auto start_sum_factorization = Clock::now();
#endif //#ifdef TIME_PROFILING

    TensorSize<3> tensor_size_C0;
    tensor_size_C0[0] = n_points_1D.flat_size(); // theta size
    tensor_size_C0[1] = 1; // alpha size
    tensor_size_C0[2] = 1; // beta size

    DynamicMultiArray<Real,3> C0(tensor_size_C0);
    const Size n_entries = tensor_size_C0.flat_size();

    DenseMatrix operator_gradu_gradv_tmp(n_basis_test.flat_size(),n_basis_trial.flat_size());
    operator_gradu_gradv.clear();

    std::array<DynamicMultiArray<Real,3>,dim> J;

    std::array< ValueTable<Real>,dim> trial_1D;
    std::array< ValueTable<Real>,dim> test_1D;
    for (Index k = 0 ; k < dim ; ++k)
    {
        for (Index l = 0 ; l < dim ; ++l)
        {
            for (Index i = 0 ; i < dim ; ++i)
            {
                if (i == k && i == l)
                {
                    trial_1D[i] = grad_phi_1D_trial[i];
                    test_1D [i] = grad_phi_1D_test [i];
                }
                else if (i == k && i != l)
                {
                    trial_1D[i] = grad_phi_1D_trial[i];
                    test_1D [i] =      phi_1D_test [i];
                }
                else if (i != k && i == l)
                {
                    trial_1D[i] =      phi_1D_trial[i];
                    test_1D [i] = grad_phi_1D_test [i];
                }
                else if (i != k && i != l)
                {
                    trial_1D[i] = phi_1D_trial[i];
                    test_1D [i] = phi_1D_test [i];
                }
            } // end loop i


            //----------------------------------------------------
            // precalculation of the J[i](theta_i,alpha_i,beta_i) terms
            // (i.e. the weigths[theta_i] * phi_trial[alpha_i] * phi_test[beta_i] )
#ifdef TIME_PROFILING
            const auto start_compute_phi1Dtest_phi1Dtrial = Clock::now();
#endif //#ifdef TIME_PROFILING

            const std::array<Real,dim> length_element_edge =
                elem_test.get_ref_space_accessor().get_coordinate_lengths();

            const auto J = evaluate_w_phi1Dtrial_phi1Dtest(
                               test_1D,
                               trial_1D,
                               elem_test.get_ref_space_accessor().get_quad_points().get_weights(),
                               length_element_edge);

#ifdef TIME_PROFILING
            const auto end_compute_phi1Dtest_phi1Dtrial = Clock::now();
            elapsed_time_compute_phi1Dtest_phi1Dtrial +=
                end_compute_phi1Dtest_phi1Dtrial- start_compute_phi1Dtest_phi1Dtrial;
#endif //#ifdef TIME_PROFILING
            //----------------------------------------------------



            Assert(n_entries == C_hat.size(),
                   ExcDimensionMismatch(n_entries,C_hat.size()));
            for (Index entry_id = 0 ; entry_id < n_entries ; ++entry_id)
                C0[entry_id] = C_hat[entry_id][k][l];

            operator_gradu_gradv_tmp.clear();


            SumFactorizationIntegrator<dim> integrate_sf;
            integrate_sf(false, //non symmetric
                         n_points_1D,
                         n_basis_trial,
                         n_basis_test,
                         J,
                         C0,
                         operator_gradu_gradv_tmp);

            operator_gradu_gradv += operator_gradu_gradv_tmp;

        }

    }



#ifdef TIME_PROFILING
    const auto end_sum_factorization = Clock::now();

    std::cout << "Elapsed seconds w * trial * test = "
              << elapsed_time_compute_phi1Dtest_phi1Dtrial.count() << std::endl;

    Duration elapsed_time_sum_factorization = end_sum_factorization - start_sum_factorization;
    std::cout << "Elapsed seconds sum-factorization = " << elapsed_time_sum_factorization.count() << std::endl;
    // Assembly of the local stiffness matrix using sum-factorization -- end
    //----------------------------------------------------


    const Duration elapsed_time_assemble = elapsed_time_sum_factorization +
                                           elapsed_time_coefficient_evaluation +
                                           elapsed_time_initialization ;
    std::cout << "Elapsed seconds assemblying = " << elapsed_time_assemble.count() << std::endl;


    const TimePoint end_assembly_stiffness_matrix = Clock::now();

    const_cast<Duration &>(this->elapsed_time_operator_gradu_gradv_)
        = end_assembly_stiffness_matrix - start_assembly_stiffness_matrix;
    std::cout << "Elapsed seconds operator gradu_gradv sum-factorization= "
              << this->elapsed_time_operator_gradu_gradv_.count() << std::endl;
#endif //#ifdef TIME_PROFILING
    // Assembly of the local stiffness matrix using sum-factorization -- end
    //----------------------------------------------------


//    Assert(false,ExcNotImplemented());
//    AssertThrow(false,ExcNotImplemented());
}


template <class PhysSpaceTest,class PhysSpaceTrial>
inline
void
EllipticOperatorsSFIntegration<PhysSpaceTest,PhysSpaceTrial>::
eval_operator_rhs_v(
    const ElemTest &elem_test,
    const ValueVector<typename PhysSpaceTrial::Value> &f,
    DenseVector &operator_rhs_v) const
{
    //TODO: (martinelli 22 Sep 2014): this function is not implemented using sum_factorization. Fix it!
    const Size n_basis_test  = elem_test .get_num_basis();

    const auto &phi_test  = elem_test.get_basis_values();
    const auto &w_meas  = elem_test.get_w_measures();


    Assert(operator_rhs_v.size() == n_basis_test,
           ExcDimensionMismatch(operator_rhs_v.size(),n_basis_test));


    const Size n_qp = f.get_num_points();
    Assert(n_qp == phi_test.get_num_points(),ExcDimensionMismatch(n_qp,phi_test.get_num_points()));

    vector<Real> f_times_w_meas(n_qp);
    for (int qp = 0; qp < n_qp; ++qp)
        f_times_w_meas[qp] = f[qp](0) * w_meas[qp];


    operator_rhs_v.clear();
    for (int i = 0; i < n_basis_test; ++i)
    {
        const auto phi_i = phi_test.get_function_view(i);
        for (int qp = 0; qp < n_qp; ++qp)
            operator_rhs_v(i) += phi_i[qp](0) * f_times_w_meas[qp];
    }
}


template <class PhysSpaceTest,class PhysSpaceTrial>
inline
void
EllipticOperatorsSFIntegration<PhysSpaceTest,PhysSpaceTrial>::
eval_operator_gradu_v(
    const ElemTest &elem_test,
    const ElemTrial &elem_trial,
    const ValueVector<typename PhysSpaceTrial::Gradient> &beta,
    DenseMatrix &operator_gradu_v) const
{
    Assert(false,ExcNotImplemented());
    AssertThrow(false,ExcNotImplemented());
}





IGA_NAMESPACE_CLOSE


#endif // #ifndef ELLIPTIC_OPERATORS_SF_INTEGRATION_H_
