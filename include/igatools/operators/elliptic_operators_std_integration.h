//-+--------------------------------------------------------------------
// Igatools a general purpose Isogeometric analysis library.
// Copyright (C) 2012-2016  by the igatools authors (see authors.txt).
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

#ifndef ELLIPTIC_OPERATORS_STD_INTEGRATION_H_
#define ELLIPTIC_OPERATORS_STD_INTEGRATION_H_

#include <igatools/base/config.h>
#include <igatools/operators/elliptic_operators.h>

#include <vector>

IGA_NAMESPACE_OPEN

#if 0
/**
 * @brief Class containing the methods for the evaluation of some elliptic operators
 * (mass matrix, stiffness matrix, etc.) on a Bezier element
 * using the <em>standard quadrature approach</em>.
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
class EllipticOperatorsStdIntegration :
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
    EllipticOperatorsStdIntegration<PhysSpaceTest,PhysSpaceTrial> &
    operator=(const EllipticOperatorsStdIntegration<PhysSpaceTest,PhysSpaceTrial> &in) = default;


    /** Move assignment operator. */
    EllipticOperatorsStdIntegration<PhysSpaceTest,PhysSpaceTrial> &
    operator=(EllipticOperatorsStdIntegration<PhysSpaceTest,PhysSpaceTrial> &&in) = default;
    ///@}


    /**
     * This function evaluates the local (i.e. element-based) matrix \f$ A_e \f$
     * for witch its entries are
     * \f[
          (A_e)_{ij} = \int_{\Omega_e} \phi^{e,\text{test}}_i(x) C(x) \phi^{e,\text{trial}}_j(x) \; d \Omega.
       \f]
     * The matrix \f$ A_e \f$ is commonly referred as <em>local mass-matrix</em>.
     */
    void eval_operator_u_v(
        const ElemTest &elem_test,
        const ElemTrial &elem_trial,
        const ValueVector<Real> &c,
        DenseMatrix &operator_u_v) const;

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
     void eval_operator_gradu_gradv(
        const ElemTest &elem_test,
        const ElemTrial &elem_trial,
        const vector<TMatrix<space_dim,space_dim>> &coeffs,
        DenseMatrix &operator_gradu_gradv) const;

    /**
     * This function evaluates the local (i.e. element-based) vector \f$ f_e \f$
     * for witch its entries are
     * \f[
          (f_e)_{i} = \int_{\Omega_e} \phi^{e,\text{test}}_i
          f(x)  \; d \Omega.
       \f]
     */
    void eval_operator_rhs_v(
        const ElemTest &elem_test,
        const ValueVector<typename PhysSpaceTrial::Value> &f,
        DenseVector &operator_rhs_v) const;

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
    void eval_operator_gradu_v(
        const ElemTest &elem_test,
        const ElemTrial &elem_trial,
        const ValueVector<typename PhysSpaceTrial::Gradient> &beta,
        DenseMatrix &operator_gradu_v) const;

};


template <class PhysSpaceTest,class PhysSpaceTrial>
inline
void
EllipticOperatorsStdIntegration<PhysSpaceTest,PhysSpaceTrial>::
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
    // Assembly of the local mass matrix using the standard quadrature -- begin
#ifdef TIME_PROFILING
    const auto start_assembly_mass_matrix = Clock::now();
#endif // #ifdef TIME_PROFILING



    const bool is_symmetric = this->test_if_same_space(elem_test,elem_trial);

    const Size n_basis_test  = elem_test .get_num_basis();
    const Size n_basis_trial = elem_trial.get_num_basis();

    const auto &phi_test  = elem_test.get_basis_values();
    const auto &phi_trial = elem_trial.get_basis_values();
    const auto &w_meas  = elem_test.get_w_measures();


    Assert(operator_u_v.get_num_rows() == n_basis_test,
           ExcDimensionMismatch(operator_u_v.get_num_rows(),n_basis_test));
    Assert(operator_u_v.get_num_cols() == n_basis_trial,
           ExcDimensionMismatch(operator_u_v.get_num_cols(),n_basis_trial));


    const Size n_qp = coeffs.size();
    Assert(n_qp == phi_test.get_num_points(),ExcDimensionMismatch(n_qp,phi_test.get_num_points()));
    Assert(n_qp == phi_trial.get_num_points(),ExcDimensionMismatch(n_qp,phi_trial.get_num_points()));

    vector<Real> coeffs_times_w_meas(n_qp);
    for (int qp = 0; qp < n_qp; ++qp)
        coeffs_times_w_meas[qp] = coeffs[qp] * w_meas[qp];


    operator_u_v.clear();

    if (!is_symmetric)
    {
        for (int i = 0; i < n_basis_test; ++i)
        {
            const auto phi_i = phi_test.get_function_view(i);
            for (int j = 0; j < n_basis_trial; ++j)
            {
                const auto phi_j = phi_trial.get_function_view(j);
                for (int qp = 0; qp < n_qp; ++qp)
                    operator_u_v(i,j) += phi_j[qp](0) * (phi_i[qp](0) * coeffs_times_w_meas[qp]);
            }
        }

    } // end if (!is_symmetric)
    else
    {
        for (int i = 0; i < n_basis_test; ++i)
        {
            const auto phi_i = phi_test.get_function_view(i);
            for (int j = i; j < n_basis_trial; ++j)
            {
                const auto phi_j = phi_trial.get_function_view(j);
                for (int qp = 0; qp < n_qp; ++qp)
                    operator_u_v(i,j) += phi_j[qp](0) * (phi_i[qp](0) * coeffs_times_w_meas[qp]);
            }
        }
        for (int i = 0; i < n_basis_test; ++i)
            for (int j = 0; j < i; ++j)
                operator_u_v(i,j) = operator_u_v(j,i);
    } // end if (is_symmetric)



#ifdef TIME_PROFILING
    const auto end_assembly_mass_matrix = Clock::now();
    const_cast<Duration &>(this->elapsed_time_operator_u_v_) = end_assembly_mass_matrix - start_assembly_mass_matrix;
    std::cout << "Elapsed seconds operator u_v standard quadrature= "
              << this->elapsed_time_operator_u_v_.count() << std::endl;
#endif // #ifdef TIME_PROFILING

    // Assembly of the local mass matrix using the standard quadrature -- begin
    //----------------------------------------------------

}


template <class PhysSpaceTest,class PhysSpaceTrial>
inline
void
EllipticOperatorsStdIntegration<PhysSpaceTest,PhysSpaceTrial>::
eval_operator_gradu_gradv(
    const ElemTest &elem_test,
    const ElemTrial &elem_trial,
    const vector<TMatrix<space_dim,space_dim>> &coeffs,
    DenseMatrix &operator_gradu_gradv) const
{
    //TODO: only the symmetric case is tested. In the non symmetric case, we need to check that
    // the physical space iterators have the same grid, map, reference space, index, etc.
    Assert(&elem_test == &elem_trial,ExcNotImplemented());

    //----------------------------------------------------
    // Assembly of the local stiffness matrix using the standard quadrature -- begin
#ifdef TIME_PROFILING
    const TimePoint start_assembly_stiffness_matrix = Clock::now();
#endif // #ifdef TIME_PROFILING




    const Size n_basis_test  = elem_test .get_num_basis();
    const Size n_basis_trial = elem_trial.get_num_basis();

    const auto &grad_phi_test  = elem_test.get_basis_gradients();
    const auto &grad_phi_trial = elem_trial.get_basis_gradients();
    const auto &w_meas  = elem_test.get_w_measures();


    Assert(operator_gradu_gradv.get_num_rows() == n_basis_test,
           ExcDimensionMismatch(operator_gradu_gradv.get_num_rows(),n_basis_test));
    Assert(operator_gradu_gradv.get_num_cols() == n_basis_trial,
           ExcDimensionMismatch(operator_gradu_gradv.get_num_cols(),n_basis_trial));


    const Size n_qp = coeffs.size();
    Assert(n_qp == grad_phi_test.get_num_points(),
           ExcDimensionMismatch(n_qp,grad_phi_test.get_num_points()));
    Assert(n_qp == grad_phi_trial.get_num_points(),
           ExcDimensionMismatch(n_qp,grad_phi_trial.get_num_points()));


//    LogStream out;
    vector<TMatrix<space_dim,space_dim>> coeffs_times_w_meas(n_qp);
    for (Index qp = 0; qp < n_qp; ++qp)
    {
        for (Index i = 0 ; i < PhysSpaceTest::space_dim ; ++i)
            for (Index j = 0 ; j < PhysSpaceTrial::space_dim ; ++j)
                coeffs_times_w_meas[qp][i][j] = coeffs[qp][i][j] * w_meas[qp];

//        out << "Coeffs at point " << qp <<"=   " << coeffs_times_w_meas[qp] <<endl;
    }

    // type of the gradients of the basis functions in the test space
    using grad_test_t = Derivatives<PhysSpaceTest::space_dim,PhysSpaceTest::range,PhysSpaceTest::rank,1>;

    operator_gradu_gradv.clear();
    for (Index i = 0; i < n_basis_test; ++i)
    {
        const auto grad_phi_i = grad_phi_test.get_function_view(i);

        vector<grad_test_t> coeffs_times_grad_phi_i(n_qp);
        for (Index qp = 0; qp < n_qp; ++qp)
        {
            const auto &C = coeffs_times_w_meas[qp];
            auto &C_grad_phi_test = coeffs_times_grad_phi_i[qp];
            for (Index i = 0 ; i < PhysSpaceTest::space_dim; ++i)
            {
                C_grad_phi_test[i] = 0.0;
                for (Index j = 0 ; j < PhysSpaceTest::space_dim; ++j)
                {
                    //TODO: check if C[i][j] or C[j][i]
                    C_grad_phi_test[i] += C[i][j] * grad_phi_i[qp][j];
                }
            }
        }

        for (Index j = 0; j < n_basis_trial; ++j)
        {
            const auto grad_phi_j = grad_phi_trial.get_function_view(j);
            for (Index qp = 0; qp < n_qp; ++qp)
                operator_gradu_gradv(i,j) += scalar_product(coeffs_times_grad_phi_i[qp],grad_phi_j[qp]);
        }
    }



#ifdef TIME_PROFILING
    const TimePoint end_assembly_stiffness_matrix = Clock::now();
    const_cast<Duration &>(this->elapsed_time_operator_gradu_gradv_)
        = end_assembly_stiffness_matrix - start_assembly_stiffness_matrix;
    std::cout << "Elapsed seconds operator gradu_gradv standard quadrature= "
              << this->elapsed_time_operator_gradu_gradv_.count() << std::endl;
#endif // #ifdef TIME_PROFILING

    // Assembly of the local mass matrix using the standard quadrature -- begin
    //----------------------------------------------------
}


template <class PhysSpaceTest,class PhysSpaceTrial>
inline
void
EllipticOperatorsStdIntegration<PhysSpaceTest,PhysSpaceTrial>::
eval_operator_rhs_v(
    const ElemTest &elem_test,
    const ValueVector<typename PhysSpaceTrial::Value> &f,
    DenseVector &operator_rhs_v) const
{
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
EllipticOperatorsStdIntegration<PhysSpaceTest,PhysSpaceTrial>::
eval_operator_gradu_v(
    const ElemTest &elem_test,
    const ElemTrial &elem_trial,
    const ValueVector<typename PhysSpaceTrial::Gradient> &beta,
    DenseMatrix &operator_gradu_v) const
{
    //TODO: only the symmetric case is tested. In the non symmetric case, we need to check that
    // the physical space iterators have the same grid, map, reference space, index, etc.
    Assert(&elem_test == &elem_trial,ExcNotImplemented());

    //----------------------------------------------------
    // Assembly of the local stiffness matrix using the standard quadrature -- begin
#ifdef TIME_PROFILING
    const TimePoint start_assembly_convection_matrix = Clock::now();
#endif // #ifdef TIME_PROFILING




    const Size n_basis_test  = elem_test .get_num_basis();
    const Size n_basis_trial = elem_trial.get_num_basis();

    const auto      &phi_test  = elem_test.get_basis_values();
    const auto &grad_phi_trial = elem_trial.get_basis_gradients();
    const auto &w_meas  = elem_test.get_w_measures();


    Assert(operator_gradu_v.get_num_rows() == n_basis_test,
           ExcDimensionMismatch(operator_gradu_v.get_num_rows(),n_basis_test));
    Assert(operator_gradu_v.get_num_cols() == n_basis_trial,
           ExcDimensionMismatch(operator_gradu_v.get_num_cols(),n_basis_trial));


    const Size n_qp = beta.get_num_points();
    Assert(n_qp == phi_test.get_num_points(),
           ExcDimensionMismatch(n_qp,phi_test.get_num_points()));
    Assert(n_qp == grad_phi_trial.get_num_points(),
           ExcDimensionMismatch(n_qp,grad_phi_trial.get_num_points()));


//    LogStream out;
    vector<typename PhysSpaceTrial::Gradient> beta_times_w_meas(n_qp);
    for (Index qp = 0; qp < n_qp; ++qp)
        beta_times_w_meas[qp] = beta[qp] * w_meas[qp];


    ValueTable<typename PhysSpaceTrial::Value> beta_dot_grad_trial(n_basis_trial,n_qp);
    for (Index i = 0; i < n_basis_trial; ++i)
    {
        const auto grad_trial_i = grad_phi_trial.get_function_view(i);
        auto beta_dot_grad_trial_i = beta_dot_grad_trial.get_function_view(i);

        for (Index qp = 0; qp < n_qp; ++qp)
            beta_dot_grad_trial_i[qp] = scalar_product(beta_times_w_meas[qp], grad_trial_i[qp]) ;
    }

    operator_gradu_v.clear();
    for (Index i = 0; i < n_basis_test; ++i)
    {
        const auto phi_i = phi_test.get_function_view(i);

        for (Index j = 0; j < n_basis_trial; ++j)
        {
            auto beta_dot_grad_trial_j = beta_dot_grad_trial.get_function_view(j);

            for (Index qp = 0; qp < n_qp; ++qp)
                operator_gradu_v(i,j) += scalar_product(phi_i[qp], beta_dot_grad_trial_j[qp]);
        }
    }


#ifdef TIME_PROFILING
    const TimePoint end_assembly_convection_matrix = Clock::now();
    const_cast<Duration &>(this->elapsed_time_operator_gradu_v_)
        = end_assembly_convection_matrix - start_assembly_convection_matrix;
    std::cout << "Elapsed seconds operator gradu_gradv standard quadrature= "
              << this->elapsed_time_operator_gradu_v_.count() << std::endl;
#endif // #ifdef TIME_PROFILING

    // Assembly of the local mass matrix using the standard quadrature -- begin
    //----------------------------------------------------
}
#endif

IGA_NAMESPACE_CLOSE


#endif // #ifndef ELLIPTIC_OPERATORS_STD_INTEGRATION_H_
