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

#include <igatools/operators/integrator_sum_factorization.h>
#include <igatools/base/exceptions.h>


IGA_NAMESPACE_OPEN







template <int dim, int r>
void
IntegratorSumFactorization<dim,r>::
operator()(
  const bool is_symmetric,
  const TensorSize<dim> &t_size_theta,
  const TensorSize<dim> &t_size_alpha,
  const TensorSize<dim> &t_size_beta,
  const SafeSTLArray<DynamicMultiArray<Real,3>,dim> &J,
  const DynamicMultiArray<Real,3> &Cpre,
  const int row_id_begin,
  const int row_id_last,
  const int col_id_begin,
  const int col_id_last,
  DenseMatrix &local_operator) const
{
  Assert(false,ExcNotImplemented());
  AssertThrow(false,ExcNotImplemented());
}




template <int dim>
void
IntegratorSumFactorization<dim,1>::
operator()(
  const bool is_symmetric,
  const TensorSize<dim> &t_size_theta,
  const TensorSize<dim> &t_size_alpha,
  const TensorSize<dim> &t_size_beta,
  const SafeSTLArray<DynamicMultiArray<Real,3>,dim> &J,
  const DynamicMultiArray<Real,3> &Cpre,
  const int row_id_begin,
  const int row_id_last,
  const int col_id_begin,
  const int col_id_last,
  DenseMatrix &local_operator) const
{
  Assert(false,ExcNotImplemented());
  AssertThrow(false,ExcNotImplemented());
}



void
IntegratorSumFactorization<0,0>::
operator()(
  const bool is_symmetric,
  const TensorSize<0> &t_size_theta,
  const TensorSize<0> &t_size_alpha,
  const TensorSize<0> &t_size_beta,
  const SafeSTLArray<DynamicMultiArray<Real,3>,0> &J,
  const DynamicMultiArray<Real,3> &Cpre,
  const int row_id_begin,
  const int row_id_last,
  const int col_id_begin,
  const int col_id_last,
  DenseMatrix &local_operator) const
{
  Assert(false,ExcNotImplemented());
  AssertThrow(false,ExcNotImplemented());
}


void
IntegratorSumFactorization<1,1>::
operator()(
  const bool is_symmetric,
  const TensorSize<1> &t_size_theta,
  const TensorSize<1> &t_size_alpha,
  const TensorSize<1> &t_size_beta,
  const SafeSTLArray<DynamicMultiArray<Real,3>,1> &J,
  const DynamicMultiArray<Real,3> &C,
  const int row_id_begin,
  const int row_id_last,
  const int col_id_begin,
  const int col_id_last,
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

void
IntegratorSumFactorization<2,2>::
operator()(
  const bool is_symmetric,
  const TensorSize<2> &t_size_theta,
  const TensorSize<2> &t_size_alpha,
  const TensorSize<2> &t_size_beta,
  const SafeSTLArray<DynamicMultiArray<Real,3>,2> &J,
  const DynamicMultiArray<Real,3> &C,
  const int row_id_begin,
  const int row_id_last,
  const int col_id_begin,
  const int col_id_last,
  DenseMatrix &local_operator) const
{
  const int n_rows = row_id_last - row_id_begin + 1;
  Assert(n_rows >= 1,ExcLowerRange(n_rows,1));
  Assert(n_rows <= local_operator.get_num_rows(),
         ExcUpperRange(n_rows,local_operator.get_num_rows()));
  Assert(t_size_beta.flat_size() == n_rows,
         ExcDimensionMismatch(t_size_beta.flat_size(),n_rows));
  Assert(row_id_begin >= 0,ExcLowerRange(row_id_begin,0));
  Assert(row_id_last <= (local_operator.get_num_rows()-1),
         ExcUpperRange(row_id_last,local_operator.get_num_rows()-1));


  const int n_cols = col_id_last - col_id_begin + 1;
  Assert(n_cols >= 1,ExcLowerRange(n_cols,1));
  Assert(n_cols <= local_operator.get_num_cols(),
         ExcUpperRange(n_cols,local_operator.get_num_cols()));
  Assert(t_size_alpha.flat_size() == n_cols,
         ExcDimensionMismatch(t_size_alpha.flat_size(),n_cols));
  Assert(col_id_begin >= 0,ExcLowerRange(col_id_begin,0));
  Assert(col_id_last <= (local_operator.get_num_cols()-1),
         ExcUpperRange(col_id_last,local_operator.get_num_cols()-1));





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
        Index beta_0_1 = beta_1*t_size_beta[0] + beta_0 + row_id_begin;

        t_id_C1[2] = beta_0;

        for (Index alpha_1 = 0 ; alpha_1 < t_size_alpha[1] ; ++alpha_1)
        {
          t_id_J[1] = alpha_1;
          for (Index alpha_0 = 0 ; alpha_0 < t_size_alpha[0] ; ++alpha_0)
          {
            Index alpha_0_1 = alpha_1*t_size_alpha[0] + alpha_0 + col_id_begin;

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
        Index beta_0_1 = beta_1*t_size_beta[0] + beta_0 + row_id_begin;

        t_id_C1[2] = beta_0;

        for (Index alpha_1 = beta_1 ; alpha_1 < t_size_alpha[1] ; ++alpha_1)
        {
          t_id_J[1] = alpha_1;
          for (Index alpha_0 = 0 ; alpha_0 < t_size_alpha[0] ; ++alpha_0)
          {
            Index alpha_0_1 = alpha_1*t_size_alpha[0] + alpha_0 + col_id_begin;

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


void
IntegratorSumFactorization<3,3>::
operator()(
  const bool is_symmetric,
  const TensorSize<3> &t_size_theta,
  const TensorSize<3> &t_size_alpha,
  const TensorSize<3> &t_size_beta,
  const SafeSTLArray<DynamicMultiArray<Real,3>,3> &J,
  const DynamicMultiArray<Real,3> &C,
  const int row_id_begin,
  const int row_id_last,
  const int col_id_begin,
  const int col_id_last,
  DenseMatrix &local_operator) const
{
  const int n_rows = row_id_last - row_id_begin + 1;
  Assert(n_rows >= 1,ExcLowerRange(n_rows,1));
  Assert(n_rows <= local_operator.get_num_rows(),
         ExcUpperRange(n_rows,local_operator.get_num_rows()));
  Assert(t_size_beta.flat_size() == n_rows,
         ExcDimensionMismatch(t_size_beta.flat_size(),n_rows));
  Assert(row_id_begin >= 0,ExcLowerRange(row_id_begin,0));
  Assert(row_id_last <= (local_operator.get_num_rows()-1),
         ExcUpperRange(row_id_last,local_operator.get_num_rows()-1));


  const int n_cols = col_id_last - col_id_begin + 1;
  Assert(n_cols >= 1,ExcLowerRange(n_cols,1));
  Assert(n_cols <= local_operator.get_num_cols(),
         ExcUpperRange(n_cols,local_operator.get_num_cols()));
  Assert(t_size_alpha.flat_size() == n_cols,
         ExcDimensionMismatch(t_size_alpha.flat_size(),n_cols));
  Assert(col_id_begin >= 0,ExcLowerRange(col_id_begin,0));
  Assert(col_id_last <= (local_operator.get_num_cols()-1),
         ExcUpperRange(col_id_last,local_operator.get_num_cols()-1));



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
          Index beta_0_1_2 = (beta_2*t_size_beta[1] + beta_1) * t_size_beta[0] + beta_0 + row_id_begin;

          t_id_C2[2] = beta_0_1 ;

          for (Index alpha_2 = 0 ; alpha_2 < t_size_alpha[2] ; ++alpha_2)
          {
            t_id_J[1] = alpha_2;
            for (Index alpha_1 = 0 ; alpha_1 < t_size_alpha[1] ; ++alpha_1)
            {
              for (Index alpha_0 = 0 ; alpha_0 < t_size_alpha[0] ; ++alpha_0)
              {
                Index alpha_0_1 = alpha_1 * t_size_alpha[0] + alpha_0 ;
                Index alpha_0_1_2 = (alpha_2*t_size_alpha[1] + alpha_1) * t_size_alpha[0] + alpha_0 + col_id_begin;

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
          Index beta_0_1_2 = b_tmp_1 + beta_0 + row_id_begin;

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
                const Index alpha_0_1_2 = a_tmp_1 + alpha_0 + col_id_begin;

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


IGA_NAMESPACE_CLOSE

