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


#include <chrono>


//#define TIME_PROFILING

IGA_NAMESPACE_OPEN



#ifdef TIME_PROFILING
using Clock = std::chrono::high_resolution_clock;
using TimePoint = std::chrono::time_point<Clock>;
using Duration = std::chrono::duration<Real>;
#endif // TIME_PROFILING




/*

template <int dim>
void
IntegratorSumFactorization<dim>::
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
//*/





void
IntegratorSumFactorization<0>::
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
IntegratorSumFactorization<1>::
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
  const int n_rows = row_id_last - row_id_begin + 1;
#ifndef NDEBUG
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
#endif

  TensorIndex<3> t_id_J;
  TensorIndex<3> t_id_C;

  if (!is_symmetric)
  {
    //--------------------------------------------------------------
    t_id_C[1] = 0;
    t_id_C[2] = 0;

    for (Index beta_0 = 0 ; beta_0 < t_size_beta[0] ; ++beta_0)
    {
      t_id_J[2] = beta_0 + row_id_begin;
      for (Index alpha_0 = 0 ; alpha_0 < t_size_alpha[0] ; ++alpha_0)
      {
        t_id_J[1] = alpha_0 + col_id_begin;

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
      t_id_J[2] = beta_0 + row_id_begin;
      for (Index alpha_0 = beta_0 ; alpha_0 < t_size_alpha[0] ; ++alpha_0)
      {
        t_id_J[1] = alpha_0 + col_id_begin;

        Real sum = 0.0;
        for (Index theta_0 = 0; theta_0 < t_size_theta[0] ; ++theta_0)
        {
          t_id_J[0] = theta_0;
          t_id_C[0] = theta_0;
          sum += C(t_id_C) * J[0](t_id_J);
        }

        local_operator(beta_0,alpha_0) += sum;
      }
    }

    // here we copy the upper triangular part of the current block on the lower triangular part
    for (int loc_row = 0 ; loc_row < n_rows ; ++loc_row)
    {
      const int r_src = loc_row + row_id_begin;
      const int c_tgt = loc_row + col_id_begin;
      for (int loc_col = loc_row+1 ; loc_col < n_rows ; ++loc_col)
      {
        const int r_tgt = loc_col + row_id_begin;
        const int c_src = loc_col + col_id_begin;
        local_operator(r_tgt,c_tgt) =
          local_operator(r_src,c_src);
      } // end loop loc_col
    } // end_loop loc_row
    //--------------------------------------------------------------

  } // end if (is_symmetric)
}

void
IntegratorSumFactorization<2>::
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
#ifndef NDEBUG
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
#endif





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
  {
    t_id_C[0] = 0;
    t_id_C[1] = 0;
    t_id_C[2] = 0;
    const Real *C_it_begin = &C(t_id_C);

    t_id_J[0] = 0;

    for (Index beta_0 = 0 ; beta_0 < t_size_beta[0] ; ++beta_0)
    {
      t_id_J [2] = beta_0;
      t_id_C1[2] = beta_0;
      for (Index alpha_0 = 0 ; alpha_0 < t_size_alpha[0] ; ++alpha_0)
      {
        t_id_J [1] = alpha_0;
        t_id_C1[1] = alpha_0;

        const Real *const J0_it_begin = &J[0](t_id_J);
        const Real *const J0_it_end = J0_it_begin + t_size_theta[0];

        t_id_C1[0] = 0;
        Real *const C1_it_begin = &C1(t_id_C1);
        Real *const C1_it_end = C1_it_begin + t_size_theta[1];

        const Real *C_it = C_it_begin;
        for (Real *C1_it = C1_it_begin; C1_it != C1_it_end ; ++C1_it, C_it += t_size_theta[0])
        {
          (*C1_it) = std::inner_product(J0_it_begin,J0_it_end,C_it,0.0);
        } //end loop theta_1

      } //end loop alpha_0
    } // end loop beta_0
  }
  //--------------------------------------------------------------


  //--------------------------------------------------------------
  t_id_C1[0] = 0;
  t_id_J[0] = 0;

  for (Index beta_1 = 0 ; beta_1 < t_size_beta[1] ; ++beta_1)
  {
    t_id_J[2] = beta_1;
    for (Index beta_0 = 0 ; beta_0 < t_size_beta[0] ; ++beta_0)
    {
      Index beta_0_1 = beta_1*t_size_beta[0] + beta_0 + row_id_begin;

      t_id_C1[2] = beta_0;

      for (Index alpha_1 = (is_symmetric ? beta_1 : 0) ;
           alpha_1 < t_size_alpha[1] ;
           ++alpha_1)
      {
        t_id_J[1] = alpha_1;
        t_id_C1[1] = 0;

        const Real *J1_it_begin = &J[1](t_id_J);

        Real *row_it = &local_operator(beta_0_1,alpha_1*t_size_alpha[0] + col_id_begin);
        const Real *const row_it_end = row_it + t_size_alpha[0];
        for (; row_it != row_it_end ; ++t_id_C1[1], ++row_it)
        {
          const Real *C1_it_begin = &C1(t_id_C1);
          const Real *const C1_it_end = C1_it_begin + t_size_theta[1];


          (*row_it) += std::inner_product(C1_it_begin,C1_it_end,J1_it_begin,0.0);

        } //end loop row_it
      } //end loop alpha_1
    } // end loop beta_0
  } // end loop beta_1

  if (is_symmetric)
  {
    // here we copy the upper triangular part of the current block on the lower triangular part
    for (int loc_row = 0 ; loc_row < n_rows ; ++loc_row)
    {
      const int r_src = loc_row + row_id_begin;
      const int c_tgt = loc_row + col_id_begin;
      for (int loc_col = loc_row+1 ; loc_col < n_rows ; ++loc_col)
      {
        const int r_tgt = loc_col + row_id_begin;
        const int c_src = loc_col + col_id_begin;
        local_operator(r_tgt,c_tgt) =
          local_operator(r_src,c_src);
      } // end loop loc_col
    } // end_loop loc_row
  }//end if (is_symmetric)
}


void
IntegratorSumFactorization<3>::
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
#ifndef NDEBUG
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
#endif



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
#ifdef TIME_PROFILING
  const auto start_l1 = Clock::now();
#endif // TIME_PROFILING

  t_id_J[0] = 0;

  t_id_C[0] = 0;
  t_id_C[1] = 0;
  t_id_C[2] = 0;

  t_id_C1[0] = 0;
  for (Index beta_0 = 0 ; beta_0 < t_size_beta[0] ; ++beta_0)
  {
    t_id_J [2] = beta_0;
    t_id_C1[2] = beta_0;
    for (Index alpha_0 = 0 ; alpha_0 < t_size_alpha[0] ; ++alpha_0)
    {
      t_id_J [1] = alpha_0;
      t_id_C1[1] = alpha_0;


      Real *C1_it = &C1(t_id_C1);

      const Real *C_it_begin = &C(t_id_C);

      const Real *const J0_it_begin = &J[0](t_id_J);
      const Real *const J0_it_end = J0_it_begin + t_size_theta[0];

      for (Index theta_2 = 0; theta_2 < t_size_theta[2] ; ++theta_2)
      {
        const Real *const C1_it_end = C1_it + t_size_theta[1];
        for (; C1_it != C1_it_end ; ++C1_it, C_it_begin += t_size_theta[0])
        {
          (*C1_it) = std::inner_product(J0_it_begin,J0_it_end,C_it_begin,0.0);
        } // end loop theta_1
      } // end loop theta_2
    } // end loop alpha_0
  } // end loop beta_0
#ifdef TIME_PROFILING
  const Duration time_l1 = Clock::now() - start_l1;
  std::cout << "Elapsed_seconds loop1 = " << time_l1.count() << std::endl;
#endif // TIME_PROFILING

  //--------------------------------------------------------------



  //--------------------------------------------------------------
#ifdef TIME_PROFILING
  const auto start_l2 = Clock::now();
#endif // TIME_PROFILING

  t_id_C2[0] = 0;
  t_id_C1[0] = 0;
  t_id_J[0] = 0;

  for (t_id_J[2] = 0 ; t_id_J[2] < t_size_beta[1] ; ++t_id_J[2])
  {
    t_id_C2[2] = t_id_J[2]*t_size_beta[0];
    for (t_id_C1[2] = 0 ; t_id_C1[2] < t_size_beta[0] ; ++t_id_C1[2], ++t_id_C2[2])
    {
      t_id_C2[1] = 0;
      for (t_id_J[1] = 0 ;
           t_id_J[1] < t_size_alpha[1] ;
           ++t_id_J[1],
           t_id_C2[1] += t_size_alpha[0])
      {
        t_id_C1[1] = 0;
        t_id_C2[1] = t_id_J[1]*t_size_alpha[0];

        const Real *const J1_it_begin = &J[1](t_id_J);

        for (t_id_C1[1] = 0 ; t_id_C1[1] < t_size_alpha[0] ; ++t_id_C1[1],++t_id_C2[1])
        {
          const Real *const C1_it_begin = &C1(t_id_C1);
          const Real *const C1_it_end = C1_it_begin + t_size_theta[1];


          Real *const C2_it_begin = &C2(t_id_C2);
          const Real *const C2_it_end = C2_it_begin + t_size_C2[0];
          for (Real *C2_it = C2_it_begin; C2_it != C2_it_end ; ++C2_it)
          {
            (*C2_it) = std::inner_product(C1_it_begin,C1_it_end,J1_it_begin,0.0);
          } // end loop C2_it
        } //end loop t_id_C1[1]
      } //end loop t_id_J[1]
    } // end loop beta_0
  } // end loop beta_1
#ifdef TIME_PROFILING
  const Duration time_l2 = Clock::now() - start_l2;
  std::cout << "Elapsed_seconds loop2 = " << time_l2.count() << std::endl;
#endif // TIME_PROFILING

  //--------------------------------------------------------------


  //--------------------------------------------------------------
#ifdef TIME_PROFILING
  const auto start_l3 = Clock::now();
#endif // TIME_PROFILING


  t_id_C2[0] = 0;
  t_id_C2[1] = 0;
  t_id_J[0] = 0;

  const int n_pts = t_size_C2[0];

  const int size_alpha_0 = t_size_alpha[0];
  const int size_alpha_1 = t_size_alpha[1];

//  bool not_symm =false;

  for (t_id_J[2] = 0 ; t_id_J[2] < t_size_beta[2] ; ++t_id_J[2])
  {
    const Index b_tmp_2 = t_id_J[2]*t_size_beta[1];

    t_id_C2[2] = 0;
    for (Index beta_1 = 0 ; beta_1 < t_size_beta[1] ; ++beta_1)
    {
      const Index beta_0_1_2 = (b_tmp_2 + beta_1) * t_size_beta[0] + row_id_begin;

      int r_id = beta_0_1_2;

      for (Index beta_0 = 0;
           beta_0 < t_size_beta[0] ;
           ++beta_0,
           ++t_id_C2[2],
           ++r_id)
      {
        Real *const row_it_begin_0 = &local_operator(r_id,0);

        for (t_id_J[1] = (is_symmetric ? t_id_J[2] : 0) ;
             t_id_J[1] < t_size_alpha[2] ;
             ++t_id_J[1])
        {
          const Real *const J2_it_begin = &J[2](t_id_J);

          const Real *C2_it_begin = &C2(t_id_C2);
          const Real *C2_it_end = C2_it_begin + n_pts;

          const int c_begin = t_id_J[1] * size_alpha_1 * size_alpha_0 + col_id_begin;

          Real *row_it_begin = row_it_begin_0 + c_begin;
          const Real *row_it_end = row_it_begin + size_alpha_0;

          const Real *const block_end = row_it_begin + size_alpha_0 * size_alpha_1;

          for (;
               row_it_begin != block_end ;
               row_it_begin += size_alpha_0,
               row_it_end   += size_alpha_0) // loop of size_alpha_1 iterations
          {
            for (Real *row_it = row_it_begin ;
                 row_it != row_it_end ;
                 C2_it_begin += n_pts,
                 C2_it_end += n_pts,
                 ++row_it) // loop of size_alpha_0 iterations
            {
              (*row_it) += std::inner_product(C2_it_begin,C2_it_end,J2_it_begin,0.0);
            } // end loop row_it
          } // end loop alpha_1
        } // end loop alpha_2
      } // end loop beta_0
    } // end loop beta_1
  } // end loop beta_2
#ifdef TIME_PROFILING
  const Duration time_l3 = Clock::now() - start_l3;
  std::cout << "Elapsed_seconds loop3 = " << time_l3.count() << std::endl;
#endif // TIME_PROFILING

  //--------------------------------------------------------------


  if (is_symmetric)
  {
#ifdef TIME_PROFILING
    const auto start_l4 = Clock::now();
#endif // TIME_PROFILING

    // here we copy the upper triangular part of the current block on the lower triangular part
    for (int loc_row = 0 ; loc_row < n_rows ; ++loc_row)
    {
      const int r_src = loc_row + row_id_begin;
      const int c_tgt = loc_row + col_id_begin;
      for (int loc_col = loc_row+1 ; loc_col < n_rows ; ++loc_col)
      {
        const int r_tgt = loc_col + row_id_begin;
        const int c_src = loc_col + col_id_begin;
        local_operator(r_tgt,c_tgt) =
          local_operator(r_src,c_src);
      } // end loop loc_col
    } // end_loop loc_row
    //*/

    //*/
    //--------------------------------------------------------------
#ifdef TIME_PROFILING
    const Duration time_l4 = Clock::now() - start_l4;
    std::cout << "Elapsed_seconds loop4 = " << time_l4.count() << std::endl;
#endif // TIME_PROFILING


  } // end if (is_symmetric)
}


IGA_NAMESPACE_CLOSE

