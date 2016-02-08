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

#ifndef INTEGRATOR_SUM_FACTORIZATION_H_
#define INTEGRATOR_SUM_FACTORIZATION_H_

#include <igatools/base/config.h>
#include <igatools/utils/dynamic_multi_array.h>
#include <igatools/utils/multi_array_utils.h>
#include <igatools/linear_algebra/dense_matrix.h>

#include <vector>

IGA_NAMESPACE_OPEN




template <int dim>
class IntegratorSumFactorization
{
public:
  void
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
    DenseMatrix &local_operator) const;
};



template <>
class IntegratorSumFactorization<0>
{
public:
  void operator()(
    const bool is_symmetric,
    const TensorSize<0> &t_size_theta,
    const TensorSize<0> &t_size_alpha,
    const TensorSize<0> &t_size_beta,
    const SafeSTLArray<DynamicMultiArray<Real,3>,0> &J,
    const DynamicMultiArray<Real,3> &C,
    const int row_id_begin,
    const int row_id_last,
    const int col_id_begin,
    const int col_id_last,
    DenseMatrix &local_operator) const;
};


template <>
class IntegratorSumFactorization<1>
{
public:
  void operator()(
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
    DenseMatrix &local_operator) const;
};

template <>
class IntegratorSumFactorization<2>
{
public:
  void operator()(
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
    DenseMatrix &local_operator) const;
};


template <>
class IntegratorSumFactorization<3>
{
public:
  void operator()(
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
    DenseMatrix &local_operator) const;
};

IGA_NAMESPACE_CLOSE



#endif // #ifndef INTEGRATOR_SUM_FACTORIZATION_H_
