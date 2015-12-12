//-+--------------------------------------------------------------------
// Igatools a general purpose Isogeometric analysis library.
// Copyright (C) 2012-2015  by the igatools authors (see authors.txt).
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

#ifndef TENSOR_INLINE_H_
#define TENSOR_INLINE_H_

#include <igatools/base/tensor.h>

IGA_NAMESPACE_OPEN

inline
Tdouble::
Tdouble()
  :
  val_(0.0)
{}

inline
Tdouble::
Tdouble(const Real val)
{
  Assert(!std::isnan(val),ExcNotANumber());
  Assert(!std::isinf(val),ExcNumberNotFinite());
  val_ = val;
}



inline
Tdouble::
Tdouble(const bool non_init)
{}



inline
Tdouble &
Tdouble::operator=(const value_t &val)
{
  Assert(!std::isnan(val),ExcNotANumber());
  Assert(!std::isinf(val),ExcNumberNotFinite());

  val_ = val;
  return (*this);
}



inline
auto
Tdouble::operator[](const int i) noexcept -> value_t &
{
  return *this;
}



inline
auto
Tdouble::
operator[](const int i) const noexcept -> const value_t &
{
  return *this;
}



inline
auto
Tdouble::
operator()(const TensorIndex<0>  &i) noexcept -> value_t &
{
  return *this;
}



inline
auto
Tdouble::
operator()(const TensorIndex<0>  &i) const noexcept -> const value_t &
{
  return *this;
}



inline
auto
Tdouble::
operator()(const int i) noexcept -> value_t &
{
  return *this;
}



inline
auto
Tdouble::
operator()(const int i) const noexcept -> const value_t &
{
  return *this;
}



inline
Tdouble &
Tdouble::
operator+=(const Real val) noexcept
{
  Assert(!std::isnan(val),ExcNotANumber());
  Assert(!std::isinf(val),ExcNumberNotFinite());

  val_ += val;
  return *this;
}



inline
Tdouble &
Tdouble::
operator-=(const Real val) noexcept
{
  Assert(!std::isnan(val),ExcNotANumber());
  Assert(!std::isinf(val),ExcNumberNotFinite());

  val_ -= val;
  return *this;
}



inline
Tdouble &
Tdouble::
operator*=(const Real val) noexcept
{
  Assert(!std::isnan(val),ExcNotANumber());
  Assert(!std::isinf(val),ExcNumberNotFinite());

  val_ *= val;
  return *this;
}



inline
Tdouble &
Tdouble::
operator/=(const Real val) noexcept
{
  Assert(!std::isnan(val),ExcNotANumber());
  Assert(!std::isinf(val),ExcNumberNotFinite());
  Assert(val != Real(0.0),ExcDivideByZero());

  val_ /= val;
  return *this;
}



inline
Real
Tdouble::
norm() const noexcept
{
  return std::abs(val_);
}



inline
Real
Tdouble::
norm_square() const noexcept
{
  return val_ * val_;
}



inline
auto
Tdouble::
flat_to_tensor_index(const int flat_index) const noexcept -> TensorIndex<0>
{
  return TensorIndex<0>();
}



inline
int
Tdouble::
tensor_to_flat_index(const TensorIndex<0> &tensor_index) const noexcept
{
  return 0;
}



inline
Size
Tdouble::
get_number_of_entries()
{
  return size;
}



inline
auto
Tdouble::
get_flat_values() const noexcept ->
SafeSTLArray<value_t, n_entries>
{
  return SafeSTLArray<value_t, n_entries> ({val_});
}




/*------ Inline functions: Tensor<  dim_, rank_, tensor_type, value_type > ---*/

template<int dim_, int rank_, class tensor_type, class value_type >
inline
Tensor< dim_, rank_, tensor_type, value_type >::
Tensor(std::initializer_list<value_type> list)
{
  Assert(list.size() == self_t::size,
         ExcDimensionMismatch(list.size(),self_t::size)) ;
  for (int i = 0; i < self_t::size; ++i)
    tensor_[i] = list.begin()[i];
}



template<int dim_, int rank_, class tensor_type, class value_type >
inline
Tensor<dim_, rank_, tensor_type, value_type> &
Tensor< dim_, rank_, tensor_type, value_type >::
operator=(std::initializer_list<value_type> list)
{
  Assert(list.size() == self_t::size,
         ExcDimensionMismatch(list.size(),self_t::size));
  for (int i = 0; i < self_t::size; ++i)
    tensor_[i] = list.begin()[i];
  return *this;
}



template<int dim_, int rank_, class tensor_type, class value_type >
inline
Tensor<dim_, rank_, tensor_type, value_type> &
Tensor< dim_, rank_, tensor_type, value_type >::
operator=(const value_type &entry_val)
{
  Assert(dim_ == 1 && rank_ == 1,
         ExcMessage("Assignment is allowed only if dim_==1 and rank_==1")) ;
  for (auto &val : tensor_)
    val = entry_val;
  return *this;
}



template<int dim_, int rank_, class tensor_type, class value_type >
inline
auto
Tensor< dim_, rank_, tensor_type, value_type >::
operator[](const int  i) -> SubTensor<self_t> &
{
  Assert(i >= 0 && i < dim_, ExcIndexRange(i, 0, dim_)) ;
  return tensor_[i] ;
}



template<int dim_, int rank_, class tensor_type, class value_type >
inline
auto
Tensor< dim_, rank_, tensor_type, value_type >::
operator[](const int  i) const -> const SubTensor<self_t> &
{
  Assert(i >= 0 && i < dim_, ExcIndexRange(i, 0, dim_)) ;
  return tensor_[i] ;
}



template<int dim_, int rank_, class tensor_type, class value_type >
inline
value_type &
Tensor< dim_, rank_, tensor_type, value_type >::
operator()(const TensorIndex<rank_> &i)
{
  Assert(tensor_to_flat_index(i) < self_t::size,
         ExcIndexRange(tensor_to_flat_index(i), 0, self_t::size)) ;
  return (reinterpret_cast<value_type *>(tensor_)[tensor_to_flat_index(i)]);
}



template<int dim_, int rank_, class tensor_type, class value_type >
inline
const value_type &
Tensor< dim_, rank_, tensor_type, value_type >::
operator()(const TensorIndex<rank_> &i) const
{
  Assert(tensor_to_flat_index(i) < self_t::size,
         ExcIndexRange(tensor_to_flat_index(i), 0, self_t::size)) ;
  return (reinterpret_cast<value_type const *>(tensor_)[tensor_to_flat_index(i)]);
}



template<int dim_, int rank_, class tensor_type, class value_type >
inline
value_type &
Tensor< dim_, rank_, tensor_type, value_type >::
operator()(const int i)
{
  Assert(i >= 0 && i < self_t::size, ExcIndexRange(i, 0, self_t::size)) ;
  return (reinterpret_cast<value_type *>(tensor_)[i]);
}



template<int dim_, int rank_, class tensor_type, class value_type >
inline
const value_type &
Tensor< dim_, rank_, tensor_type, value_type >::
operator()(const int i) const
{
  Assert(i >= 0 && i < self_t::size, ExcIndexRange(i, 0, self_t::size)) ;
  return (reinterpret_cast<value_type const *>(tensor_)[i]);
}



template<int dim_, int rank_, class tensor_type, class value_type >
inline
Tensor< dim_, rank_, tensor_type, value_type > &
Tensor< dim_, rank_, tensor_type, value_type >::
operator=(const Real value)
{
  Assert(value==0.0 || (dim_ == 1 && rank_ == 1),
         ExcMessage("Assignment with non-zero value is allowed only if dim_==1 and rank_==1")) ;

  Assert(!std::isnan(value),ExcNotANumber());
  Assert(!std::isinf(value),ExcNumberNotFinite());

  for (auto &tensor_component : tensor_)
    tensor_component = value ;

  return *this;
}



template<int dim_, int rank_, class tensor_type, class value_type >
inline
Tensor< dim_, rank_, tensor_type, value_type > &
Tensor< dim_, rank_, tensor_type, value_type >::
operator+=(const Tensor< dim_, rank_, tensor_type, value_type > &tensor)
{
  for (int i = 0 ; i < dim_ ; i++)
    tensor_[i] += tensor[i] ;

  return *this  ;
}



template<int dim_, int rank_, class tensor_type, class value_type >
inline
Tensor< dim_, rank_, tensor_type, value_type > &
Tensor< dim_, rank_, tensor_type, value_type >::
operator-=(const Tensor< dim_, rank_, tensor_type, value_type > &tensor)
{
  for (int i = 0 ; i < dim_ ; i++)
    tensor_[i] -= tensor.tensor_[i] ;

  return *this ;
}



template<int dim_, int rank_, class tensor_type, class value_type >
inline
Tensor< dim_, rank_, tensor_type, value_type > &
Tensor< dim_, rank_, tensor_type, value_type >::
operator*=(const Real value)
{
  Assert(!std::isnan(value),ExcNotANumber());
  Assert(!std::isinf(value),ExcNumberNotFinite());

  for (auto &tensor_component : tensor_)
    tensor_component *= value ;

  return *this ;
}



template<int dim_, int rank_, class tensor_type, class value_type >
inline
Tensor< dim_, rank_, tensor_type, value_type > &
Tensor< dim_, rank_, tensor_type, value_type >::
operator/=(const Real value)
{
  Assert(!std::isnan(value),ExcNotANumber());
  Assert(!std::isinf(value),ExcNumberNotFinite());
  Assert(value != Real(0.0),ExcDivideByZero());

  for (auto &tensor_component : tensor_)
    tensor_component /= value ;

  return *this ;
}


template<int dim_, int rank_, class tensor_type, class value_type >
inline
Tensor< dim_, rank_, tensor_type, value_type > &
Tensor< dim_, rank_, tensor_type, value_type >::
operator-()
{
  for (auto &tensor_component : tensor_)
    tensor_component = -tensor_component ;

  return *this ;
}



template<int dim_, int rank_, class tensor_type, class value_type >
inline
Real
Tensor< dim_, rank_, tensor_type, value_type >::
norm() const noexcept
{
  return std::sqrt(norm_square());
}



template<int dim_, int rank_, class tensor_type, class value_type >
inline
Real
Tensor< dim_, rank_, tensor_type, value_type >::
norm_square() const noexcept
{
  Real s = 0.;
  for (int i=0; i<dim_; ++i)
    s += tensor_[i].norm_square();

  return s;
}



template<int dim_, int rank_, class tensor_type, class value_type >
inline
auto
Tensor< dim_, rank_, tensor_type, value_type >::
flat_to_tensor_index(const int flat_index) const noexcept -> TensorIndex<rank_>
{
  Assert(flat_index < self_t::size,
  ExcIndexRange(flat_index, 0, self_t::size)) ;

  TensorIndex<rank_> tensor_index;

  int l = flat_index;
  for (int i = 0; i < rank_ ; ++i)
  {
    const int w = pow(dim_,rank_-i-1);
    tensor_index[i] = l / w;
    l %= w;
  }

  return tensor_index;
}



template<int dim_, int rank_, class tensor_type, class value_type >
inline
Size
Tensor< dim_, rank_, tensor_type, value_type >::
get_number_of_entries()
{
  return n_entries;
}



template<int dim_, int rank_, class tensor_type, class value_type >
inline
int
Tensor< dim_, rank_, tensor_type, value_type >::
tensor_to_flat_index(const TensorIndex<rank_> &tensor_index) const noexcept
{
  int flat_index = 0;
  for (int i = 0; i < rank_; ++i)
    flat_index += pow(dim_,i) * tensor_index[rank_-1-i];

  return flat_index;
}



template<int dim_, int rank_, class tensor_type, class value_type >
inline
auto
Tensor< dim_, rank_, tensor_type, value_type >::
get_flat_values() const noexcept ->
SafeSTLArray<typename Tdouble::value_t, n_entries>
{
  SafeSTLArray<typename Tdouble::value_t, n_entries> values;
  for (int i = 0; i < dim_; ++i)
  {
    const auto new_values= tensor_[i].get_flat_values();
    std::copy (new_values.cbegin(), new_values.cend(),
               values.begin() + i * value_type::n_entries);
  }

  return values;
}


//--------------  Related functions   ------------------------------------------
template<class T>
inline
EnableIf<T::is_tensor, T>
operator+(const T &A, const T &B) noexcept
{
  T R(A);
  R += B;
  return R;
}



template<class T>
inline
EnableIf<T::is_tensor,T>
operator-(const T &A, const T &B) noexcept
{
  T R(A);
  R -= B;
  return R;
}



template<class T>
inline
EnableIf<T::is_tensor,T>
operator*(const T &A, const Real scalar) noexcept
{
  Assert(!std::isnan(scalar),ExcNotANumber());
  Assert(!std::isinf(scalar),ExcNumberNotFinite());

  T R(A);
  R *= scalar;
  return (R);
}



template<class T>
inline
EnableIf<T::is_tensor,T>
operator*(const Real scalar, const T &A) noexcept
{
  return A * scalar;
}



template<class T>
inline
EnableIf<T::is_tensor,T>
operator/(const T &A, const Real scalar) noexcept
{
  Assert(!std::isnan(scalar), ExcNotANumber());
  Assert(!std::isinf(scalar), ExcNumberNotFinite());
  Assert(scalar !=0 , ExcDivideByZero());

  T R(A);
  R /= scalar;
  return R;
}

IGA_NAMESPACE_CLOSE

#endif /* TENSOR_INLINE_H_ */
