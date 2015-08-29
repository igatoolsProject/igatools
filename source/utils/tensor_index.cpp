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


#include <igatools/utils/tensor_index.h>
#include <igatools/base/exceptions.h>

// TODO (pauletti, Apr 22, 2015): the inlining in release mode should be uniform
// through out the library
#ifndef NDEBUG
#include <igatools/utils/tensor_index-inline.h>
#endif

IGA_NAMESPACE_OPEN

template <int rank>
TensorIndex<rank>::
TensorIndex(const Size val) noexcept
{
  Assert(val >= 0, ExcLowerRange(val,0));
  for (auto &idx : (*this))
    idx = val;
}



template <int rank>
TensorIndex<rank>::
TensorIndex(const SafeSTLArray<Index,rank> &arr) noexcept
  :
  SafeSTLArray<Index,rank>(arr)
{
#ifndef NDEBUG
  for (const auto &idx : (*this))
    Assert(idx >= 0,ExcLowerRange(idx,0));
#endif
}



template <int rank>
TensorIndex<rank>::
TensorIndex(std::initializer_list<Index> list) noexcept
{
  if (rank > 0)
  {
    Assert(list.size() == rank, ExcDimensionMismatch(list.size(),rank));
    std::copy(list.begin(), list.end(), this->begin());
  }

#ifndef NDEBUG
  for (const auto &idx : (*this))
    Assert(idx >= 0,ExcLowerRange(idx,0));
#endif
}



template <int rank>
std::size_t
TensorIndex<rank>::
memory_consumption() const
{
  return sizeof(static_cast<const SafeSTLArray<Index,rank> &>(*this));
}



template <int rank>
TensorIndex<rank>
operator+(const TensorIndex<rank> &index_a,const TensorIndex<rank> &index_b)
{
  TensorIndex<rank> tensor_index;
  for (int i = 0 ; i < rank ; ++i)
    tensor_index[i] = index_a[i] + index_b[i];

  return tensor_index;
}



template <int rank>
TensorIndex<rank>
operator+(const TensorIndex<rank> &index,const Index j)
{
  Assert(j>=0,ExcLowerRange(j,0));
  TensorIndex<rank> tensor_index;
  for (int i = 0 ; i < rank ; ++i)
    tensor_index[i] = index[i] + j;

  return tensor_index;
}



template <int rank>
TensorIndex<rank>
operator-(const TensorIndex<rank> &index,const Index j)
{
  Assert(j>=0,ExcLowerRange(j,0));
  TensorIndex<rank> tensor_index;
  for (int i = 0 ; i < rank ; ++i)
  {
    tensor_index[i] = index[i] - j;
    Assert(tensor_index[i] >= 0, ExcLowerRange(tensor_index[i],0));
  }

  return tensor_index;
}



template <int rank>
LogStream &
operator<<(LogStream &out, const TensorIndex<rank> &tensor_index)
{
  out << "[";
  if (rank > 0)
  {
    out << tensor_index[0];
    for (int i = 1 ; i < rank ; ++i)
      out << "," << tensor_index[i];
  }
  out << "]";

  return out;
}



#if 0
template<>
SafeSTLVector<TensorIndex<1> >
tensor_range(TensorIndex<1> first, TensorIndex<1> last)
{
  Assert(first <= last, ExcMessage("first bigger than last"));
  SafeSTLVector<TensorIndex<1>> result(last[0]-first[0]);
  for (int i=first[0]; i<last[0]; ++i)
    result[i-first[0]][0] = i;

  return result;
}



template<>
SafeSTLVector<TensorIndex<0> >
tensor_range(TensorIndex<0> first, TensorIndex<0> last)
{
  Assert(false, ExcNotImplemented());
  SafeSTLVector<TensorIndex<0>> result;

  return result;
}
#endif


IGA_NAMESPACE_CLOSE

#include <igatools/utils/tensor_index.inst>
