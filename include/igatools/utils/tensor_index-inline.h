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


// QualityAssurance: martinelli, 21 Jan 2014

#ifndef TENSOR_INDEX_INLINE_H_
#define TENSOR_INDEX_INLINE_H_

#include <igatools/utils/tensor_index.h>

IGA_NAMESPACE_OPEN



template <int rank>
inline
TensorIndex<rank> &
TensorIndex<rank>::
operator +=(const TensorIndex<rank> &ti) noexcept
{
    for (int i = 0; i < rank; ++i)
        (*this)[i] += ti[i];
    return *this;
}

template <int rank>
inline
TensorIndex<rank> &
TensorIndex<rank>::
operator -=(const TensorIndex<rank> &ti) noexcept
{
    for (int i = 0; i < rank; ++i)
    {
        (*this)[i] -= ti[i];
        Assert((*this)[i] >= 0,ExcLowerRange((*this)[i],0));
    }
    return *this;
}

template <int rank>
inline
TensorIndex<rank> &
TensorIndex<rank>::
operator -=(const int j) noexcept
{
    Assert(j >= 0, ExcLowerRange(j,0));
    for (auto &idx : (*this))
    {
        Assert(idx >= j,ExcIndexRange(j,0,idx+1));
        idx -= j;
    }

    return *this;
}

template <int rank>
inline
TensorIndex<rank> &
TensorIndex<rank>::
operator +=(const int j) noexcept
{
    Assert(j >= 0, ExcLowerRange(j,0));
    for (auto &idx : (*this))
        idx += j;

    return *this;
}


template <int rank>
Index &
TensorIndex<rank>::
operator[](const Index i)
{
    Assert(i >= 0 && i < rank, ExcIndexRange(i,0,rank));
    return std::array<Index,rank>::operator[](i);
}



template <int rank>
const Index &
TensorIndex<rank>::
operator[](const Index i) const
{
    Assert(i >= 0 && i < rank, ExcIndexRange(i,0,rank));
    return std::array<Index,rank>::operator[](i);
}



IGA_NAMESPACE_CLOSE


#endif // TENSOR_INDEX_INLINE_H_
