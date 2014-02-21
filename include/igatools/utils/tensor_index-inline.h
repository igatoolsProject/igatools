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
TensorIndex<rank>::
TensorIndex(Index val) noexcept
{
    Assert(val >= 0, ExcLowerRange(val,0));
    for (int i= 0 ; i < rank ; ++i)
        (*this)[i] = val;
}


template <int rank>
inline
TensorIndex<rank>::
TensorIndex(const std::array<int,rank> &arr) noexcept
:
std::array<int,rank>::array(arr)
{
#ifndef NDEBUG
    for (int i = 0; i < rank; ++i)
        Assert((*this)[i] >= 0,ExcLowerRange((*this)[i],0));
#endif
}

template <int rank>
inline
TensorIndex<rank>::
TensorIndex(std::initializer_list<Index> list) noexcept
{
    Assert(list.size() == rank, ExcDimensionMismatch(list.size(),rank));
    std::copy(list.begin(), list.end(), this->begin());

#ifndef NDEBUG
    for (int i = 0; i < rank; ++i)
        Assert((*this)[i] >= 0,ExcLowerRange((*this)[i],0));
#endif
}


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
operator -=(const int j) noexcept
{
    for (int i = 0; i < rank; ++i)
    {
        Assert(j >= 0, ExcLowerRange(j,0));
        Assert((*this)[i] >= j,ExcIndexRange(j,0,(*this)[i]+1));
        (*this)[i] -= j;
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
    for (int i = 0; i < rank; ++i)
    {
        Assert(j >= 0, ExcLowerRange(j,0));
        (*this)[i] += j;
    }
    return *this;
}



IGA_NAMESPACE_CLOSE


#endif // TENSOR_INDEX_INLINE_H_
