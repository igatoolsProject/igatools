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

#ifndef MULTI_ARRAY_UTILS_INLINE_H_
#define MULTI_ARRAY_UTILS_INLINE_H_

#include <igatools/utils/multi_array_utils.h>

IGA_NAMESPACE_OPEN


template <int rank>
inline
TensorIndex<rank>
MultiArrayUtils<rank>::
flat_to_tensor_index(const Index flat_index,
                     const TensorIndex<rank> &weight) noexcept
{
    /*
     * To compute the tensor index we basically use
     * the integer division algorithm.
     * The last index is how many times weight[last index] units
     * fits in flat_index and continue this weight with the
     * remainder.
     */
    TensorIndex<rank> tensor_index;

    int l = flat_index;
    for (int i = rank-1; i > 0 ; --i)
    {
        tensor_index[i] = l / weight[i-1];
        l %= weight[i-1];
    }
    tensor_index[0] = l;

    return tensor_index;
}

template <>
inline
TensorIndex<0>
MultiArrayUtils<0>::
flat_to_tensor_index(const Index flat_index,
                     const TensorIndex<0> &weight) noexcept
{
    TensorIndex<0> tensor_index;
    return tensor_index;
}

template <>
inline
Index
MultiArrayUtils<0>::
tensor_to_flat_index(const TensorIndex<0> &tensor_index,
                     const TensorIndex<0> &weight) noexcept
{
    return 0;
}

template <int rank>
inline
Index
MultiArrayUtils<rank>::
tensor_to_flat_index(const TensorIndex<rank> &tensor_index,
                     const TensorIndex<rank> &weight) noexcept
{
    Index flat_index = tensor_index[0];
    for (int i = 1; i < rank; ++i)
        flat_index += weight[i-1] * tensor_index[i];

    return flat_index;
}



template <int rank>
inline
TensorIndex<rank>
MultiArrayUtils<rank>::
compute_weight(const TensorSize<rank> &size) noexcept
{
    TensorIndex<rank> weight;

    weight[0] = size[0];
    for (int i = 1; i < rank; ++i)
        weight[i] = weight[i-1] * size[i];

    return weight;
}

template <>
inline
TensorIndex<0>
MultiArrayUtils<0>::
compute_weight(const TensorSize<0> &size) noexcept
{
    TensorIndex<0> weight;
    return weight;
}


template <int rank>
inline
Size
MultiArrayUtils<rank>::
size(const TensorIndex<rank> &extend) noexcept
{
    Size res=1;
    for (int i = 0; i < rank; ++i)
    {
        Assert(extend[i] > 0, ExcLowerRange(extend[i],1)) ;
        res *= extend[i];
    }
    return res;
}

template <int rank>
inline
TensorIndex<rank>
MultiArrayUtils<rank>::
get_tensor_entry(const std::array< vector<Index>, rank> &data,
                 const Index flat_index, const TensorIndex<rank> &weight) noexcept
{
    TensorIndex<rank> entry;
    auto tensor_index = flat_to_tensor_index(flat_index,weight);
    for (int i = 0; i < rank; ++i)
        entry[i]=data[i][tensor_index[i]];

    return entry;
}


IGA_NAMESPACE_CLOSE


#endif // MULTI_ARRAY_UTILS_INLINE_H_
