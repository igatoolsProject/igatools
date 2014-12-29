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


#include <igatools/utils/tensor_sized_container.h>
#include <igatools/base/exceptions.h>
#include <igatools/utils/multi_array_utils.h>

#ifndef NDEBUG
#include <igatools/utils/tensor_sized_container-inline.h>
#endif


IGA_NAMESPACE_OPEN



template <int rank>
TensorSizedContainer<rank>::
TensorSizedContainer()
    :
    TensorSizedContainer(0)
{}

template <int rank>
TensorSizedContainer<rank>::
TensorSizedContainer(const Size size)
    :
    TensorSizedContainer<rank>(TensorSize<rank>(size))
{}

template <int rank>
TensorSizedContainer<rank>::
TensorSizedContainer(const TensorSize<rank> &size)
    :
    size_(size),
    weight_(MultiArrayUtils<rank>::compute_weight(size_))
{}



template <int rank>
TensorSize<rank>
TensorSizedContainer<rank>::
tensor_size() const
{
    return size_;
}


template <int rank>
Size
TensorSizedContainer<rank>::
flat_size() const
{
    return size_.flat_size();
}


template <int rank>
void
TensorSizedContainer<rank>::
reset_size(const TensorSize<rank> &size)
{
    size_ = size;
    weight_ = MultiArrayUtils<rank>::compute_weight(size_) ;
}

template<int rank>
void
TensorSizedContainer<rank>::
print_info(LogStream &out) const
{
    out << "Size: " << size_;
    out << "  Weights: "<< weight_;
}


IGA_NAMESPACE_CLOSE


#include <igatools/utils/tensor_sized_container.inst>
