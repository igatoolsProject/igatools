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


#include <igatools/utils/tensor_size.h>
#include <igatools/base/exceptions.h>

// QualityAssurance: martinelli, 21 Jan 2014

IGA_NAMESPACE_OPEN

template <int rank>
TensorSize<rank>::
TensorSize(Size val) noexcept
{
    Assert(val >= 0, ExcLowerRange(val,0));
    for ( auto & size_dir : size_directions_)
    	size_dir = val;
}



template <int rank>
TensorSize<rank>::
TensorSize(const std::array<Size,rank> &arr) noexcept
:
size_directions_(arr)
{
#ifndef NDEBUG
    for ( const auto & size_dir : size_directions_)
        Assert(size_dir >= 0, ExcLowerRange(size_dir,0));
#endif
}



template <int rank>
TensorSize<rank>::
TensorSize(const TensorIndex<rank> &arr) noexcept
{
	std::copy(arr.begin(), arr.end(), size_directions_.begin());
#ifndef NDEBUG
    for ( const auto & size_dir : size_directions_)
        Assert(size_dir >= 0, ExcLowerRange(size_dir,0));
#endif
}



template <int rank>
inline
TensorSize<rank>::
TensorSize(std::initializer_list<Size> list) noexcept
{
    Assert(list.size() == rank, ExcDimensionMismatch(list.size(),rank));
    std::copy(list.begin(), list.end(), size_directions_.begin());

#ifndef NDEBUG
    for ( const auto & size_dir : size_directions_)
        Assert(size_dir >= 0, ExcLowerRange(size_dir,0));
#endif
}



template <int rank>
Size
TensorSize<rank>::
flat_size() const noexcept
{
    Size res=1;
    for ( const auto & size_dir : size_directions_)
        res *= size_dir;
    return res;
}



template <int rank>
Size &
TensorSize<rank>::
operator()(const Index i)
{
    Assert(i >= 0 && i < rank, ExcIndexRange(i,0,rank));
    return size_directions_[i];
}



template <int rank>
const Size &
TensorSize<rank>::
operator()(const Index i) const
{
    Assert(i >= 0 && i < rank, ExcIndexRange(i,0,rank));
    return size_directions_[i];
}



template <int rank>
bool
operator==(const TensorSize<rank> &size1,const TensorSize<rank> &size2)
{
    bool res = true;
    for (Index i = 0 ; i < rank ; ++i)
        if (size1(i) != size2(i))
        {
            res = false;
            break;
        }

    return res;
}

IGA_NAMESPACE_CLOSE

#include <igatools/utils/tensor_size.inst>
