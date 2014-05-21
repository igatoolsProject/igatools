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


#ifndef CARTESIAN_PRODUCT_ARRAY_INLINE_H_
#define CARTESIAN_PRODUCT_ARRAY_INLINE_H_


#include<igatools/utils/cartesian_product_array.h>

IGA_NAMESPACE_OPEN

template< class T, int rank>
inline
const std::vector<T> &
CartesianProductArray<T,rank>::
get_data_direction(const int i) const
{
    Assert(i >= 0 && i < rank, ExcIndexRange(i, 0, rank)) ;
    return data_[i];
}


IGA_NAMESPACE_CLOSE
#endif // #ifndef CARTESIAN_PRODUCT_ARRAY_INLINE_H_

