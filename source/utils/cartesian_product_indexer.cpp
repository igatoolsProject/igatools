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

// QualityAssurance: martinelli, 29 Jan 2014


#include <igatools/utils/cartesian_product_indexer.h>
#include <igatools/base/exceptions.h>


IGA_NAMESPACE_OPEN



template <int rank>
CartesianProductIndexer<rank>::
CartesianProductIndexer(const TensorSize<rank> &num_indices_direction)
    :
    DynamicMultiArray<TensorIndex<rank>,rank>(num_indices_direction)
{
    const Size flat_size = this->flat_size();

    for (Index flat_id = 0 ; flat_id < flat_size ; ++flat_id)
        (*this)(flat_id) = this->flat_to_tensor(flat_id);
}


template <int rank>
TensorIndex<rank>
CartesianProductIndexer<rank>::
get_tensor_index(const Index flat_index) const
{
    Assert((flat_index >=0) && (flat_index < this->flat_size()),
           ExcIndexRange(flat_index,0,this->flat_size()));
    return this->data_[flat_index] ;
}


template <int rank>
Size
CartesianProductIndexer<rank>::
get_num_indices() const
{
    return this->flat_size() ;
}



IGA_NAMESPACE_CLOSE

#include <igatools/utils/cartesian_product_indexer.inst>

