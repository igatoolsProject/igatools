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
#include <igatools/utils/cartesian_product_array.h>


IGA_NAMESPACE_OPEN



template <int rank>
CartesianProductIndexer<rank>::
CartesianProductIndexer(const TensorSize<rank> &num_indices_direction)
{
    CartesianProductArray<Index,rank> direction_indices(num_indices_direction) ;

    for (int i = 0 ; i < rank ; ++i)
    {
        Assert(num_indices_direction(i) >= 1,
               ExcLowerRange(num_indices_direction(i),1));

        const Size n_indices_dir = num_indices_direction(i) ;

        for (Index id = 0 ; id < n_indices_dir ; ++id)
            direction_indices.entry(i,id) = id ;
    }

    tensor_indices_ = direction_indices.get_flat_cartesian_product();
}


template <int rank>
TensorIndex<rank>
CartesianProductIndexer<rank>::
get_tensor_index(const Index flat_index) const
{
    Assert((flat_index >=0) && (flat_index < tensor_indices_.size()),
           ExcIndexRange(flat_index,0,tensor_indices_.size()));
    return tensor_indices_[flat_index] ;
}

template <int rank>
Size
CartesianProductIndexer<rank>::
get_num_indices() const
{
    return tensor_indices_.size() ;
}


template <int rank>
void
CartesianProductIndexer<rank>::
print_info_(LogStream &out) const
{
    out << tensor_indices_ << std::endl;
}


IGA_NAMESPACE_CLOSE

#include <igatools/utils/cartesian_product_indexer.inst>

