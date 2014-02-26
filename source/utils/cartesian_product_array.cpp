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

#include <igatools/utils/cartesian_product_array.h>
#include <igatools/base/exceptions.h>


using std::array;
using std::vector;

IGA_NAMESPACE_OPEN

template< class T, int rank>
CartesianProductArray<T,rank>::
CartesianProductArray()
    :
    CartesianProductArray(0)
{}


template< class T, int rank>
CartesianProductArray<T,rank>::
CartesianProductArray(std::initializer_list<std::initializer_list<T>> list)
{
    Assert(list.size() == rank,
           ExcDimensionMismatch(list.size(),rank)) ;

    TensorSize<rank> size;
    for (int i=0; i<rank; ++i)
    {
        data_[i] = list.begin()[i];
        size(i) = data_[i].size();
    }
    TensorSizedContainer<rank>::reset_size(size);
}


template< class T, int rank>
CartesianProductArray<T,rank>::
CartesianProductArray(const TensorSize<rank> size)
    :
    TensorSizedContainer<rank>(size)
{
    for (int i=0; i<rank; ++i)
    {
        Assert(size(i) >= 0, ExcLowerRange(size(i),0)) ;
        data_[i].resize(size(i));
    }
}


template< class T, int rank>
CartesianProductArray<T,rank>::
CartesianProductArray(const Size size)
    :
    CartesianProductArray(TensorSize<rank>(filled_array<Size,rank>(size)))
{}


template< class T, int rank>
CartesianProductArray<T,rank>::
CartesianProductArray(const array<vector<T>,rank> &data_directions)
{
    TensorSize<rank> size;
    for (int i=0; i<rank; ++i)
    {
        data_[i] = data_directions[i];
        size(i) = data_[i].size();
    }
    TensorSizedContainer<rank>::reset_size(size);
}


template< class T, int rank>
void CartesianProductArray<T,rank>::
resize(const TensorSize<rank> &size)
{
    TensorSizedContainer<rank>::reset_size(size);

    for (int i = 0 ; i < rank ; i++)
    {
        Assert(size(i) >= 1, ExcLowerRange(size(i), 1)) ;
        data_[i].resize(size(i)) ;
    }
}


template< class T, int rank>
T &
CartesianProductArray<T,rank>::
entry(const int i, const int j)
{
    Assert(i >= 0 && i < rank, ExcIndexRange(i,0,rank));
    Assert(j >= 0 && j < this->tensor_size()(i),
           ExcIndexRange(j,0,this->tensor_size()(i)));

    return data_[i][j];
}

template< class T, int rank>
void
CartesianProductArray<T,rank>::
copy_data_direction(const int i, const vector<T> &data)
{
    Assert(i>=0 && i<rank, ExcIndexRange(i,0,rank));
    Assert(data.size()>0, ExcLowerRange(data.size(),0));
    data_[i] = data;
    TensorSize<rank> size = this->tensor_size();
    if (data_[i].size() != size(i))
    {
        size(i) = data_[i].size();
        TensorSizedContainer<rank>::reset_size(size);
    }
}






template< class T, int rank>
auto
CartesianProductArray<T,rank>::
cartesian_product(const TensorIndex<rank> &index) const -> point_t
{
    point_t result;
    for (int i = 0; i < rank; ++i)
        result[i] = data_[i][index[i]];
    return result;
}



template< class T, int rank>
auto
CartesianProductArray<T,rank>::
get_flat_cartesian_product() const -> vector<point_t>
{
    const Size flat_size = this->flat_size();
    vector<point_t> result(flat_size);
    for (Size i = 0; i < flat_size; ++i)
    {
        const auto comp_index = this->flat_to_tensor(i);
        result[i] = this->cartesian_product(comp_index);
    }

    return result;
}



template< class T, int rank>
auto
CartesianProductArray<T,rank>::
get_sub_product(const TensorIndex<rank-1> &index) const -> SubProduct
{
    SubProduct sub_data;
    for (int i=0; i<rank-1; ++i)
    {
        Assert(index[i]<rank, ExcIndexRange(index[i],0,rank));
        sub_data.copy_data_direction(i,data_[index[i]]);
    }

    return sub_data;
}

#if 0
template <class T, int rank>
CartesianProductArray<T,rank+1>
CartesianProductArray<T,rank>::
insert(const int index, const std::vector<T> &new_vector) const
{
    Assert(index<rank+1, ExcIndexRange(index,0,rank+1));

    TensorSize<rank+1> size;

    for (int i=0, j=0; i<rank+1; ++i)
    {
        if (i == index)
            size(i) = new_vector.size();
        else
        {
            size(i) = this->tensor_size()(j);
            ++j;
        }
    }

    CartesianProductArray<T,rank+1> product(size);

    for (int i=0, j=0; i<rank+1; ++i)
    {
        if (i == index)
            product.copy_data_direction(i,new_vector);
        else
        {
            product.copy_data_direction(i,this->data_[j]);
            ++j;
        }
    }
    return product;
}
#endif

template< class T, int rank>
void
CartesianProductArray<T,rank>::
print_info(LogStream &out) const
{
    for (int i = 0 ; i < rank ; i++)
        out << data_[i] << std::endl;
}





IGA_NAMESPACE_CLOSE

#ifndef NDEBUG
#include <igatools/utils/cartesian_product_array-inline.h>
#endif


#include <igatools/utils/cartesian_product_array.inst>
