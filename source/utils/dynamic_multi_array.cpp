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


#include <igatools/utils/dynamic_multi_array.h>
#include <igatools/utils/multi_array_utils.h>

using std::vector;

IGA_NAMESPACE_OPEN



template<class T, int rank>
DynamicMultiArray<T,rank>::
DynamicMultiArray(const Size dim)
    :
    MultiArray<vector<T>,rank>(dim)
{
    this->data_.resize(this->flat_size());
}


template<class T, int rank>
DynamicMultiArray<T,rank>::
DynamicMultiArray(const TensorSize<rank> &dim)
    :
    MultiArray<vector<T>,rank>(dim)
{
    this->data_.resize(this->flat_size());
}


template<class T, int rank>
void
DynamicMultiArray<T,rank>::
resize(const Size dim)
{
    TensorSizedContainer<rank>::reset_size(TensorSize<rank>(dim));
    this->data_.resize(this->flat_size());
}


template<class T, int rank>
void
DynamicMultiArray<T,rank>::
resize(const TensorSize<rank> &dim)
{
    TensorSizedContainer<rank>::reset_size(dim);
    this->data_.resize(this->flat_size());
}







template<class T, int rank>
std::vector<T>
DynamicMultiArray<T,rank>::
get_flat_view(const TensorIndex<rank> &start, const TensorIndex<rank> &inc) const
{
    TensorSize<rank> loc_size;
    for (int i = 0 ; i < rank ; ++i)
        loc_size(i) = inc[i];

    auto loc_weight = MultiArrayUtils<rank>::compute_weight(loc_size);
    const Size size = loc_size.flat_size();
    Assert(size == MultiArrayUtils<rank>::size(inc),
           ExcDimensionMismatch(size,MultiArrayUtils<rank>::size(inc)));

    std::vector<T> result(size);
    for (int i = 0; i < size; ++i)
    {
        auto tensor_index = MultiArrayUtils<rank>::flat_to_tensor_index(i,loc_weight);
        tensor_index += start;
        result[i] = (*this)(tensor_index);
    }

    return result;
}



template<class T, int rank>
inline
Conditional< (rank > 0),DynamicMultiArray<T,rank-1>,DynamicMultiArray<T,0> >
DynamicMultiArray<T,rank>::
get_slice(const int direction, const Index index) const
{
//  LogStream out ;
//  out << *this << std::endl ;

    const auto tensor_size = this->tensor_size();

    Assert(direction >= 0 && direction < rank,
           ExcIndexRange(direction,0,rank));
    Assert(index >= 0 && index < tensor_size(direction),
           ExcIndexRange(index,0,tensor_size(direction)));

    const int rank_slice = (rank>0)?rank-1:0;
    TensorSize<rank_slice> sizes_slice;
    for (Index i = 0 ; i < direction ; ++i)
        sizes_slice(i) = tensor_size(i);

    for (Index i = direction+1 ; i < rank ; ++i)
        sizes_slice(i-1) = tensor_size(i);


    TensorIndex<rank> tensor_id;
    tensor_id[direction] = index;

    DynamicMultiArray<T,rank_slice> slice(sizes_slice);
    const Size n_entries_slice = slice.flat_size();
    for (Index flat_id_slice = 0 ; flat_id_slice < n_entries_slice ; ++flat_id_slice)
    {
        const auto tensor_id_slice = slice.flat_to_tensor(flat_id_slice);

        for (Index i = 0 ; i < direction ; ++i)
            tensor_id[i] = tensor_id_slice[i];

        for (Index i = direction+1 ; i < rank ; ++i)
            tensor_id[i] = tensor_id_slice[i-1];

        slice(flat_id_slice) = this->operator()(tensor_id);
    }


    return slice;
}



template<class T, int rank>
inline
void
DynamicMultiArray<T,rank>::
copy_slice(const int direction, const Index index,
           const Conditional<(rank > 0),DynamicMultiArray<T,rank-1>,DynamicMultiArray<T,0>> &slice)
{
    Assert(direction >= 0 && direction < rank,
           ExcIndexRange(direction,0,rank));


#ifndef NDEBUG
    const auto tensor_size = this->tensor_size();
    const auto sizes_slice = slice.tensor_size();
    for (Index i = 0 ; i < direction ; ++i)
        Assert(tensor_size(i) == sizes_slice(i),
               ExcDimensionMismatch(tensor_size(i),sizes_slice(i)));


    for (Index i = direction+1 ; i < rank ; ++i)
        Assert(tensor_size(i) == sizes_slice(i-1),
               ExcDimensionMismatch(tensor_size(i),sizes_slice(i-1)));
#endif

    TensorIndex<rank> tensor_id;
    tensor_id[direction] = index;

    const Size n_entries_slice = slice.flat_size();
    for (Index flat_id_slice = 0 ; flat_id_slice < n_entries_slice ; ++flat_id_slice)
    {
        const auto tensor_id_slice = slice.flat_to_tensor(flat_id_slice);

        for (Index i = 0 ; i < direction ; ++i)
            tensor_id[i] = tensor_id_slice[i];

        for (Index i = direction+1 ; i < rank ; ++i)
            tensor_id[i] = tensor_id_slice[i-1];

        this->operator()(tensor_id) = slice(flat_id_slice);
    }


//  AssertThrow(false,ExcNotImplemented());
}



template<class T, int rank>
LogStream &operator<<(LogStream &out, const DynamicMultiArray<T,rank> &data)
{
    Size size = data.flat_size();

    for (Index i = 0 ; i < size ; ++i)
        out << data(i) << " ";
    return out;
}


template<class T, int rank>
EnableIf<std::is_floating_point<T>::value  &&(rank>=0),DynamicMultiArray<T,rank> >
operator*(const Real a, const DynamicMultiArray<T,rank> &B)
{
    DynamicMultiArray<T,rank> res(B);

    for (auto & r : res)
        r *= a;

    return res;
}


template<class T, int rank>
EnableIf< std::is_floating_point<T>::value  &&(rank>=0),DynamicMultiArray<T,rank>>
        operator+(const DynamicMultiArray<T,rank> &A, const DynamicMultiArray<T,rank> &B)
{
    for (Index i = 0 ; i < rank ; ++i)
        Assert(A.tensor_size()[i] == B.tensor_size()[i],
               ExcDimensionMismatch(A.tensor_size()[i],B.tensor_size()[i]));

    DynamicMultiArray<T,rank> res(A);

    const Size n_entries = res.flat_size();
    for (Index i = 0 ; i < n_entries ; ++i)
        res(i) += B(i);

    return res;
}


IGA_NAMESPACE_CLOSE



#include <igatools/utils/dynamic_multi_array.inst>
