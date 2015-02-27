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

// QualityAssurance: martinelli, 31 Jan 2014

#ifndef MULTI_ARRAY_INLINE_H_
#define MULTI_ARRAY_INLINE_H_


#include<igatools/utils/multi_array.h>

IGA_NAMESPACE_OPEN

template<class STLContainer, int rank>
inline
MultiArray<STLContainer,rank>::
MultiArray()
    :
    TensorSizedContainer<rank>(0)
{}

template<class STLContainer, int rank>
inline
MultiArray<STLContainer,rank>::
MultiArray(const int dim)
    :
    TensorSizedContainer<rank>(dim)
{}


template<class STLContainer, int rank>
inline
MultiArray<STLContainer,rank>::
MultiArray(const TensorSize<rank> &dim)
    :
    TensorSizedContainer<rank>(dim)
{}


template<class STLContainer, int rank>
inline
auto
MultiArray<STLContainer,rank>::
get_data() const -> const STLContainer &
{
    return this->data_;
}


template<class STLContainer, int rank>
inline
auto
MultiArray<STLContainer,rank>::
operator[](const Index i) -> reference
{
    Assert((0<=i)  &&(i<this->flat_size()),
    ExcIndexRange(i,0,this->flat_size()));
    return this->data_[i];
}

template<class STLContainer, int rank>
inline
auto
MultiArray<STLContainer,rank>::
operator[](const Index i) const -> const_reference
{
    Assert((0<=i) &&(i<this->flat_size()),
           ExcIndexRange(i,0,this->flat_size()));
    return this->data_[i];
}
//*/

template<class STLContainer, int rank>
inline
auto
MultiArray<STLContainer,rank>::
operator()(const TensorIndex<rank> &tensor_index) -> reference
{
#ifndef NDEBUG
    const TensorSize<rank> tensor_size = this->tensor_size();

    for (int i = 0 ; i < rank ; ++i)
        Assert((tensor_index[i] >= 0) &&
        (tensor_index[i] < tensor_size[i]),
        ExcIndexRange(tensor_index[i],0,tensor_size[i]));
#endif

    Assert((this->tensor_to_flat(tensor_index) >= 0) &&
    (this->tensor_to_flat(tensor_index) < this->flat_size()),
    ExcIndexRange(this->tensor_to_flat(tensor_index),0,this->flat_size()));

    return this->data_[this->tensor_to_flat(tensor_index)];
}


template<class STLContainer, int rank>
inline
auto
MultiArray<STLContainer,rank>::
operator()(const TensorIndex<rank> &tensor_index) const -> const_reference
{
    return const_cast<MultiArray<STLContainer,rank> &>(*this)(tensor_index);
}




template<class STLContainer, int rank>
inline
auto
MultiArray<STLContainer,rank>::
cbegin() const -> const_iterator
{
    return const_iterator(*this,0);
}

template<class STLContainer, int rank>
inline
auto
MultiArray<STLContainer,rank>::
cend() const -> const_iterator
{
    return const_iterator(*this,IteratorState::pass_the_end);
}

template<class STLContainer, int rank>
inline
auto
MultiArray<STLContainer,rank>::
begin() const -> const_iterator
{
    return cbegin();
}

template<class STLContainer, int rank>
inline
auto
MultiArray<STLContainer,rank>::
end() const -> const_iterator
{
    return cend();
}


template<class STLContainer, int rank>
inline
auto
MultiArray<STLContainer,rank>::
begin() -> iterator
{
    return iterator(*this,0);
}

template<class STLContainer, int rank>
inline
auto
MultiArray<STLContainer,rank>::
end() -> iterator
{
    return iterator(*this,IteratorState::pass_the_end);
}




template<class STLContainer, int rank>
inline
void
MultiArray<STLContainer,rank>::
fill_progression(const value_type &init)
{
    auto val = init;
    for (auto &d : data_)
        d = val++;
}



template<class STLContainer, int rank>
inline
void
MultiArray<STLContainer,rank>::
fill(const value_type &value)
{
    std::fill(this->data_.begin(),this->data_.end(),value);
}



template<class STLContainer, int rank>
inline
ContainerView<MultiArray<STLContainer,rank> >
MultiArray<STLContainer,rank>::
get_view()
{
    return ContainerView<MultiArray<STLContainer,rank>>(this->begin(),this->end());
}

template<class STLContainer, int rank>
inline
ConstContainerView<MultiArray<STLContainer,rank> >
MultiArray<STLContainer,rank>::
get_const_view() const
{
    return ConstContainerView<MultiArray<STLContainer,rank>>(this->cbegin(),this->cend());
}


template<class STLContainer, int rank>
inline
ContainerView<STLContainer>
MultiArray<STLContainer,rank>::
get_flat_view()
{
    return ContainerView<STLContainer>(this->data_.begin(),this->data_.end());
}

template<class STLContainer, int rank>
inline
ConstContainerView<STLContainer>
MultiArray<STLContainer,rank>::
get_flat_const_view() const
{
    return ConstContainerView<STLContainer>(this->data_.cbegin(),this->data_.cend());
}


template<class STLContainer, int rank>
inline
void
MultiArray<STLContainer,rank>::
print_info(LogStream &out) const
{
    TensorSizedContainer<rank>::print_info(out);
    out << std::endl;
    t_print_info<STLContainer>(out);
}


template<class STLContainer, int rank>
inline
void
MultiArray<STLContainer,rank>::
resize(const Size dim)
{
    this->resize(TensorSize<rank>(dim));
}


template<class STLContainer, int rank>
inline
void
MultiArray<STLContainer,rank>::
resize(const TensorSize<rank> &dim)
{
    this->reset_size(dim);
    this->data_.resize(this->flat_size());
}


template<class STLContainer, int rank>
inline
void
MultiArray<STLContainer,rank>::
resize(const TensorSize<rank> &dim, const typename STLContainer::value_type &val)
{
    this->reset_size(dim);
    this->data_.resize(this->flat_size(),val);
}



IGA_NAMESPACE_CLOSE


#endif // #ifndef MULTI_ARRAY_INLINE_H_


