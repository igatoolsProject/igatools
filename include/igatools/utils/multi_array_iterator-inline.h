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


// QualityAssurance: martinelli, 31 Jan 2014


#ifndef MULTI_ARRAY_ITERATOR_INLINE_H_
#define MULTI_ARRAY_ITERATOR_INLINE_H_

#include <igatools/utils/multi_array_iterator.h>

IGA_NAMESPACE_OPEN


template <class Container>
inline
MultiArrayIterator<Container>::
MultiArrayIterator()
    :
    container_(nullptr),
    id_(IteratorState::invalid),
    stride_(0)
{}

template <class Container>
inline
MultiArrayIterator<Container>::
MultiArrayIterator(Container &container,const Index id,const Index stride)
    :
    container_(&container),
    id_(id),
    stride_(stride)
{
    //the iterator must not be built id the container is empty!
    Assert(container_->flat_size() > 0,ExcEmptyObject());

    Assert(id_ <= container_->flat_size(),ExcIteratorPastEnd());
    if (id_ == container_->flat_size())
        id_ = IteratorState::pass_the_end;
}

template <class Container>
inline
MultiArrayIterator<Container> &
MultiArrayIterator<Container>::
operator++()
{
    Assert(id_ != IteratorState::pass_the_end,ExcIteratorPastEnd());

    id_ += stride_ ;
    if (id_ >= container_->flat_size())
        id_ = IteratorState::pass_the_end;

    return (*this);
}

template <class Container>
inline
auto
MultiArrayIterator<Container>::
operator*() const -> const_reference
{
    Assert(id_ != IteratorState::pass_the_end,ExcIteratorPastEnd());
//    return container_->get_data()[id_];
    return (*container_)(id_);
}

template <class Container>
inline
auto
MultiArrayIterator<Container>::
operator*() -> reference
{
    Assert(id_ != IteratorState::pass_the_end,ExcIteratorPastEnd());
    return (*container_)(id_);
}

template <class Container>
inline
auto
MultiArrayIterator<Container>::
operator[](const Index i) const -> const_reference
{
    Assert(id_ + i*stride_ < container_->flat_size(),ExcIteratorPastEnd());
    return (*container_)(id_ + i*stride_);
}

template <class Container>
inline
auto
MultiArrayIterator<Container>::
operator[](const Index i) -> reference
{
    Assert(id_ + i*stride_ < container_->flat_size(),ExcIteratorPastEnd());
    return (*container_)(id_ + i*stride_);
}


template <class Container>
inline
auto
MultiArrayIterator<Container>::
operator->() const -> const pointer
{
    Assert(id_ != IteratorState::pass_the_end,ExcIteratorPastEnd());
//    return &container_->get_data()[id_];
    return &(*container_)(id_);
}

template <class Container>
inline
auto
MultiArrayIterator<Container>::
operator->() -> pointer
{
    Assert(id_ != IteratorState::pass_the_end,ExcIteratorPastEnd());
    return &(*container_)(id_);
}

template <class Container>
inline
bool
MultiArrayIterator<Container>::
operator==(const MultiArrayIterator<Container> &it) const
{
    return (id_ == it.id_ && stride_ == it.stride_ && container_ == it.container_)?true:false;
}

template <class Container>
inline
bool
MultiArrayIterator<Container>::
operator!=(const MultiArrayIterator<Container> &it) const
{
    return (id_ != it.id_ || stride_ != it.stride_ || container_ != it.container_)?true:false;
}


template <class Container>
inline
bool
MultiArrayIterator<Container>::
operator<(const MultiArrayIterator<Container> &it) const
{
    // check if the iterators are comparable
    Assert(stride_ == it.stride_ && container_ == it.container_,
           ExcInvalidIterator());

    bool ret = false ;
    if (id_ != IteratorState::pass_the_end)
    {
        if (it.id_ == IteratorState::pass_the_end || id_ < it.id_)
            ret = true;
    }
    return ret;
}

template <class Container>
inline
bool
MultiArrayIterator<Container>::
operator<=(const MultiArrayIterator<Container> &it) const
{
    return ((*this) < it || (*this) == it)?true:false;
}


template <class Container>
inline
MultiArrayIterator<Container>
MultiArrayIterator<Container>::
operator+(const Index n) const
{
    Assert(id_ + n*stride_ < container_->flat_size(),ExcIteratorPastEnd());
    return MultiArrayIterator<Container>(*container_,id_ + n*stride_,stride_);
}

template <class Container>
inline
MultiArrayIterator<Container>
MultiArrayIterator<Container>::
operator-(const Index n) const
{
    Assert(id_ - n*stride_ >=0, ExcLowerRange(id_ - n*stride_,0));
    return MultiArrayIterator<Container>(*container_,id_ - n*stride_,stride_);
}

#if 0
template <class Container>
inline
auto
MultiArrayIterator<Container>::
operator-(const MultiArrayIterator<Container> &a) const -> difference_type
{
    // check if the iterators are comparable
    Assert(a.stride_ == stride_ && a.container_ == container_,
           ExcInvalidIterator());

    Assert(a <= (*this), ExcInvalidIterator());

    const int position_difference = id_ - a.id_;

    const difference_type n = position_difference / a.stride_;

    Assert(a+n == (*this), ExcMessage("Iterator a cannot advance to (*this)."));
    return n;
}
#endif

IGA_NAMESPACE_CLOSE



#endif // MULTI_ARRAY_ITERATOR_INLINE_H_
