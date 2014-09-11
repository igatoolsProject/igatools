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

#include <algorithm>

IGA_NAMESPACE_OPEN


template <class Container>
inline
MultiArrayIteratorBase<Container>::
MultiArrayIteratorBase()
    :
    container_(nullptr),
    id_(IteratorState::invalid),
    stride_(0)
{}

template <class Container>
inline
MultiArrayIteratorBase<Container>::
MultiArrayIteratorBase(Container &container,const Index id,const Index stride)
    :
    container_(&container),
    id_(id),
    stride_(stride)
{
    Assert(container_ != nullptr,ExcNullPtr());
    Assert(id_ >= 0 || IteratorState::pass_the_end,ExcInvalidState());

#ifndef NDEBUG
    if (container_->flat_size() != 0)
    {
        Assert(stride_ > 0 && stride <= container_->flat_size(),
               ExcIndexRange(stride_,1,container_->flat_size()+1));
    }
#endif


    Assert(id_ <= container_->flat_size(),ExcIteratorPastEnd());
    if (container_->flat_size() > 0)
    {
        if (id_ == container_->flat_size())
            id_ = IteratorState::pass_the_end;
    }
    else
        id_ = IteratorState::invalid;
}



template <class Container>
inline
MultiArrayIteratorBase<Container> &
MultiArrayIteratorBase<Container>::
operator++()
{
    Assert(id_ != IteratorState::pass_the_end,ExcIteratorPastEnd());
    Assert(id_ != IteratorState::invalid,ExcInvalidIterator());

    id_ += stride_ ;

    if (id_ >= container_->flat_size())
        id_ = IteratorState::pass_the_end;
    return (*this);
}

template <class Container>
inline
MultiArrayIteratorBase<Container>
MultiArrayIteratorBase<Container>::
operator++(int)
{
    MultiArrayIteratorBase<Container> tmp(*this);
    operator++();
    return tmp;
}


template <class Container>
inline
auto
MultiArrayIteratorBase<Container>::
operator*() const -> const reference
{
    Assert(id_ != IteratorState::pass_the_end,ExcIteratorPastEnd());
    Assert(id_ != IteratorState::invalid,ExcInvalidIterator());
    return (*container_)[id_];
}

template <class Container>
inline
auto
MultiArrayIteratorBase<Container>::
operator*() -> reference
{
    Assert(id_ != IteratorState::pass_the_end,ExcIteratorPastEnd());
    Assert(id_ != IteratorState::invalid,ExcInvalidIterator());
    return (*container_)[id_];
}

template <class Container>
inline
auto
MultiArrayIteratorBase<Container>::
operator[](const Index i) const -> const reference
{
    Assert(i >= 0,ExcLowerRange(i,0));
    Assert(id_ + i*stride_ < container_->flat_size(),ExcIteratorPastEnd());
    Assert(id_ != IteratorState::invalid,ExcInvalidIterator());
    return (*container_)[id_ + i*stride_];
}

template <class Container>
inline
auto
MultiArrayIteratorBase<Container>::
operator[](const Index i) -> reference
{
    Assert(i >= 0,ExcLowerRange(i,0));
    Assert(id_ + i*stride_ < container_->flat_size(),ExcIteratorPastEnd());
    Assert(id_ != IteratorState::invalid,ExcInvalidIterator());
    return (*container_)[id_ + i *stride_];
}


template <class Container>
inline
auto
MultiArrayIteratorBase<Container>::
operator->() const -> const pointer
{
    Assert(id_ != IteratorState::pass_the_end,ExcIteratorPastEnd());
    Assert(id_ != IteratorState::invalid,ExcInvalidIterator());
//    return &container_->get_data()[id_];
    return &(*container_)[id_];
}

template <class Container>
inline
auto
MultiArrayIteratorBase<Container>::
operator->() -> pointer
{
    Assert(id_ != IteratorState::pass_the_end,ExcIteratorPastEnd());
    Assert(id_ != IteratorState::invalid,ExcInvalidIterator());
    return &(*container_)[id_];
}

template <class Container>
inline
bool
MultiArrayIteratorBase<Container>::
operator==(const MultiArrayIteratorBase<Container> &it) const
{
    return (id_ == it.id_ && stride_ == it.stride_ && container_ == it.container_)?true:false;
}

template <class Container>
inline
bool
MultiArrayIteratorBase<Container>::
operator!=(const MultiArrayIteratorBase<Container> &it) const
{
    return (id_ != it.id_ || stride_ != it.stride_ || container_ != it.container_)?true:false;
}


template <class Container>
inline
bool
MultiArrayIteratorBase<Container>::
operator<(const MultiArrayIteratorBase<Container> &it) const
{
    Assert(id_ != IteratorState::invalid,ExcInvalidIterator());
    Assert(it.id_ != IteratorState::invalid,ExcInvalidIterator());

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
MultiArrayIteratorBase<Container>::
operator<=(const MultiArrayIteratorBase<Container> &it) const
{
    return ((*this) < it || (*this) == it)?true:false;
}


template <class Container>
inline
MultiArrayIteratorBase<Container>
MultiArrayIteratorBase<Container>::
operator+(const Index n) const
{
    Assert(n >= 0,ExcLowerRange(n,0));
//    Assert(id_ + n*stride_ < container_->flat_size(),ExcIteratorPastEnd());
    Assert(id_ != IteratorState::pass_the_end,ExcInvalidIterator());
    Assert(id_ != IteratorState::invalid,ExcInvalidIterator());

    const Index pos = id_ + n*stride_;
    if (pos < container_->flat_size())
        return MultiArrayIteratorBase<Container>(*container_,pos,stride_);
    else
        return MultiArrayIteratorBase<Container>(*container_,IteratorState::pass_the_end,stride_);
}

template <class Container>
inline
MultiArrayIteratorBase<Container>
MultiArrayIteratorBase<Container>::
operator-(const Index n) const
{
    Assert(n >= 0,ExcLowerRange(n,0));
    Assert(id_ != IteratorState::pass_the_end,ExcInvalidIterator());
    Assert(id_ != IteratorState::invalid,ExcInvalidIterator());
    Assert(id_ - n*stride_ >=0, ExcLowerRange(id_ - n*stride_,0));
    Assert(id_ != IteratorState::invalid,ExcInvalidIterator());
    return MultiArrayIteratorBase<Container>(*container_,id_ - n*stride_,stride_);
}

template <class Container>
inline
auto
MultiArrayIteratorBase<Container>::
operator-(const MultiArrayIteratorBase<Container> &a) const -> difference_type
{
    // check if the iterators are comparable
    Assert(a.stride_ == stride_ && a.container_ == container_,
           ExcInvalidIterator());

    Assert(a <= (*this), ExcInvalidIterator());

    difference_type position_difference;
    if (id_ != IteratorState::pass_the_end)
    {
        position_difference = id_ - a.id_;
    }
    else
    {
        if (a.id_ != IteratorState::pass_the_end)
            position_difference = container_->flat_size() - a.id_;
        else
            position_difference = 0;
    }

    Assert(position_difference >= 0,ExcLowerRange(position_difference,0));

    const difference_type n = position_difference / a.stride_;
    Assert(n >= 0,ExcLowerRange(n,0));

    Assert(a+n == (*this), ExcMessage("Iterator a cannot advance to (*this)."));
    return n;
}



template <class Container>
inline
Index
MultiArrayIteratorBase<Container>::
get_id() const
{
    return id_;
}

template <class Container>
inline
Index
MultiArrayIteratorBase<Container>::
get_stride() const
{
    return stride_;
}




template <class Container>
inline
MultiArrayConstIterator<Container>::
MultiArrayConstIterator(const Container &container,const Index id,const Index stride)
    :
    MultiArrayIteratorBase<Container>(const_cast<Container &>(container),id,stride)
{}


template <class Container>
inline
MultiArrayIterator<Container>::
MultiArrayIterator(Container &container,const Index id,const Index stride)
    :
    MultiArrayConstIterator<Container>(container,id,stride)
{}









IGA_NAMESPACE_CLOSE



#endif // MULTI_ARRAY_ITERATOR_INLINE_H_
