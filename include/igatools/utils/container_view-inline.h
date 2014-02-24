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



#ifndef CONTAINER_VIEW_INLINE_H_
#define CONTAINER_VIEW_INLINE_H_

#include <igatools/utils/container_view.h>



IGA_NAMESPACE_OPEN




template <class Container>
inline
ContainerView<Container>::
ContainerView(const iterator begin, const iterator end)
    : begin_(begin), end_(end)
{
    Assert(begin_ <= end_, ExcInvalidIterator());
}

template <class Container>
inline
auto
ContainerView<Container>::
begin() -> iterator
{
    return begin_;
}

template <class Container>
inline
auto
ContainerView<Container>::
begin() const -> const_iterator
{
    return begin_;
}

template <class Container>
inline
auto
ContainerView<Container>::
end() -> iterator
{
    return end_;
}

template <class Container>
inline
auto
ContainerView<Container>::
end() const -> const_iterator
{
    return end_;
}

template <class Container>
inline
auto
ContainerView<Container>::
operator[](const Index n) -> reference
{

    Assert(begin_+n < end_, ExcIteratorPastEnd());
    return begin_[n];
}

template <class Container>
inline
auto
ContainerView<Container>::
operator[](const Index n) const -> const_reference
{
    Assert(begin_+n < end_, ExcIteratorPastEnd());
    return begin_[n];
}

template <class Container>
inline
ConstContainerView<Container>::
ConstContainerView(const const_iterator begin, const const_iterator end)
    : begin_(begin), end_(end)
{
    Assert(begin_ <= end_, ExcInvalidIterator());
}


template <class Container>
inline
auto
ConstContainerView<Container>::
begin() const -> const_iterator
{
    return begin_;
}

template <class Container>
inline
auto
ConstContainerView<Container>::
end() const -> const_iterator
{
    return end_;
}

template <class Container>
inline
auto
ConstContainerView<Container>::
operator[](const Index n) const -> const_reference
{
    Assert(begin_+n < end_, ExcIteratorPastEnd());
    return begin_[n];
}


IGA_NAMESPACE_CLOSE

#endif //#ifndef CONTAINER_VIEW_INLINE_H_


