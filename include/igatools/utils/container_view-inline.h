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

template <class IteratorType>
ViewData<IteratorType>::
ViewData(const IteratorType begin, const IteratorType end)
    :
    begin_(begin),
    end_(end)
{
    Assert(begin_ < end_, ExcInvalidIterator());
//    Assert(begin_ < end_ || begin_ == end_, ExcInvalidIterator());
}

template <class IteratorType>
inline
Size
ViewData<IteratorType>::
get_num_entries() const
{
    return end_ - begin_;
}


template <class IteratorType>
inline
ViewData<IteratorType> &
ViewData<IteratorType>::
operator=(const ViewData<IteratorType> &view_data)
{
    if (this != &view_data)
    {
        begin_ = view_data.begin_;
        end_   = view_data.end_;
    }
    return *this;
}


template <class Iterator,class ConstIterator>
inline
View<Iterator,ConstIterator>::
View(const iterator begin, const iterator end)
    :
    ViewData<Iterator>(begin,end)
{}


template <class Iterator,class ConstIterator>
inline
auto
View<Iterator,ConstIterator>::
begin() -> iterator
{
    return this->begin_;
}

template <class Iterator,class ConstIterator>
inline
auto
View<Iterator,ConstIterator>::
begin() const -> const_iterator
{
    return ConstIterator(this->begin_);
}

template <class Iterator,class ConstIterator>
inline
auto
View<Iterator,ConstIterator>::
cbegin() const -> const_iterator
{
    return ConstIterator(this->begin_);
}

template <class Iterator,class ConstIterator>
inline
auto
View<Iterator,ConstIterator>::
end() -> iterator
{
    return this->end_;
}

template <class Iterator,class ConstIterator>
inline
auto
View<Iterator,ConstIterator>::
end() const -> const_iterator
{
    return ConstIterator(this->end_);
}

template <class Iterator,class ConstIterator>
inline
auto
View<Iterator,ConstIterator>::
cend() const -> const_iterator
{
    return ConstIterator(this->end_);
}

template <class Iterator,class ConstIterator>
inline
auto
View<Iterator,ConstIterator>::
operator[](const Index n) -> reference
{
#ifndef NDEBUG
    auto tmp = this->begin_;
    Assert(tmp+n < this->end_, ExcIteratorPastEnd());
#endif
    return this->begin_[n];
}

template <class Iterator,class ConstIterator>
inline
auto
View<Iterator,ConstIterator>::
operator[](const Index n) const -> const reference
{
#ifndef NDEBUG
    auto tmp = this->begin_;
    Assert(tmp+n < this->end_, ExcIteratorPastEnd());
#endif
//    return const_cast<const Iterator &>(this->begin_)[n];
    return this->begin_[n];
}

template <class Iterator,class ConstIterator>
inline
ConstView<Iterator,ConstIterator>::
ConstView(const const_iterator begin, const const_iterator end)
    :
    ViewData<ConstIterator>(begin,end)
{}


template <class Iterator,class ConstIterator>
inline
ConstView<Iterator,ConstIterator>::
ConstView(const View<Iterator,ConstIterator> &view)
    :
    ViewData<ConstIterator>(view.begin(),view.end())
{
//  this->begin_ = view.begin();
//  this->end_   = view.end();
}

template <class Iterator,class ConstIterator>
inline
auto
ConstView<Iterator,ConstIterator>::
begin() const -> const_iterator
{
    return this->begin_;
}

template <class Iterator,class ConstIterator>
inline
auto
ConstView<Iterator,ConstIterator>::
end() const -> const_iterator
{
    return this->end_;
}

template <class Iterator,class ConstIterator>
inline
auto
ConstView<Iterator,ConstIterator>::
operator[](const Index n) const -> const reference
{
    Assert(this->begin_+n < this->end_, ExcIteratorPastEnd());
    return this->begin_[n];
}


IGA_NAMESPACE_CLOSE

#endif //#ifndef CONTAINER_VIEW_INLINE_H_


