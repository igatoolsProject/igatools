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



#ifndef CONCATENATED_FORWARD_ITERATOR_INLINE_H_
#define CONCATENATED_FORWARD_ITERATOR_INLINE_H_

#include <igatools/utils/concatenated_forward_iterator.h>
#include <igatools/base/exceptions.h>

IGA_NAMESPACE_OPEN


template <class ViewType,class DerivedClass>
ConcatenatedForwardIteratorData<ViewType,DerivedClass>::
ConcatenatedForwardIteratorData()
    :
    range_id_(IteratorState::invalid)
{}


template <class ViewType,class DerivedClass>
ConcatenatedForwardIteratorData<ViewType,DerivedClass>::
ConcatenatedForwardIteratorData(
    const std::vector<ViewType> &ranges,
    const Index index)
    :
    ranges_(ranges)
{
    const int n_ranges = ranges.size();
    Assert(n_ranges != 0 , ExcEmptyObject());

#ifndef NDEBUG
    for (int i = 0 ; i < n_ranges ; ++i)
    {
        //TODO (mm): maybe it is better to assert  ranges_[i].begin() < ranges_[i].end()
        Assert(ranges_[i].begin() != ranges_[i].end(),ExcInvalidIterator());
    }
#endif

    Assert(index == 0 || index == IteratorState::pass_the_end,ExcInvalidIterator());
    if (index == 0)
    {
        iterator_current_ = ranges_.front().begin();
        range_id_ = 0;
    }
    else if (index == IteratorState::pass_the_end)
    {
        iterator_current_ = ranges_.back().end();
        range_id_ = n_ranges - 1;
    }
}


template <class ViewType,class DerivedClass>
inline
auto
ConcatenatedForwardIteratorData<ViewType,DerivedClass>::
get_ranges() const -> std::vector<ViewType>
{
    return this->ranges_;
}




template <class ViewType,class DerivedClass>
inline
bool
ConcatenatedForwardIteratorData<ViewType,DerivedClass>::
operator==(const ConcatenatedForwardIteratorData<ViewType,DerivedClass> &it) const
{
    // check the equality of the size
    bool same_size = (this->ranges_.size() == it.ranges_.size());
    Assert(same_size,ExcMessage("Iterators are not comparable."));


    const int n_ranges = this->ranges_.size();
    bool ranges_are_equal = true;
    for (int i = 0 ; i < n_ranges ; ++i)
        if (this->ranges_[i].begin() != it.ranges_[i].begin() ||
            this->ranges_[i].end()   != it.ranges_[i].end())
        {
            ranges_are_equal = false;
            break;
        }
    Assert(ranges_are_equal,ExcMessage("Iterators are not comparable."));


    return (same_size &&
            ranges_are_equal &&
            this->range_id_ == it.range_id_ &&
            this->iterator_current_ == it.iterator_current_);
}

template <class ViewType,class DerivedClass>
inline
bool
ConcatenatedForwardIteratorData<ViewType,DerivedClass>::
operator!=(const ConcatenatedForwardIteratorData<ViewType,DerivedClass> &it) const
{
    return !(*this == it);
}

template <class ViewType,class DerivedClass>
inline
DerivedClass &
ConcatenatedForwardIteratorData<ViewType,DerivedClass>::
as_derived_class()
{
    return static_cast<DerivedClass &>(*this);
}

template <class ViewType,class DerivedClass>
inline
const DerivedClass &
ConcatenatedForwardIteratorData<ViewType,DerivedClass>::
as_derived_class() const
{
    return static_cast<const DerivedClass &>(*this);
}




template <class ViewType,class DerivedClass>
inline
DerivedClass &
ConcatenatedForwardIteratorData<ViewType,DerivedClass>::
operator++()
{
    if (this->range_id_ < this->ranges_.size()-1)
    {
        // if the current iterator is before the end, advance one position
        if (this->iterator_current_ != this->ranges_[this->range_id_].end())
            ++this->iterator_current_;

        // if the current iterator is already at the end of one iterator,
        // point to the first element of the next one
        if (this->iterator_current_ == this->ranges_[this->range_id_].end())
            this->iterator_current_ = this->ranges_[++this->range_id_].begin();
    }
    else
    {
        Assert(this->iterator_current_ != this->ranges_.back().end(),ExcIteratorPastEnd());
        ++this->iterator_current_;
    }

    return this->as_derived_class();
}



template <class ViewType,class DerivedClass>
inline
void
ConcatenatedForwardIteratorData<ViewType,DerivedClass>::
print_info(LogStream &out) const
{
    using std::endl;
    std::string tab("   ");

    out << "ConcatenatedForwardIteratorData infos:" << endl;
    out.push(tab);

    out << "Num. ranges = " << this->ranges_.size() << endl;
    int i = 0 ;
    for (const auto &r : this->ranges_)
    {
        out << "Range[" << i << "].begin() = " << &r.begin() << "   ";
        out << "Range[" << i << "].end() = " << &r.end();
        out << endl;
        ++i;
    }
    out << "range_id_ = " << this->range_id_ << endl;

    out.pop();
}


template <class ViewType,class DerivedClass>
inline
auto
ConcatenatedForwardIteratorData<ViewType,DerivedClass>::
operator*() const -> const value_type &
{
    Assert(this->iterator_current_ != this->ranges_.back().end(),ExcIteratorPastEnd());
    return *this->iterator_current_;
}

template <class ViewType,class DerivedClass>
inline
auto
ConcatenatedForwardIteratorData<ViewType,DerivedClass>::
operator->() const -> const value_type *
{
    return &(this->operator*());
}


template <class Iterator>
inline
ConcatenatedForwardConstIterator<Iterator>::
ConcatenatedForwardConstIterator(
    const std::vector<ConstView<Iterator>> &ranges,
    const Index index)
    :
    ConcatenatedForwardIteratorData<
    ConstView<Iterator>,
    ConcatenatedForwardConstIterator<Iterator>
    >(ranges,index)
{}




template <class ViewType>
inline
ConcatenatedForwardIterator<ViewType>::
ConcatenatedForwardIterator(
    const std::vector<ViewType> &ranges,
    const Index index)
    :
    ConcatenatedForwardIteratorData<
    ViewType,
    ConcatenatedForwardIterator<ViewType>
    >(ranges,index)
{}


template <class ViewType>
inline
auto
ConcatenatedForwardIterator<ViewType>::
operator*() -> value_type &
{
    Assert(this->iterator_current_ != this->ranges_.back().end(),ExcIteratorPastEnd());
    return *this->iterator_current_;
}

template <class ViewType>
inline
auto
ConcatenatedForwardIterator<ViewType>::
operator->() -> value_type *
{
    return &(this->operator*());
}

IGA_NAMESPACE_CLOSE


#endif // #ifndef CONCATENATED_FORWARD_ITERATOR_INLINE_H_
