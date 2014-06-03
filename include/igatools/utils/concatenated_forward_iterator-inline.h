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

template <class Iterator>
inline
ConcatenatedForwardIterator<Iterator>::
ConcatenatedForwardIterator()
{}


template <class Iterator>
inline
ConcatenatedForwardIterator<Iterator>::
ConcatenatedForwardIterator(
    const std::vector<std::pair<Iterator,Iterator>> &ranges,
    const Index index)
    :
    ranges_(ranges)
{
    const int n_ranges = ranges.size();
    Assert(n_ranges != 0 , ExcEmptyObject());

#ifndef NDEBUG
    for (int i = 0 ; i < n_ranges ; ++i)
        Assert(ranges_[i].first != ranges_[i].second,ExcInvalidIterator());
#endif

    Assert(index == 0 || index == IteratorState::pass_the_end,ExcInvalidIterator());
    if (index == 0)
    {
        iterator_current_ = ranges_.front().first;
        range_id_ = 0;
    }
    else if (index == IteratorState::pass_the_end)
    {
        iterator_current_ = ranges_.back().second;
        range_id_ = n_ranges - 1;
    }
    /*
        LogStream out;
        this->print_info(out);
    //*/
}



template <class Iterator>
inline
auto
ConcatenatedForwardIterator<Iterator>::
operator*() -> value_type &
{
    Assert(iterator_current_ != ranges_.back().second,ExcIteratorPastEnd());
    return *iterator_current_;
}

template <class Iterator>
inline
auto
ConcatenatedForwardIterator<Iterator>::
operator*() const -> const value_type &
{
    Assert(iterator_current_ != ranges_.back().second,ExcIteratorPastEnd());
    return *iterator_current_;
}

template <class Iterator>
inline
auto
ConcatenatedForwardIterator<Iterator>::
operator->() const -> const value_type *
{
    return &(this->operator*());
}


template <class Iterator>
inline
auto
ConcatenatedForwardIterator<Iterator>::
operator++() -> ConcatenatedForwardIterator<Iterator> &
{
    if (range_id_ < ranges_.size()-1)
    {
        // if the current iterator is before the end, advance one position
        if (iterator_current_ != ranges_[range_id_].second)
            ++iterator_current_;

        // if the current iterator is already at the end of one iterator,
        // point to the first element of the next one
        if (iterator_current_ == ranges_[range_id_].second)
            iterator_current_ = ranges_[++range_id_].first;
    }
    else
    {
        Assert(iterator_current_ != ranges_.back().second,ExcIteratorPastEnd());
        ++iterator_current_;
    }

    return *this;
}

template <class Iterator>
inline
bool
ConcatenatedForwardIterator<Iterator>::
operator==(const ConcatenatedForwardIterator<Iterator> &it) const
{
    // check the equality of the size
    bool same_size = (ranges_.size() == it.ranges_.size());
    Assert(same_size,ExcMessage("Iterators are not comparable."));


    const int n_ranges = ranges_.size();
    bool ranges_are_equal = true;
    for (int i = 0 ; i < n_ranges ; ++i)
        if (ranges_[i].first  != it.ranges_[i].first ||
            ranges_[i].second != it.ranges_[i].second)
        {
            ranges_are_equal = false;
            break;
        }
    Assert(ranges_are_equal,ExcMessage("Iterators are not comparable."));


    return (same_size &&
            ranges_are_equal &&
            range_id_ == it.range_id_ &&
            iterator_current_ == it.iterator_current_);
}

template <class Iterator>
inline
bool
ConcatenatedForwardIterator<Iterator>::
operator!=(const ConcatenatedForwardIterator<Iterator> &it) const
{
    return !(*this == it);
}


template <class Iterator>
inline
auto
ConcatenatedForwardIterator<Iterator>::
get_ranges() const -> std::vector<std::pair<Iterator,Iterator>>
{
    return ranges_;
}


template <class Iterator>
inline
void
ConcatenatedForwardIterator<Iterator>::
print_info(LogStream &out) const
{
    using std::endl;
    std::string tab("   ");

    out << "ConcatenatedForwardIterator infos:" << endl;
    out.push(tab);

    out << "Num. ranges = " << ranges_.size() << endl;
    int i = 0 ;
    for (const auto &r : ranges_)
    {
        out << "Range[" << i << "].first = " << &r.first << "   ";
        out << "Range[" << i << "].second = " << &r.second;
        out << endl;
        ++i;
    }
    out << "range_id_ = " << range_id_ << endl;

    out.pop();
}


IGA_NAMESPACE_CLOSE


#endif // #ifndef CONCATENATED_FORWARD_ITERATOR_INLINE_H_
