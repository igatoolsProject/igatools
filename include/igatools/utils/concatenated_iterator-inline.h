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



#ifndef CONCATENATED_ITERATOR_INLINE_H_
#define CONCATENATED_ITERATOR_INLINE_H_

#include <igatools/utils/concatenated_iterator.h>
#include <igatools/base/exceptions.h>

IGA_NAMESPACE_OPEN


template <class ViewType,class DerivedClass>
ConcatenatedIteratorData<ViewType,DerivedClass>::
ConcatenatedIteratorData()
    :
    range_id_(IteratorState::invalid)
{}


template <class ViewType,class DerivedClass>
ConcatenatedIteratorData<ViewType,DerivedClass>::
ConcatenatedIteratorData(
    const SafeSTLVector<ViewType> &ranges,
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
        Assert(ranges_[i].begin() < ranges_[i].end(),
               ExcInvalidIterator());
//        Assert(ranges_[i].begin() < ranges_[i].end() || ranges_[i].begin() == ranges_[i].end(),
//              ExcInvalidIterator());
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
ConcatenatedIteratorData<ViewType,DerivedClass>::
get_ranges() const -> SafeSTLVector<ViewType>
{
    return this->ranges_;
}

template <class ViewType,class DerivedClass>
inline
int
ConcatenatedIteratorData<ViewType,DerivedClass>::
get_range_id() const
{
    return range_id_;
}

template <class ViewType,class DerivedClass>
inline
auto
ConcatenatedIteratorData<ViewType,DerivedClass>::
get_iterator_current() const -> Iterator
{
    return iterator_current_;
}



template <class ViewType,class DerivedClass>
inline
bool
ConcatenatedIteratorData<ViewType,DerivedClass>::
is_comparable(const ConcatenatedIteratorData<ViewType,DerivedClass> &it) const
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

    return (same_size && ranges_are_equal);
}

template <class ViewType,class DerivedClass>
inline
bool
ConcatenatedIteratorData<ViewType,DerivedClass>::
operator==(const ConcatenatedIteratorData<ViewType,DerivedClass> &it) const
{
    Assert(this->is_comparable(it), ExcMessage("Iterators are not comparable."));

    return (this->range_id_ == it.range_id_ &&
            this->iterator_current_ == it.iterator_current_);
}

template <class ViewType,class DerivedClass>
inline
bool
ConcatenatedIteratorData<ViewType,DerivedClass>::
operator<(const ConcatenatedIteratorData<ViewType,DerivedClass> &it) const
{
    Assert(this->is_comparable(it), ExcMessage("Iterators are not comparable."));

    return (this->range_id_ < it.range_id_ ||
            (this->range_id_ == it.range_id_ && this->iterator_current_ < it.iterator_current_));
}


template <class ViewType,class DerivedClass>
bool
ConcatenatedIteratorData<ViewType,DerivedClass>::
operator<=(const ConcatenatedIteratorData<ViewType,DerivedClass> &a) const
{
    return ((*this) == a || (*this) < a);
}

template <class ViewType,class DerivedClass>
inline
bool
ConcatenatedIteratorData<ViewType,DerivedClass>::
operator!=(const ConcatenatedIteratorData<ViewType,DerivedClass> &it) const
{
    return !(*this == it);
}

template <class ViewType,class DerivedClass>
inline
DerivedClass &
ConcatenatedIteratorData<ViewType,DerivedClass>::
as_derived_class()
{
    return static_cast<DerivedClass &>(*this);
}

template <class ViewType,class DerivedClass>
inline
const DerivedClass &
ConcatenatedIteratorData<ViewType,DerivedClass>::
as_derived_class() const
{
    return static_cast<const DerivedClass &>(*this);
}




template <class ViewType,class DerivedClass>
inline
DerivedClass &
ConcatenatedIteratorData<ViewType,DerivedClass>::
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
DerivedClass &
ConcatenatedIteratorData<ViewType,DerivedClass>::
operator+(const int n)
{
    Assert(n>=0,ExcLowerRange(n,0));
    for (int i=0 ; i < n ; ++i, ++(*this));

    return this->as_derived_class();
}


template <class ViewType,class DerivedClass>
inline
Size
ConcatenatedIteratorData<ViewType,DerivedClass>::
operator-(const ConcatenatedIteratorData<ViewType,DerivedClass> &a) const
{
    Assert(this->is_comparable(a), ExcMessage("Iterators are not comparable."));
    Assert(a <= (*this),ExcInvalidIterator());

    Size n_entries = 0;
    for (Index i = a.range_id_ ; i < this->range_id_ ; ++i)
        n_entries += this->ranges_[i].get_num_entries();

    n_entries -= a.iterator_current_ - a.ranges_[a.range_id_].begin();

    for (auto it = ranges_[range_id_].begin() ; it != this->iterator_current_ ; ++it)
        ++n_entries;

    return n_entries;
}


template <class ViewType,class DerivedClass>
inline
void
ConcatenatedIteratorData<ViewType,DerivedClass>::
get_range_id_and_entry_id_in_range(const Index id, Index &rng_id, Index &entry_id_rng) const
{
    Assert(!this->ranges_.empty(),ExcEmptyObject());

//        using std::cout;
//        using std::endl;

    entry_id_rng = 0;

    // find the view that holds the data
    Index id_first = 0;
    Index id_last = - 1;
    for (const auto rng : this->ranges_)
    {
        id_first = id_last + 1;
        id_last += rng.get_num_entries() ;

//          cout << "id_first=" << id_first << "   id_last=" <<id_last << "   num entries = " << rng.get_num_entries() << endl;
        if (id >= id_first && id <= id_last)
        {
            entry_id_rng = id - id_first;
            break;
        }

        rng_id += 1;
    }
}

template <class ViewType,class DerivedClass>
inline
auto
ConcatenatedIteratorData<ViewType,DerivedClass>::
get_entry_const_reference(const Index id) const -> const typename ViewType::reference
{
    Index rng_id = 0 ;
    Index dof_id_rng = 0;
    this->get_range_id_and_entry_id_in_range(id,rng_id,dof_id_rng);

    return this->ranges_[rng_id][dof_id_rng];
}

template <class ViewType,class DerivedClass>
inline
auto
ConcatenatedIteratorData<ViewType,DerivedClass>::
get_entry_reference(const Index id) -> typename ViewType::reference
{
    Index rng_id = 0 ;
    Index dof_id_rng = 0;
    this->get_range_id_and_entry_id_in_range(id,rng_id,dof_id_rng);

    return this->ranges_[rng_id][dof_id_rng];
}


template <class ViewType,class DerivedClass>
inline
void
ConcatenatedIteratorData<ViewType,DerivedClass>::
print_info(LogStream &out) const
{
    using std::endl;
    std::string tab("   ");

    out << "ConcatenatedIteratorData infos:" << endl;
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
ConcatenatedIteratorData<ViewType,DerivedClass>::
operator*() const -> const value_type &
{
    Assert(this->iterator_current_ != this->ranges_.back().end(),ExcIteratorPastEnd());
    return *this->iterator_current_;
}

template <class ViewType,class DerivedClass>
inline
auto
ConcatenatedIteratorData<ViewType,DerivedClass>::
operator->() const -> const value_type *
{
    return &(this->operator*());
}


template <class ViewType,class ConstViewType>
inline
ConcatenatedConstIterator<ViewType,ConstViewType>::
ConcatenatedConstIterator(
    const SafeSTLVector<ConstViewType> &ranges,
    const Index index)
    :
    ConcatenatedIteratorData<
    ConstViewType,
    ConcatenatedConstIterator<ViewType,ConstViewType>
    >(ranges,index)
{}

template <class ViewType,class ConstViewType>
inline
ConcatenatedConstIterator<ViewType,ConstViewType>::
ConcatenatedConstIterator(const ConcatenatedIterator<ViewType> &it)
{
    for (const auto &rng : it.get_ranges())
        this->ranges_.emplace_back(ConstViewType(rng)) ;

    this->range_id_ = it.get_range_id();
    this->iterator_current_ = it.get_iterator_current();
}


template <class ViewType,class ConstViewType>
inline
auto
ConcatenatedConstIterator<ViewType,ConstViewType>::
operator[](const Index id) const -> const typename ConstViewType::reference
{
    return this->get_entry_const_reference(id);
}



template <class ViewType>
inline
ConcatenatedIterator<ViewType>::
ConcatenatedIterator(
    const SafeSTLVector<ViewType> &ranges,
    const Index index)
    :
    ConcatenatedIteratorData<
    ViewType,
    ConcatenatedIterator<ViewType>
    >(ranges,index)
{}


template <class ViewType>
inline
auto
ConcatenatedIterator<ViewType>::
operator*() -> value_type &
{
    Assert(this->iterator_current_ != this->ranges_.back().end(),ExcIteratorPastEnd());
    return *this->iterator_current_;
}

template <class ViewType>
inline
auto
ConcatenatedIterator<ViewType>::
operator->() -> value_type *
{
    return &(this->operator*());
}


template <class ViewType>
inline
auto
ConcatenatedIterator<ViewType>::
operator[](const Index id) -> typename ViewType::reference
{
    return this->get_entry_reference(id);
}

template <class ViewType>
inline
auto
ConcatenatedIterator<ViewType>::
operator[](const Index id) const -> const typename ViewType::reference
{
    return this->get_entry_const_reference(id);
}


IGA_NAMESPACE_CLOSE


#endif // #ifndef CONCATENATED_ITERATOR_INLINE_H_
