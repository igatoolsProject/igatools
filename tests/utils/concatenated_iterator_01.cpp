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
/*
 * martinelli
 * 30 May 2014
 *
 */

#include "../tests.h"

#include <igatools/base/exceptions.h>
#include <vector>

using std::vector;

template <class Iterator>
class ConcatenatedIterator
    : public std::iterator<std::forward_iterator_tag, typename Iterator::value_type>
{
public:
    using value_type = typename Iterator::value_type;

    void push_back(const Iterator &begin,const Iterator &end)
    {
        Assert(begin < end,ExcInvalidIterator());
        ranges_.push_back(std::make_pair(begin,end));
    }


    /** @name Dereferencing operators */
    ///@{
    /**
     *  Dereferencing operator, returns a
     *  const reference to the value_type.
     */
    const value_type &operator*() const;

    /**
     *  Dereferencing operator, returns a
     *  reference to the value_type.
     */
    value_type &operator*();

    /**
     *  Dereferencing operator, returns a
     *  const pointer to the value_type.
     */
    const value_type *operator->() const;

    /**
     *  Dereferencing operator, returns a
     *  pointer to the value_type.
     */
    value_type *operator->();
    ///@}

    /** @name Comparison operators */
    ///@{
    /** Compare for equality.*/
    bool operator== (const ConcatenatedIterator<Iterator> &) const;

    /** Compare for inequality.*/
    bool operator!= (const ConcatenatedIterator<Iterator> &) const;
    ///@}

    /** @name Advance operator */
    ///@{
    /**
     *  Prefix <tt>++</tt> operator: <tt>++i</tt>. This
     *  operator advances the iterator to
     *  the next element and returns
     *  a reference to <tt>*this</tt>.
     */
    ConcatenatedIterator<Iterator> &operator++();
    ///@}


    void init()
    {
        Assert(!ranges_.empty(),ExcEmptyObject());

        iterator_current_ = ranges_[0].first;
    }

private:

    std::vector<std::pair<Iterator,Iterator>> ranges_;

    int range_id_ = 0;

    Iterator iterator_current_;
};


template <class Iterator>
auto
ConcatenatedIterator<Iterator>::
operator*() -> value_type &
{
    return *iterator_current_;
}

template <class Iterator>
auto
ConcatenatedIterator<Iterator>::
operator*() const -> const value_type &
{
    return *iterator_current_;
}

template <class Iterator>
auto
ConcatenatedIterator<Iterator>::
operator++() -> ConcatenatedIterator<Iterator> &
{
    if (range_id_ < ranges_.size()-1)
    {
        // if the current iterator is before the end, advance one position
        if (iterator_current_ < ranges_[range_id_].second)
            ++iterator_current_;

        // if the current iterator is already at the end of one iterator,
        // point to the first element of the next one
        if (iterator_current_ == ranges_[range_id_].second)
            iterator_current_ = ranges_[++range_id_].first;
    }
    else
    {
        Assert(iterator_current_ != ranges_[range_id_].second,ExcIteratorPastEnd());
        ++iterator_current_;
    }

    return *this;
}

int main()
{
    vector<int> v0 = {1,2,3,4};
    vector<int> v1 = {5,6,7,8,9};
    vector<int> v2 = {10,11,12};

    using VecIterator = vector<int>::iterator;

    ConcatenatedIterator<VecIterator> concatenated_iterator;
    concatenated_iterator.push_back(v0.begin(),v0.end());
    concatenated_iterator.push_back(v1.begin(),v1.end());
    concatenated_iterator.push_back(v2.begin(),v2.end());


    concatenated_iterator.init();

    using std::cout;
    using std::endl;
    for (int i = 0 ; i < 12 ; ++i)
    {
        cout << "i = " << i << "     value = " << *concatenated_iterator << std::endl;
        ++concatenated_iterator;
    }
}
