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



#ifndef CONCATENATED_FORWARD_ITERATOR_H_
#define CONCATENATED_FORWARD_ITERATOR_H_

#include <igatools/base/config.h>
#include <igatools/base/logstream.h>

#include <vector>


IGA_NAMESPACE_OPEN


/**
 * This class represents a forward iterator made by the concatenation of several forward iterator.
 *
 * @author M. Martinelli
 * @date 03 June 2014
 */
template <class Iterator>
class ConcatenatedForwardIterator
    : public std::iterator<std::forward_iterator_tag, typename Iterator::value_type>
{
public:
    using value_type = typename Iterator::value_type;
    using reference = typename Iterator::reference;


    /** @name Constructors & destructor */
    ///@{
    /**
     * Default constructor. It does nothing.
     */
    ConcatenatedForwardIterator();

    /**
     * Constructor.
     */
    ConcatenatedForwardIterator(
        const std::vector<std::pair<Iterator,Iterator>> &ranges,
        const Index index);



    /** Copy constructor. */
    ConcatenatedForwardIterator(const ConcatenatedForwardIterator<Iterator> &it) = default;

    /** Move constructor. */
    ConcatenatedForwardIterator(ConcatenatedForwardIterator<Iterator> &&it) = default;

    /** Destructor */
    ~ConcatenatedForwardIterator() = default ;
    ///@}


    /** @name Assignment operators */
    ///@{
    /** Copy assignment operator. */
    ConcatenatedForwardIterator<Iterator> &operator=(
        const ConcatenatedForwardIterator<Iterator> &it) = default;

    /** Move assignment operator. */
    ConcatenatedForwardIterator<Iterator> &operator=(
        ConcatenatedForwardIterator<Iterator> &&it) = default;

    ///@}


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
    bool operator==(const ConcatenatedForwardIterator<Iterator> &it) const;

    /** Compare for inequality.*/
    bool operator!=(const ConcatenatedForwardIterator<Iterator> &it) const;
    ///@}

    /** @name Advance operator */
    ///@{
    /**
     *  Prefix <tt>++</tt> operator: <tt>++i</tt>. This
     *  operator advances the iterator to
     *  the next element and returns
     *  a reference to <tt>*this</tt>.
     */
    ConcatenatedForwardIterator<Iterator> &operator++();
    ///@}



    std::vector<std::pair<Iterator,Iterator>> get_ranges() const;

    /** Prints some information. Mostly used for debug and testing. */
    void print_info(LogStream &out) const;

private:

    std::vector<std::pair<Iterator,Iterator>> ranges_;

    int range_id_ = 0;

    Iterator iterator_current_;
};



IGA_NAMESPACE_CLOSE


#include <igatools/utils/concatenated_forward_iterator-inline.h>

#endif // #ifndef CONCATENATED_FORWARD_ITERATOR_H_
