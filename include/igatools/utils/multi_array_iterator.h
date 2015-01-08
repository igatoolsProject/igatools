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

#ifndef MULTI_ARRAY_ITERATOR_H_
#define MULTI_ARRAY_ITERATOR_H_

#include <igatools/base/config.h>
#include <igatools/utils/tensor_size.h>
#include <igatools/utils/tensor_index.h>

#include <iterator>

IGA_NAMESPACE_OPEN

/**
 * @brief Random access iterator for the MultiArray container
 * (and for the derived classes DynamicMultiArray and StaticMultiArray).
 *
 * This iterator can move in a non-contiguous way along the @p Container entries,
 * specifying (using the proper the constructor parameter) a @p stride value greater
 * than @p 1.
 *
 * @ingroup multi_array_containers
 *
 * @author M. Martinelli
 * @date 04 Feb 2014
 */
template <class Container>
class MultiArrayIteratorBase
    : public std::iterator<
    std::random_access_iterator_tag,
    typename Container::value_type >
{
public:
    using parent_t = std::iterator<
                     std::random_access_iterator_tag,
                     typename Container::value_type >;

    /**
     * Type for the values referenced by the iterator.
     * If the @p Container type is declared to be @p const then also
     * the @p Value type will be @p const.
     */
    using value_type = typename Container::value_type;

    /** Type for the reference. */
    using reference = typename Container::reference;


    /** Type for the pointer. */
    using pointer = value_type *;

    using difference_type = typename std::iterator_traits<parent_t>::difference_type;


    /** @name Constructors and destructor */
    ///@{

    /** Default constructor.*/
    MultiArrayIteratorBase();

    /** Copy constructor. */
    MultiArrayIteratorBase(const MultiArrayIteratorBase<Container> &in) = default;

    /** Move constructor.*/
    MultiArrayIteratorBase(MultiArrayIteratorBase<Container> &&in) = default;


    /**
     * Construct a iterator on a tensor sized container
     * like type pointed to by the @p container, the @p index of the
     * object pointed to by the iterator and the @p stride used
     * to advance to the next entry in the container.
     *
     * @note In order to successfully build the iterator, the @p container must not be empty,
     * otherwise, in DEBUG mode, an exception will be raised.
     */
    MultiArrayIteratorBase(Container &container,const Index id,const Index stride=1);

    /** Destructor. */
    ~MultiArrayIteratorBase() = default;
    ///@}


    /** @name Assignment operators */
    ///@{
    /** Copy assignment operator.*/
    MultiArrayIteratorBase<Container> &operator=(
        const MultiArrayIteratorBase<Container> &it) = default;

    /** Move assignment operator.*/
    MultiArrayIteratorBase<Container> &operator=(
        MultiArrayIteratorBase<Container> &&) = default;

    ///@}


    /** @name Advance operators */
    ///@{
    /**
     *  Prefix <tt>++</tt> operator: <tt>++i</tt>. This
     *  operator advances the iterator to
     *  the next entry and returns
     *  a reference to <tt>*this</tt>.
     */
    MultiArrayIteratorBase<Container> &operator++();

    /**
     *  Postfix <tt>++</tt> operator: <tt>i++</tt>.
     */
    MultiArrayIteratorBase<Container> operator++(int);


    /**
     * This operator returns an iterator pointing to the entry obtained by advancing
     * @p n positions from the current one.
     *
     * @note The positions are ``counted'' considering the stride used to build the iterator so,
     * if the original iterator was pointing the position identified by @p id, the returned iterator
     * will point the position identified by <tt>id + n*stride</tt>.
     */
    MultiArrayIteratorBase<Container> operator+(const Index n) const;


    /**
     * This operator returns an iterator pointing to the entry obtained by going back
     * @p n positions from the current one.
     *
     * @note The positions are ``counted'' considering the stride used to build the iterator so,
     * if the original iterator was pointing the position identified by @p id, the returned iterator
     * will point the position identified by <tt>id - n*stride</tt>.
     */
    MultiArrayIteratorBase<Container> operator-(const Index n) const;

    /**
     * Returns <tt>n</tt> such that <tt>a + n == (*this)</tt>, where <tt>(*this) == a + ( (*this) - a)</tt>.
     *
     */
    difference_type
    operator-(const MultiArrayIteratorBase<Container> &a) const;
    ///@}

    /** @name Dereferencing operators */
    ///@{
    /** Dereferencing operator (const version).*/
    const reference operator* () const;

    /** Dereferencing operator (non-const version).*/
    reference operator* ();

    /** Dereferencing operator (const version).*/
    const pointer operator-> () const;

    /** Dereferencing operator (non-const version).*/
    pointer operator->();


    /** Dereferencing operator (const version). Returns the i-th entry in the iterator. */
    const reference operator[](const Index i) const;

    /** Dereferencing operator (non-const version). Returns the i-th entry in the iterator. */
    reference operator[](const Index i);

    ///@}


    /** @name Comparison operators. */
    ///@{
    /** Compare for equality. */
    bool operator==(const MultiArrayIteratorBase<Container> &) const;


    /** Compare for inequality.*/
    bool operator!=(const MultiArrayIteratorBase<Container> &) const;


    /**
     * Returns true if,
     * given two iterators (on the same @p Container object and with the same @p stride),
     * the left one is pointing to previous entry w.r.t. the right one.
     */
    bool operator<(const MultiArrayIteratorBase<Container> &) const;

    /** Compare for less-than or equality. */
    bool operator<=(const MultiArrayIteratorBase<Container> &) const;
    ///@}



    /** Returns the flat index of the entry in the container pointed by the iterator. */
    Index get_id() const;

    /** Returns the stride used to advance to the next entry in the container. */
    Index get_stride() const;

private:

    /** Container upon which the iterator is built for. */
    Container *container_;

    /** Flat index of the entry in the container pointed by the iterator. */
    Index id_;

    /** Stride used to advance to the next entry in the container. */
    Index stride_;

};




template <class Container>
class MultiArrayConstIterator
    : public MultiArrayIteratorBase<Container>
{
public:
    MultiArrayConstIterator() = default;

    /**
     * Construct a iterator on a tensor sized container
     * like type pointed to by the @p container, the @p index of the
     * object pointed to by the iterator and the @p stride used
     * to advance to the next entry in the container.
     *
     * @note In order to successfully build the iterator, the @p container must not be empty,
     * otherwise, in DEBUG mode, an exception will be raised.
     */
    MultiArrayConstIterator(const Container &container,const Index id,const Index stride=1);
};


template <class Container>
class MultiArrayIterator
    : public MultiArrayConstIterator<Container>
{
public:
    MultiArrayIterator() = default;

    /**
     * Construct a iterator on a tensor sized container
     * like type pointed to by the @p container, the @p index of the
     * object pointed to by the iterator and the @p stride used
     * to advance to the next entry in the container.
     *
     * @note In order to successfully build the iterator, the @p container must not be empty,
     * otherwise, in DEBUG mode, an exception will be raised.
     */
    MultiArrayIterator(Container &container,const Index id,const Index stride=1);
};


IGA_NAMESPACE_CLOSE


#include <igatools/utils/multi_array_iterator-inline.h>

#endif // #ifndef MULTI_ARRAY_ITERATOR_H_

