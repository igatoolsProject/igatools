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

#ifndef GRID_FORWARD_ITERATOR_H_
#define GRID_FORWARD_ITERATOR_H_

#include <igatools/base/config.h>

#include <iterator>

IGA_NAMESPACE_OPEN


/**
 * @brief Forward iterator on objects that have a "grid-like" structure.
 *
 * It fulfills the requirements of a
 * forward iterator as intended by the Standard Template Library.
 * See the C++ documentation for further
 * details of iterator specification and usage.
 *
 * The object pointed to the GridForwardIterator is called <em>accessor</em> and
 * its type is passed as template argument <tt>Accessor</tt> of the GridForwardIterator.
 *
 * The accessor is an object that can fetch and use data stored in objects that have
 * a "grid-like" structure. The type of the object with this "grid-like" structure,
 * associated to the type <tt>Accessor</tt> can be retrieved with the type
 * <tt>Accessor::AccessorOfType</tt>.
 *
 * Using the accessor, the structure of the "grid-like" classes is hidden
 * from the application program.
 *
 * <h3>Purpose</h3>
 *
 * Iterators are used whenever a loop over all (some) elements
 * is to be performed. These loops can then be coded like this:
 * @code
   auto elem = grid.begin();
   const auto elem_end = grid.end();
   for (; elem!=end; ++elem)
     if (elem->at_boundary())
        elem->vertex(k);
  @endcode
 * Note the usage of <tt>++i</tt> instead of <tt>i++</tt> since this
 * does not involve temporaries and copying. It is recommended to use
 * a fixed value <tt>end</tt> inside the loop instead of
 * <tt>patch.end()</tt>, since the creation and copying of these
 * iterators is rather expensive compared to normal pointers.
 *
 * The same previous loop can be performed or using the C++11 syntax called
 * <em>range-based for loop</em>
 * @code
   for (const auto & elem : grid) // the elem type is: const Accessor&
     if (elem.at_boundary())
        elem.vertex(k);
  @endcode
 *
 *
 *
 *
 *
 * Iterators are not much slower than operating directly on the data
 * structures, since they perform the loops that you had to handcode
 * yourself anyway. Most iterator and accessor functions are inlined.
 *
 * The main functionality of iterators, resides in the <tt>++</tt> operator.
 * This move the iterator forward
 * just as if it were a pointer into an array.
 * Here, this operation is
 * not so easy, since it may include dealing with a tensor product
 * structure or something else. But this is completely
 * hidden from the user, though you can still create an iterator
 * pointing to an arbitrary element.  Actually, the operation of
 * moving iterator forward is not done in the iterator
 * classes, but rather in the <tt>Accessor</tt> classes.
 * Since these are passed
 * as template arguments, you can write your own versions here to add
 * more functionality.
 *
 * Furthermore, the iterators described here satisfy the requirement of
 * input and forward iterators as stated by the C++ standard and
 * the STL documentation. It is therefore possible to use the
 * functions from the algorithm section of the C++ standard,
 * e.g. <em>count_if</em>.
 *
 * <h3>Implementation</h3>
 *
 * The iterator class itself does not have much functionality. It only
 * becomes useful when assigned an <tt>Accessor</tt> (the template
 * parameter), which really does the access to data. An <tt>Accessor</tt>
 * has to fulfill some requirements:
 *
 * <ul>
 * <li> It must have a type called <tt>AccessorOfType</tt> representing
 * the type of the "grid-like" container from which the accessor is getting the data.
 *
 * <li> It must have a member named
 * <tt>present_index</tt> storing the address of the element in the
 * <tt>AccessorOfType</tt> object presently pointed to. These data have to be
 * accessible by all grid iterator.
 *
 * <li> It must have a constructor which takes a <tt>const AccessorOfType*</tt> argument
 * representing the accessed container and
 * an <tt>Index</tt> argument, denoting the index within the "grid-like container".
 *
 * <li> It must have void operator <tt>++</tt> that implements the advance of the accessor
 * whitin the container.
 * </ul>
 * Then the iterator is able to do what it is supposed to. All of the
 * necessary functions are implemented in the <tt>Accessor</tt> base
 * class, but you may write your own version (non-virtual, since we
 * use templates) to add functionality.
 *
 *
 *
 * <h3>Past-the-end iterators</h3>
 *
 * There is a representation of past-the-end-pointers, denoted by special
 * values of the member variable @p present_index:
 * <ul>
 * <li> <tt>present_index != IteratorState::past_the_end</tt>, then the object is valid
 * <li> <tt>present_index == IteratorState::past_the_end</tt>, then the iterator points
 * past the end; in all other cases, the iterator is considered invalid.
 * </ul>
 *
 * @see IteratorState
 *
 * @tparam Accessor Type of the accessor.
 *
 * @author M.Martinelli, S.Pauletti
 * @date 2012,2013,2014
 */
template <typename Accessor>
class GridForwardIterator
    : public std::iterator<std::forward_iterator_tag, Accessor>
{
public:
    /** Type of the accessor. */
    using AccessorType = Accessor;

    /** @name Constructors & destructor */
    ///@{
    /**
     * Default constructor. Not allowed to be used.
     */
    GridForwardIterator() = delete;


    /**
     * Construct a forward iterator on a cartesian grid container
     * like type pointed to by patch and the index of the
     * object pointed to by the iterator.
     */
    GridForwardIterator(typename Accessor::AccessorOfType &patch,
                        const Index index);


    /** Copy constructor. */
    GridForwardIterator(const GridForwardIterator<Accessor> &it) = default;


    /** Move constructor. */
    GridForwardIterator(GridForwardIterator<Accessor> &&it) = default;


    /** Destructor */
    ~GridForwardIterator() = default ;
    ///@}


    /** @name Assignment operators */
    ///@{
    /** Copy assignment operator. */
    GridForwardIterator<Accessor> &
    operator=(const GridForwardIterator<Accessor> &) = default;


    /** Move assignment operator. */
    GridForwardIterator<Accessor> &
    operator=(GridForwardIterator<Accessor> &&) = default;
    ///@}

    /** @name Dereferencing operators */
    ///@{
    /**
     *  Dereferencing operator, returns a
     *  const reference to the accessor.
     */
    const Accessor &operator*() const;


    /**
     *  Dereferencing operator, returns a
     *  reference to the accessor.
     */
    Accessor &operator*();


    /**
     *  Dereferencing operator, returns a
     *  const pointer to the accessor.
     */
    const Accessor *operator->() const;


    /**
     *  Dereferencing operator, returns a
     *  pointer to the accessor.
     */
    Accessor *operator->();
    ///@}


    /** @name Comparison operators */
    ///@{
    /** Compare for equality.*/
    bool operator== (const GridForwardIterator &) const;

    /** Compare for inequality.*/
    bool operator!= (const GridForwardIterator &) const;
    ///@}

    /** @name Advance operator */
    ///@{
    /**
     *  Prefix <tt>++</tt> operator: <tt>++i</tt>. This
     *  operator advances the iterator to
     *  the next element and returns
     *  a reference to <tt>*this</tt>.
     */
    GridForwardIterator<Accessor> &operator++ ();
    ///@}


protected:
    /** Object holding the Real data. */
    Accessor accessor_ ;

};



IGA_NAMESPACE_CLOSE


#endif /* PATCH_ITERATORS_H_ */
