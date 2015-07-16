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



#ifndef CONTAINER_VIEW_H_
#define CONTAINER_VIEW_H_

#include <igatools/base/config.h>
#include <igatools/base/exceptions.h>
#include <igatools/utils/multi_array_iterator.h>



IGA_NAMESPACE_OPEN

/**
 *
 * @todo Document this class
 */
template <class IteratorType>
class ViewData
{
public:
    /** Returns the number of entries represented by the ViewData. */
    Size get_num_entries() const;

protected:
    /** @name Assignement operators.*/
    ///@{
    /**
     * Default constructor. It does nothing.
     */
    ViewData() = default;

    /**
     * Builds a ViewData object from two IteratorType objects: one pointing to the beginning
     * of the view, and the other pointing to one-pass-end.
     * @note In Debug mode is checked if the relation <tt>begin < end</tt> is satisfied.
     * If it is not, an assertion will be raised.
     */
    explicit ViewData(const IteratorType begin, const IteratorType end);

    /**
     * Copy constructor.
     */
    ViewData(const ViewData<IteratorType> &view_data) = default;

    /**
     * Move constructor.
     */
    ViewData(ViewData<IteratorType> &&view_data) = default;
    ///@}



    /** @name Assignement operators.*/
    ///@{
    /**
     * Copy assignment operator.
     */
    ViewData<IteratorType> &operator=(const ViewData<IteratorType> &view_data);

    /**
     * Move assignment operator.
     */
    ViewData<IteratorType> &operator=(ViewData<IteratorType> &&view_data) = default;
    ///@}
protected:

    /** Iterator pointing to the first element in the view. */
    IteratorType begin_;

    /** Iterator pointing to one-past-end element in the view. */
    IteratorType end_;
};


template <class Iterator, class ConstIterator>
class NonConstView : public ViewData<Iterator>
{
public:
    /** Type of the iterator. */
    using iterator = Iterator;

    /** Type of the const iterator. */
    using const_iterator = ConstIterator;

    /** Type of the reference. */
    using reference = typename iterator::reference;


    /** @name Constructor & destructor */
    ///@{

    /**
     * Default constructor. It does nothing.
     */
    NonConstView() = default;

    /**
     * Construct a view defined by the iterator @p begin pointing to the first element,
     * and the iterator @p end pointing to one-past-the-end element
     * satisfying the chosen criteria.
     */
    explicit NonConstView(const iterator begin, const iterator end);


    /** Copy constructor. */
    NonConstView(const NonConstView<Iterator,ConstIterator> &view) = default;

    /** Move constructor. */
    NonConstView(NonConstView<Iterator,ConstIterator> &&view) = default;

    /** Destructor. */
    ~NonConstView() = default;
    ///@}

    /** @name Assignment operators */
    ///@{

    /** Copy assignment operator. */
    NonConstView<Iterator,ConstIterator> &operator=(const NonConstView<Iterator,ConstIterator> &view) = default;

    /** Move assignment operator. */
    NonConstView<Iterator,ConstIterator> &operator=(NonConstView<Iterator,ConstIterator> &&view) = default;
    ///@}

    /** @name Dealing with the iterator */
    ///@{
    /** Returns an iterator pointing to the first element in the view. */
    iterator begin();

    /** Returns a const-iterator pointing to the first element in the view. */
    const_iterator begin() const;

    /** Returns a const-iterator pointing to the first element in the view. */
    const_iterator cbegin() const;

    /** Returns an iterator pointing to one-pass-the end element in the view. */
    iterator end();

    /** Returns a const-iterator pointing to one-pass-the end element in the view. */
    const_iterator end() const;

    /** Returns a const-iterator pointing to one-pass-the end element in the view. */
    const_iterator cend() const;
    ///@}

    /** @name Dereference offset operators */
    ///@{
    /** Return a reference to the <tt>n</tt>-th element in the view. */
    reference operator[](const Index n);

    /** Return a const reference to the <tt>n</tt>-th element in the view. */
    const reference operator[](const Index n) const;
    ///@}
};


template <class Iterator, class ConstIterator>
class ConstView : public ViewData<ConstIterator>
{
public:
    /** Type of the const iterator. */
    using iterator = ConstIterator;

    /** Type of the const iterator. */
    using const_iterator = ConstIterator;

    /** Type of the reference. */
    using reference = typename iterator::reference;

    /** @name Constructor & destructor */
    ///@{

    /**
     * Default constructor. It does nothing.
     */
    ConstView() = default;

    /**
     * Construct a view defined by the const iterator @p begin pointing to the first element,
     * and the const iterator @p end pointing to one-past-the-end element
     * satisfying the chosen criteria.
     */
    explicit ConstView(const const_iterator begin, const const_iterator end);

    /**
     * Construct a ConstView from a View.
     */
    explicit ConstView(const NonConstView<Iterator,ConstIterator> &view);

    /** Copy constructor. */
    ConstView(const ConstView<Iterator,ConstIterator> &view) = default;

    /** Move constructor. */
    ConstView(ConstView<Iterator,ConstIterator> &&view) = default;

    /** Destructor. */
    ~ConstView() = default;
    ///@}

    /** @name Assignment operators */
    ///@{
    /** Copy assignment operator. */
    ConstView<Iterator,ConstIterator> &operator=(const ConstView<Iterator,ConstIterator> &view) = default;

    /** Move assignment operator. */
    ConstView<Iterator,ConstIterator> &operator=(ConstView<Iterator,ConstIterator> &&view) = default;
    ///@}


    /** @name Dealing with the iterator */
    ///@{
    /** Returns a const-iterator pointing to the first element in the view. */
    const_iterator begin() const;

    /** Returns a const-iterator pointing to one-pass-the end element in the view. */
    const_iterator end() const;


    /** Returns a const-iterator pointing to the first element in the view. */
    const_iterator cbegin() const;

    /** Returns a const-iterator pointing to one-pass-the end element in the view. */
    const_iterator cend() const;
    ///@}

    /** @name Dereference offset operators */
    ///@{
    /** Return a const reference to the <tt>n</tt>-th element in the view. */
    const reference operator[](const Index n) const;
    ///@}
};



/**
 * @brief This class represents a const "view" of a container of type <tt>Container</tt>.
 *
 * For the class documentation see the analogous non-const version ContainerView.
 *
 * @see ContainerView
 * @author M.Martinelli
 * @date 2014
 */
template <class Container>
class ConstContainerView : public ConstView<typename Container::iterator,typename Container::const_iterator>
{
public:
    using ConstView<typename Container::iterator,typename Container::const_iterator>::ConstView;
//    ConstContainerView(const ConstContainerView<Container> &view) = default;
//    ConstContainerView(ConstContainerView<Container> &&view) = default;
};


/**
 * @brief This class represents a "view" of a container of type <tt>Container</tt>.
 *
 * A "view" is a special iterator that operates on entries
 * that can be grouped in accordance with some criteria.
 * In order to do so, the "view" needs to be constructed from two Container::iterator object
 * (one pointing to the first element and the other pointing to one-pass-the-end) that
 * satisfy the chosen criteria.
 *
 * For example a criteria could be that from an element pointed by the view and the next one
 * corresponds to a given number of elements in the container
 * (i.e. the elements are iterated with a <em>stride > 1</em>).
 *
 * @see ValueTable
 *
 * @author M.Martinelli
 * @date 2014
 */
template <class Container>
class ContainerView : public NonConstView<typename Container::iterator,typename Container::const_iterator>
{
public:
    using NonConstView<typename Container::iterator,typename Container::const_iterator>::NonConstView;
//  ContainerView(const ContainerView<Container> &view) = default;
//  ContainerView(ContainerView<Container> &&view) = default;
};


IGA_NAMESPACE_CLOSE

#endif //#ifndef CONTAINER_VIEW_H_



#include <igatools/utils/container_view-inline.h>



