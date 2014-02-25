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



#ifndef CONTAINER_VIEW_H_
#define CONTAINER_VIEW_H_

#include <igatools/base/config.h>
#include <igatools/base/exceptions.h>
#include <igatools/utils/multi_array_iterator.h>



IGA_NAMESPACE_OPEN


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
class ContainerView
{
public:
    /** Type of the iterator. */
    using iterator = typename Container::iterator;

    /** Type of the const iterator. */
    using const_iterator = typename Container::const_iterator;

    /** Type of the reference. */
    using reference = typename iterator::reference;

    /** Type of the const reference. */
    using const_reference = typename const_iterator::reference;

    /** @name Constructor & destructor */
    ///@{
    /**
     * Construct a view defined by the iterator @p begin pointing to the first element,
     * and the iterator @p end pointing to one-pass-the-end element
     * satisfying the chosen criteria.
     */
    explicit ContainerView(const iterator begin, const iterator end);


    /** Copy constructor. */
    ContainerView(const ContainerView<Container> &view) = default;

    /** Move constructor. */
    ContainerView(ContainerView<Container> &&view) = default;

    /** Destructor. */
    ~ContainerView() = default;
    ///@}

    /** @name Assignment operators */
    ///@{
    /** Copy assignment operator. */
    ContainerView<Container> &operator=(const ContainerView<Container> &view) = default;

    /** Move assignment operator. */
    ContainerView<Container> &operator=(ContainerView<Container> &&view) = default;
    ///@}

    /** @name Dealing with the iterator */
    ///@{
    /** Returns an iterator pointing to the first element in the view. */
    iterator begin();

    /** Returns a const-iterator pointing to the first element in the view. */
    const_iterator begin() const;

    /** Returns an iterator pointing to one-pass-the end element in the view. */
    iterator end();

    /** Returns a const-iterator pointing to one-pass-the end element in the view. */
    const_iterator end() const;
    ///@}

    /** @name Dereference offset operators */
    ///@{
    /** Return a reference to the <tt>n</tt>-th element in the view. */
    reference operator[](const Index n);

    /** Return a const reference to the <tt>n</tt>-th element in the view. */
    const_reference operator[](const Index n) const;
    ///@}

private:
    /** Iterator pointing to the first element in the view. */
    iterator begin_;

    /** Iterator pointing to one-past-end element in the view. */
    iterator end_;
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
class ConstContainerView
{
public:
    /** Type of the const iterator. */
    using const_iterator = typename Container::const_iterator;

    /** Type of the const reference. */
    using const_reference = typename const_iterator::reference;

    /** @name Constructor & destructor */
    ///@{
    /**
     * Construct a view defined by the iterator @p begin pointing to the first element,
     * and the iterator @p end pointing to one-pass-the-end element
     * satisfying the chosen criteria.
     */
    explicit ConstContainerView(const iterator begin, const iterator end);


    /** Copy constructor. */
    ConstContainerView(const ConstContainerView<Container> &view) = default;

    /** Move constructor. */
    ConstContainerView(ConstContainerView<Container> &&view) = default;

    /** Destructor. */
    ~ConstContainerView() = default;
    ///@}

    /** @name Assignment operators */
    ///@{
    /** Copy assignment operator. */
    ConstContainerView<Container> &operator=(const ConstContainerView<Container> &view) = default;

    /** Move assignment operator. */
    ConstContainerView<Container> &operator=(ConstContainerView<Container> &&view) = default;
    ///@}


    /** @name Dealing with the iterator */
    ///@{
    /** Returns a const-iterator pointing to the first element in the view. */
    const_iterator begin() const;

    /** Returns a const-iterator pointing to one-pass-the end element in the view. */
    const_iterator end() const;
    ///@}

    /** @name Dereference offset operators */
    ///@{
    /** Return a const reference to the <tt>n</tt>-th element in the view. */
    const_reference operator[](const Index n) const;
    ///@}

private:
    /** Iterator pointing to the first element in the view. */
    const_iterator begin_;

    /** Iterator pointing to one-past-end element in the view. */
    const_iterator end_;
};


IGA_NAMESPACE_CLOSE

#endif //#ifndef CONTAINER_VIEW_H_



#include <igatools/utils/container_view-inline.h>



