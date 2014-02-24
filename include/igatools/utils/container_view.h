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

//TODO document this file


/**
 * @author M.Martinelli
 * @date 2014
 * @todo document this class
 */
template <class Container>
class ContainerView
{
public:
    using iterator = typename Container::iterator;
    using const_iterator = typename Container::const_iterator;
    using reference = typename iterator::reference;
    using const_reference = typename const_iterator::reference;

public:
    ContainerView(const iterator begin, const iterator end);

    iterator begin();

    const_iterator begin() const;

    iterator end();

    const_iterator end() const;

    reference operator[](const Index n);

    const_reference operator[](const Index n) const;

private:
    iterator begin_;
    iterator end_;
};


/**
 * @author M.Martinelli
 * @date 2014
 * @todo document this class
 */
template <class Container>
class ConstContainerView
{
public:
    using const_iterator = typename Container::const_iterator;
    using const_reference = typename const_iterator::reference;

public:
    ConstContainerView(const const_iterator begin, const const_iterator end);

    const_iterator begin() const;

    const_iterator end() const;

    const_reference operator[](const Index n) const;

private:
    const_iterator begin_;
    const_iterator end_;
};


IGA_NAMESPACE_CLOSE

#endif //#ifndef CONTAINER_VIEW_H_



#include <igatools/utils/container_view-inline.h>



