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

#include <igatools/geometry/cartesian_grid.h>
#include <igatools/geometry/grid_iterator.h>

using std::shared_ptr;

IGA_NAMESPACE_OPEN

template <class Element>
GridIteratorBase<Element>::
GridIteratorBase(std::shared_ptr<ContainerType> grid,
                 const ListIt &index,
                 const PropId &prop)
    :
    elem_(std::make_shared<Element>(Element(grid, index, prop)))
{}



template <class Element>
GridIteratorBase<Element>::
GridIteratorBase(const GridIteratorBase<Element> &it,
                 const CopyPolicy &copy_policy)
{
    if (copy_policy == CopyPolicy::deep)
    {
        elem_->deep_copy_from(*(it.elem_));
    }
    else if (copy_policy == CopyPolicy::shallow)
    {
        elem_->shallow_copy_from(*(it.elem_));
    }
    else
    {
        Assert(false,ExcNotImplemented());
        AssertThrow(false,ExcNotImplemented());
    }
}



template <class Element>
GridIteratorBase<Element> &
GridIteratorBase<Element>::
operator++()
{
    elem_->operator++();
    return *this;
}



template <class Element>
bool
GridIteratorBase<Element>::
operator==(const GridIteratorBase<Element> &i) const
{
    return *elem_ == *i.elem_;
}



template <class Element>
bool
GridIteratorBase<Element>::
operator>(const GridIteratorBase<Element> &i) const
{
    return (*elem_ > *i.elem_);
}

template <class Element>
bool
GridIteratorBase<Element>::
operator<(const GridIteratorBase<Element> &i) const
{
    return (*elem_ < *i.elem_);
}



template <class Element>
bool
GridIteratorBase<Element>::
operator!=(const GridIteratorBase<Element> &i) const
{
    return elem_->operator != (*(i.elem_));
}

#if 0
template <class Element>
Index
GridIteratorBase<Element>::
get_flat_index() const
{
    return elem_->get_flat_index();
}

template <class Element>
auto
GridIteratorBase<Element>::
get_tensor_index() const -> TensIndex
{
    return elem_->get_tensor_index();
}
#endif



template <class Element>
Element &
GridIterator<Element>::
operator * ()
{
    return *this->elem_;
}



template <class Element>
Element *
GridIterator<Element>::
operator -> ()
{
    return this->elem_.get();
}



template <class Element>
const Element &
GridIterator<Element>::
operator * () const
{
    return *this->elem_;
}



template <class Element>
const Element *
GridIterator<Element>::
operator -> () const
{
    return this->elem_.get();
}

IGA_NAMESPACE_CLOSE
