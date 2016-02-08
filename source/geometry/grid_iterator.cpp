//-+--------------------------------------------------------------------
// Igatools a general purpose Isogeometric analysis library.
// Copyright (C) 2012-2016  by the igatools authors (see authors.txt).
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

#include <igatools/geometry/grid.h>
#include <igatools/geometry/grid_iterator.h>

using std::shared_ptr;

IGA_NAMESPACE_OPEN
#if 0
template <class Element>
GridIterator<Element>::
GridIterator(std::shared_ptr<ContainerType> container,
             const ListIt &index,
             const PropId &prop)
  :
  elem_(container->create_element(index, prop))
{}

template <class Element>
GridIterator<Element>::
GridIterator(std::unique_ptr<Element> &&elem)
  :
  elem_(std::move(elem))
{
  Assert(elem_ != nullptr,ExcNullPtr());
}



template <class Element>
GridIterator<Element> &
GridIterator<Element>::
operator++()
{
  elem_->operator++();
  return *this;
}



template <class Element>
bool
GridIterator<Element>::
operator==(const GridIterator<Element> &i) const
{
  return *elem_ == *i.elem_;
}





template <class Element>
bool
GridIterator<Element>::
operator!=(const GridIterator<Element> &i) const
{
  return elem_->operator != (*(i.elem_));
}



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
#endif
IGA_NAMESPACE_CLOSE
