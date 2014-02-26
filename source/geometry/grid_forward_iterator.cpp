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

#include <igatools/geometry/grid_forward_iterator.h>

IGA_NAMESPACE_OPEN

template <typename Accessor>
inline
GridForwardIterator<Accessor>::
GridForwardIterator(
    typename Accessor::AccessorOfType &patch,
    const Index index)
    :
    accessor_(patch, index)
{}




template <typename Accessor>
inline
GridForwardIterator<Accessor> &GridForwardIterator<Accessor>::operator ++ ()
{
    //TODO: is an assert needed here? If so, implement.
    //Assert ();
    ++accessor_;
    return *this;
}



template <typename Accessor>
inline
const Accessor &
GridForwardIterator<Accessor>::operator * () const
{
    return accessor_;
}



template <typename Accessor>
inline
Accessor &
GridForwardIterator<Accessor>::operator * ()
{
    return accessor_;
}



template <typename Accessor>
inline
const Accessor *
GridForwardIterator<Accessor>::operator -> () const
{
    return &(this->operator* ());
}



template <typename Accessor>
inline
Accessor *
GridForwardIterator<Accessor>::operator -> ()
{
    return &(this->operator* ());
}



template <typename Accessor>
inline
bool
GridForwardIterator<Accessor>::operator == (const GridForwardIterator<Accessor> &i) const
{
    return accessor_ == i.accessor_;
}



template <typename Accessor>
inline
bool
GridForwardIterator<Accessor>::operator != (const GridForwardIterator<Accessor> &i) const
{
    return accessor_.operator != (i.accessor_);
}

IGA_NAMESPACE_CLOSE

#include <igatools/geometry/grid_forward_iterator.inst>
