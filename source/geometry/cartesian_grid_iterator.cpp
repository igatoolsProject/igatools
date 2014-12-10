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

#include <igatools/geometry/cartesian_grid_iterator.h>

IGA_NAMESPACE_OPEN

template <typename Accessor>
CartesianGridIteratorBase<Accessor>::
CartesianGridIteratorBase(std::shared_ptr<ContainerType> grid,
                          const Index index)
    :
    accessor_(grid, index)
{}



template <typename Accessor>
CartesianGridIteratorBase<Accessor>::
CartesianGridIteratorBase(std::shared_ptr<ContainerType> grid,
                          const TensorIndex<dim> &index)
    :
    accessor_(grid, index)
{}

template <typename Accessor>
CartesianGridIteratorBase<Accessor>::
CartesianGridIteratorBase(const CartesianGridIteratorBase<Accessor> &it,const CopyPolicy &copy_policy)
    :
    accessor_(it.accessor_,copy_policy)
{}


template <typename Accessor>
CartesianGridIteratorBase<Accessor>::
CartesianGridIteratorBase(const Accessor &acc,const CopyPolicy &copy_policy)
    :
    accessor_(acc,copy_policy)
{}



template <typename Accessor>
bool
CartesianGridIteratorBase<Accessor>::
jump(const TensorIndex<dim> &increment)
{
    return accessor_.jump(increment);
}

template <typename Accessor>
void
CartesianGridIteratorBase<Accessor>::
move_to(const Index flat_index)
{
    accessor_.move_to(flat_index);
}


template <typename Accessor>
void
CartesianGridIteratorBase<Accessor>::
move_to(const TensorIndex<dim> &tensor_index)
{
    accessor_.move_to(tensor_index);
}


template <typename Accessor>
CartesianGridIteratorBase<Accessor> &
CartesianGridIteratorBase<Accessor>::
operator++()
{
    ++accessor_;
    return *this;
}

#if 0

template <typename Accessor>
const Accessor &
CartesianGridIteratorBase<Accessor>::
operator * () const
{
    return accessor_;
}





template <typename Accessor>
const Accessor *
CartesianGridIteratorBase<Accessor>::
operator -> () const
{
    return &(this->operator* ());
}

#endif





template <typename Accessor>
bool
CartesianGridIteratorBase<Accessor>::
operator==(const CartesianGridIteratorBase<Accessor> &i) const
{
    return accessor_ == i.accessor_;
}


template <typename Accessor>
bool
CartesianGridIteratorBase<Accessor>::
operator>(const CartesianGridIteratorBase<Accessor> &i) const
{
    return (accessor_ > i.accessor_);
}

template <typename Accessor>
bool
CartesianGridIteratorBase<Accessor>::
operator<(const CartesianGridIteratorBase<Accessor> &i) const
{
    return (accessor_ < i.accessor_);
}



template <typename Accessor>
bool
CartesianGridIteratorBase<Accessor>::
operator!=(const CartesianGridIteratorBase<Accessor> &i) const
{
    return accessor_.operator != (i.accessor_);
}

template <typename Accessor>
Index
CartesianGridIteratorBase<Accessor>::
get_flat_index() const
{
    return accessor_.get_flat_index();
}

template <typename Accessor>
auto
CartesianGridIteratorBase<Accessor>::
get_tensor_index() const -> TensorIndex<dim>
{
    return accessor_.get_tensor_index();
}



template <typename Accessor>
Accessor &
CartesianGridIterator<Accessor>::
operator * ()
{
    return this->accessor_;
}


template <typename Accessor>
Accessor *
CartesianGridIterator<Accessor>::
operator -> ()
{
    return &(this->operator* ());
}

template <typename Accessor>
const Accessor &
CartesianGridIterator<Accessor>::
operator * () const
{
    return this->accessor_;
}


template <typename Accessor>
const Accessor *
CartesianGridIterator<Accessor>::
operator -> () const
{
    return &(this->operator* ());
}

IGA_NAMESPACE_CLOSE
