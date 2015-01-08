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

#include <igatools/geometry/cartesian_grid_iterator.h>

using std::shared_ptr;

IGA_NAMESPACE_OPEN


template <class Accessor>
CartesianGridIteratorBase<Accessor>::
CartesianGridIteratorBase(std::shared_ptr<ContainerType> grid,
                          const Index index)
    :
    CartesianGridIteratorBase(grid->create_element(index))
{}


#if 0
template <class Accessor>
CartesianGridIteratorBase<Accessor>::
CartesianGridIteratorBase(std::shared_ptr<ContainerType> grid,
                          const TensorIndex<dim> &index)
    :
    CartesianGridIteratorBase(grid,grid->tensor_to_flat(index))
{}
#endif

template <class Accessor>
CartesianGridIteratorBase<Accessor>::
CartesianGridIteratorBase(std::shared_ptr<Accessor> accessor_ptr)
    :
    accessor_(accessor_ptr)
{
    Assert(accessor_ != nullptr,ExcNullPtr());
}

template <class Accessor>
CartesianGridIteratorBase<Accessor>::
CartesianGridIteratorBase(const CartesianGridIteratorBase<Accessor> &it,const CopyPolicy &copy_policy)
{
    if (copy_policy == CopyPolicy::deep)
    {
        accessor_ = it.accessor_->clone();
    }
    else if (copy_policy == CopyPolicy::shallow)
    {
        accessor_ = it.accessor_;
    }
    else
    {
        Assert(false,ExcNotImplemented());
        AssertThrow(false,ExcNotImplemented());
    }
}


template <class Accessor>
bool
CartesianGridIteratorBase<Accessor>::
jump(const TensorIndex<dim> &increment)
{
    return accessor_->jump(increment);
}

template <class Accessor>
void
CartesianGridIteratorBase<Accessor>::
move_to(const Index flat_index)
{
    accessor_->move_to(flat_index);
}


template <class Accessor>
void
CartesianGridIteratorBase<Accessor>::
move_to(const TensorIndex<dim> &tensor_index)
{
    accessor_->move_to(tensor_index);
}


template <class Accessor>
CartesianGridIteratorBase<Accessor> &
CartesianGridIteratorBase<Accessor>::
operator++()
{
    ++(*accessor_);
    return *this;
}






template <class Accessor>
bool
CartesianGridIteratorBase<Accessor>::
operator==(const CartesianGridIteratorBase<Accessor> &i) const
{
    return *accessor_ == *i.accessor_;
}


template <class Accessor>
bool
CartesianGridIteratorBase<Accessor>::
operator>(const CartesianGridIteratorBase<Accessor> &i) const
{
    return (*accessor_ > *i.accessor_);
}

template <class Accessor>
bool
CartesianGridIteratorBase<Accessor>::
operator<(const CartesianGridIteratorBase<Accessor> &i) const
{
    return (*accessor_ < *i.accessor_);
}



template <class Accessor>
bool
CartesianGridIteratorBase<Accessor>::
operator!=(const CartesianGridIteratorBase<Accessor> &i) const
{
    return accessor_->operator != (*(i.accessor_));
}

template <class Accessor>
Index
CartesianGridIteratorBase<Accessor>::
get_flat_index() const
{
    return accessor_->get_flat_index();
}

template <class Accessor>
auto
CartesianGridIteratorBase<Accessor>::
get_tensor_index() const -> TensorIndex<dim>
{
    return accessor_->get_tensor_index();
}



#if 0
template <class Accessor>
CartesianGridIterator<Accessor>::
CartesianGridIterator(std::shared_ptr<ContainerType> grid,const Index index)
    :
    CartesianGridIteratorBase<Accessor>(grid,index)
{}

template <class Accessor>
CartesianGridIterator<Accessor>::
CartesianGridIterator(std::shared_ptr<ContainerType> grid,
                      const TensorIndex<ContainerType::dim> &index)
    :
    CartesianGridIteratorBase<Accessor>(grid,index)
{}


template <class Accessor>
CartesianGridIterator<Accessor>::
CartesianGridIterator(const CartesianGridIterator<Accessor> &it,const CopyPolicy &copy_policy)
    :
    CartesianGridIteratorBase<Accessor>(it,copy_policy)
{}
#endif

template <class Accessor>
Accessor &
CartesianGridIterator<Accessor>::
operator * ()
{
    return *this->accessor_;
}


template <class Accessor>
Accessor *
CartesianGridIterator<Accessor>::
operator -> ()
{
    return this->accessor_.get();
}

template <class Accessor>
const Accessor &
CartesianGridIterator<Accessor>::
operator * () const
{
    return *this->accessor_;
}


template <class Accessor>
const Accessor *
CartesianGridIterator<Accessor>::
operator -> () const
{
    return this->accessor_.get();
}

#if 0
template <class Accessor>
CartesianGridConstIterator<Accessor>::
CartesianGridConstIterator(std::shared_ptr<ContainerType> grid,const Index index)
    :
    CartesianGridIteratorBase<Accessor>(grid,index)
{}

template <class Accessor>
CartesianGridConstIterator<Accessor>::
CartesianGridConstIterator(std::shared_ptr<ContainerType> grid,
                           const TensorIndex<ContainerType::dim> &index)
    :
    CartesianGridIteratorBase<Accessor>(grid,index)
{}


template <class Accessor>
CartesianGridConstIterator<Accessor>::
CartesianGridConstIterator(const CartesianGridConstIterator<Accessor> &it,const CopyPolicy &copy_policy)
    :
    CartesianGridIteratorBase<Accessor>(it,copy_policy)
{}
#endif

template <class Accessor>
const Accessor &
CartesianGridConstIterator<Accessor>::
operator * () const
{
    return *this->accessor_;
}



template <class Accessor>
const Accessor *
CartesianGridConstIterator<Accessor>::
operator -> () const
{
    return this->accessor_.get();
}




IGA_NAMESPACE_CLOSE
