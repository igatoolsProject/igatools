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

using std::shared_ptr;

IGA_NAMESPACE_OPEN

template <class Accessor, class Allocator>
CartesianGridIteratorBase<Accessor,Allocator>::
CartesianGridIteratorBase(std::shared_ptr<ContainerType> grid,
                          const Index index)
    :
    accessor_(new Accessor(grid, index))
{}



template <class Accessor, class Allocator>
CartesianGridIteratorBase<Accessor,Allocator>::
CartesianGridIteratorBase(std::shared_ptr<ContainerType> grid,
                          const TensorIndex<dim> &index)
    :
    accessor_(new Accessor(grid, index))
{}

template <class Accessor, class Allocator>
CartesianGridIteratorBase<Accessor,Allocator>::
CartesianGridIteratorBase(const CartesianGridIteratorBase<Accessor,Allocator> &it,const CopyPolicy &copy_policy)
{
    if (copy_policy == CopyPolicy::deep)
    {
        accessor_ = shared_ptr<Accessor>(new Accessor(*it.accessor_));
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

#if 0
template <class Accessor, class Allocator>
CartesianGridIteratorBase<Accessor,Allocator>::
CartesianGridIteratorBase(const Accessor &acc,const CopyPolicy &copy_policy)
    :
    accessor_(acc,copy_policy)
{}
#endif


template <class Accessor, class Allocator>
bool
CartesianGridIteratorBase<Accessor,Allocator>::
jump(const TensorIndex<dim> &increment)
{
    return accessor_->jump(increment);
}

template <class Accessor, class Allocator>
void
CartesianGridIteratorBase<Accessor,Allocator>::
move_to(const Index flat_index)
{
    accessor_->move_to(flat_index);
}


template <class Accessor, class Allocator>
void
CartesianGridIteratorBase<Accessor,Allocator>::
move_to(const TensorIndex<dim> &tensor_index)
{
    accessor_->move_to(tensor_index);
}


template <class Accessor, class Allocator>
CartesianGridIteratorBase<Accessor,Allocator> &
CartesianGridIteratorBase<Accessor,Allocator>::
operator++()
{
    ++(*accessor_);
    return *this;
}






template <class Accessor, class Allocator>
bool
CartesianGridIteratorBase<Accessor,Allocator>::
operator==(const CartesianGridIteratorBase<Accessor,Allocator> &i) const
{
    return *accessor_ == *i.accessor_;
}


template <class Accessor, class Allocator>
bool
CartesianGridIteratorBase<Accessor,Allocator>::
operator>(const CartesianGridIteratorBase<Accessor,Allocator> &i) const
{
    return (*accessor_ > *i.accessor_);
}

template <class Accessor, class Allocator>
bool
CartesianGridIteratorBase<Accessor,Allocator>::
operator<(const CartesianGridIteratorBase<Accessor,Allocator> &i) const
{
    return (*accessor_ < *i.accessor_);
}



template <class Accessor, class Allocator>
bool
CartesianGridIteratorBase<Accessor,Allocator>::
operator!=(const CartesianGridIteratorBase<Accessor,Allocator> &i) const
{
    return accessor_->operator != (*(i.accessor_));
}

template <class Accessor, class Allocator>
Index
CartesianGridIteratorBase<Accessor,Allocator>::
get_flat_index() const
{
    return accessor_->get_flat_index();
}

template <class Accessor, class Allocator>
auto
CartesianGridIteratorBase<Accessor,Allocator>::
get_tensor_index() const -> TensorIndex<dim>
{
    return accessor_->get_tensor_index();
}




template <class Accessor, class Allocator>
CartesianGridIterator<Accessor,Allocator>::
CartesianGridIterator(std::shared_ptr<ContainerType> grid,const Index index)
    :
    CartesianGridIteratorBase<Accessor,Allocator>(grid,index)
{}

template <class Accessor, class Allocator>
CartesianGridIterator<Accessor,Allocator>::
CartesianGridIterator(std::shared_ptr<ContainerType> grid,
                      const TensorIndex<ContainerType::dim> &index)
    :
    CartesianGridIteratorBase<Accessor,Allocator>(grid,index)
{}


template <class Accessor, class Allocator>
CartesianGridIterator<Accessor,Allocator>::
CartesianGridIterator(const CartesianGridIterator<Accessor,Allocator> &it,const CopyPolicy &copy_policy)
    :
    CartesianGridIteratorBase<Accessor,Allocator>(it,copy_policy)
{}


template <class Accessor, class Allocator>
Accessor &
CartesianGridIterator<Accessor,Allocator>::
operator * ()
{
    return *this->accessor_;
}


template <class Accessor, class Allocator>
Accessor *
CartesianGridIterator<Accessor,Allocator>::
operator -> ()
{
    return this->accessor_.get();
}

template <class Accessor, class Allocator>
const Accessor &
CartesianGridIterator<Accessor,Allocator>::
operator * () const
{
    return *this->accessor_;
}


template <class Accessor, class Allocator>
const Accessor *
CartesianGridIterator<Accessor,Allocator>::
operator -> () const
{
    return this->accessor_.get();
}


template <class Accessor, class Allocator>
CartesianGridConstIterator<Accessor,Allocator>::
CartesianGridConstIterator(std::shared_ptr<ContainerType> grid,const Index index)
    :
    CartesianGridIteratorBase<Accessor,Allocator>(grid,index)
{}

template <class Accessor, class Allocator>
CartesianGridConstIterator<Accessor,Allocator>::
CartesianGridConstIterator(std::shared_ptr<ContainerType> grid,
                           const TensorIndex<ContainerType::dim> &index)
    :
    CartesianGridIteratorBase<Accessor,Allocator>(grid,index)
{}


template <class Accessor, class Allocator>
CartesianGridConstIterator<Accessor,Allocator>::
CartesianGridConstIterator(const CartesianGridConstIterator<Accessor,Allocator> &it,const CopyPolicy &copy_policy)
    :
    CartesianGridIteratorBase<Accessor,Allocator>(it,copy_policy)
{}


template <class Accessor, class Allocator>
const Accessor &
CartesianGridConstIterator<Accessor,Allocator>::
operator * () const
{
    return *this->accessor_;
}



template <class Accessor, class Allocator>
const Accessor *
CartesianGridConstIterator<Accessor,Allocator>::
operator -> () const
{
    return this->accessor_.get();
}




IGA_NAMESPACE_CLOSE
