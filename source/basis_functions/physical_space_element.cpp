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

#include <igatools/basis_functions/physical_space_element.h>
#include <igatools/base/exceptions.h>


using std::shared_ptr;

IGA_NAMESPACE_OPEN

template<int dim_,int range_,int rank_,int codim_>
PhysicalSpaceElement<dim_,range_,rank_,codim_>::
PhysicalSpaceElement(const std::shared_ptr<ContainerType> phys_space,
                     const Index index)
    :
    parent_t(phys_space,index),
    ref_space_element_(phys_space->get_reference_space()->create_element(index)),
    push_fwd_element_(shared_ptr<PfElemAccessor>(new PfElemAccessor(phys_space->get_map_func(), index)))
{
    Assert(ref_space_element_ != nullptr, ExcNullPtr());
    Assert(push_fwd_element_ != nullptr, ExcNullPtr());
}



template<int dim_,int range_,int rank_,int codim_>
PhysicalSpaceElement<dim_,range_,rank_,codim_>::
PhysicalSpaceElement(const PhysicalSpaceElement<dim_,range_,rank_,codim_> &in,
                     const CopyPolicy &copy_policy)
    :
    parent_t(in,copy_policy)
{
    if (copy_policy == CopyPolicy::shallow)
    {
        ref_space_element_ = in.ref_space_element_;
        push_fwd_element_ = in.push_fwd_element_;
    }
    else
    {
        ref_space_element_ =
            shared_ptr<RefElemAccessor>(new RefElemAccessor(*in.ref_space_element_));
        push_fwd_element_ =
            shared_ptr<PfElemAccessor>(new PfElemAccessor(*in.push_fwd_element_));
    }

    Assert(false,ExcNotTested());
}


template<int dim_,int range_,int rank_,int codim_>
auto
PhysicalSpaceElement<dim_,range_,rank_,codim_>::
clone() const -> std::shared_ptr<self_t>
{
//    auto elem = std::shared_ptr<self_t>(
//                    new self_t(*this,CopyPolicy::deep));
    auto elem = std::make_shared<self_t>(*this,CopyPolicy::deep);
    Assert(elem != nullptr, ExcNullPtr());
    return elem;
}


template<int dim_,int range_,int rank_,int codim_>
void
PhysicalSpaceElement<dim_,range_,rank_,codim_>::
copy_from(const PhysicalSpaceElement<dim_,range_,rank_,codim_> &element,
          const CopyPolicy &copy_policy)
{
    Assert(false,ExcNotImplemented());
//    SpaceElementAccessor<PhysSpace>::copy_from(element,copy_policy);
//
//    PhysSpace::PushForwardType::ElementAccessor::copy_from(element,copy_policy);
//
//    if (copy_policy == CopyPolicy::deep)
//        ref_space_element_->deep_copy_from(element.ref_space_element_);
//    else if (copy_policy == CopyPolicy::shallow)
//        ref_space_element_->deep_copy_from(element.ref_space_element_);
//    else
//    {
//        Assert(false,ExcNotImplemented());
//    }
}

template<int dim_,int range_,int rank_,int codim_>
void
PhysicalSpaceElement<dim_,range_,rank_,codim_>::
deep_copy_from(const PhysicalSpaceElement<dim_,range_,rank_,codim_> &element)
{
    Assert(false,ExcNotImplemented());
    //this->copy_from(element,CopyPolicy::deep);
}


template<int dim_,int range_,int rank_,int codim_>
void
PhysicalSpaceElement<dim_,range_,rank_,codim_>::
shallow_copy_from(const PhysicalSpaceElement<dim_,range_,rank_,codim_> &element)
{
    Assert(false,ExcNotImplemented());
//    this->copy_from(element,CopyPolicy::shallow);
}


template<int dim_,int range_,int rank_,int codim_>
template <int k>
auto
PhysicalSpaceElement<dim_,range_,rank_,codim_>::
get_points(const int j) const -> ValueVector<PhysPoint>
{
    return push_fwd_element_->template get_values<_Value,k>(j);
}

template<int dim_,int range_,int rank_,int codim_>
auto
PhysicalSpaceElement<dim_,range_,rank_,codim_>::
get_element_points() const -> ValueVector<PhysPoint>
{
    return this->template get_points<dim>(0);
}


template<int dim_,int range_,int rank_,int codim_>
auto
PhysicalSpaceElement<dim_,range_,rank_,codim_>::
get_element_w_measures() const -> ValueVector<Real>
{
    return this->template get_w_measures<dim>(0);
}

template<int dim_,int range_,int rank_,int codim_>
Index
PhysicalSpaceElement<dim_,range_,rank_,codim_>::
get_flat_index() const
{
    return parent_t::get_flat_index();
}

template<int dim_,int range_,int rank_,int codim_>
auto
PhysicalSpaceElement<dim_,range_,rank_,codim_>::
get_tensor_index() const -> TensorIndex<dim>
{
    return parent_t::get_tensor_index();
}


template<int dim_,int range_,int rank_,int codim_>
void
PhysicalSpaceElement<dim_,range_,rank_,codim_>::
move_to(const Index flat_index)
{
    this->as_cartesian_grid_element_accessor().move_to(flat_index);
    ref_space_element_->move_to(flat_index);
    push_fwd_element_->move_to(flat_index);
}


template<int dim_,int range_,int rank_,int codim_>
int
PhysicalSpaceElement<dim_,range_,rank_,codim_>::
get_num_basis() const
{
    return ref_space_element_->get_num_basis();
}

template<int dim_,int range_,int rank_,int codim_>
auto
PhysicalSpaceElement<dim_,range_,rank_,codim_>::
get_ref_space_element() const -> const RefElemAccessor &
{
    return *ref_space_element_;
}

template<int dim_,int range_,int rank_,int codim_>
auto
PhysicalSpaceElement<dim_,range_,rank_,codim_>::
get_ref_space_element() -> RefElemAccessor &
{
    return *ref_space_element_;
}

template<int dim_,int range_,int rank_,int codim_>
auto
PhysicalSpaceElement<dim_,range_,rank_,codim_>::
get_grid() const -> const std::shared_ptr<const CartesianGrid<dim> >
{
    return this->get_ref_space_element().get_grid();
}

template<int dim_,int range_,int rank_,int codim_>
auto
PhysicalSpaceElement<dim_,range_,rank_,codim_>::
get_push_forward_accessor() const -> const PfElemAccessor &
{
    return *push_fwd_element_;
}

template<int dim_,int range_,int rank_,int codim_>
auto
PhysicalSpaceElement<dim_,range_,rank_,codim_>::
get_push_forward_accessor() -> PfElemAccessor &
{
    return *push_fwd_element_;
}



template<int dim_,int range_,int rank_,int codim_>
void
PhysicalSpaceElement<dim_,range_,rank_,codim_>::
print_info(LogStream &out) const
{
    out.begin_item("Reference space:");
    ref_space_element_->print_info(out);
    out.end_item();

    out.begin_item("Pushforward:");
    push_fwd_element_->print_info(out);
    out.end_item();
}

template<int dim_,int range_,int rank_,int codim_>
void
PhysicalSpaceElement<dim_,range_,rank_,codim_>::
print_cache_info(LogStream &out) const
{
    out.begin_item("Reference space:");
    ref_space_element_->print_cache_info(out);
    out.end_item();

    out.begin_item("Pushforward:");
    push_fwd_element_->print_cache_info(out);
    out.end_item();
}

IGA_NAMESPACE_CLOSE

//#include <igatools/basis_functions/physical_space_element.inst>
