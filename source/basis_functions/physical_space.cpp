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

#include <igatools/basis_functions/physical_space.h>
#include <igatools/geometry/mapping_slice.h>
using std::vector;
using std::array;
using std::shared_ptr;
using std::const_pointer_cast;

IGA_NAMESPACE_OPEN

template <class RefSpace_, class PushForward_>
PhysicalSpace<RefSpace_,PushForward_>::
PhysicalSpace(
    shared_ptr<RefSpace> ref_space,
    shared_ptr<PushForwardType> push_forward,
    const Index id)
    :
    BaseSpace(ref_space->get_grid()),
    ref_space_(ref_space),
    push_forward_(push_forward),
    id_(id)
{
//TODO(pauletti, Jan 18, 2014): put static assert on h_div, h_curl range and rank
    Assert(ref_space_ != nullptr, ExcNullPtr());
    Assert(push_forward_ != nullptr, ExcNullPtr());

    Assert(ref_space_->get_grid() == push_forward_->get_mapping()->get_grid(),
           ExcMessage("Reference space and mapping grids are not the same."))

}




template <class RefSpace_, class PushForward_>
auto
PhysicalSpace<RefSpace_,PushForward_>::
clone() const -> shared_ptr<self_t>
{
    Assert(ref_space_ != nullptr, ExcNullPtr());
    Assert(push_forward_ != nullptr, ExcNullPtr());

    return shared_ptr<self_t>(
        new self_t(
            shared_ptr<RefSpace>(new RefSpace(*ref_space_)),
            shared_ptr<PushForwardType>(new PushForwardType(*push_forward_)),
            id_)
    );
};



template <class RefSpace_, class PushForward_>
auto
PhysicalSpace<RefSpace_,PushForward_>::
create(
    shared_ptr<RefSpace> ref_space,
    shared_ptr<PushForwardType> push_forward,
    const Index id) -> shared_ptr<self_t>
{
    Assert(ref_space != nullptr, ExcNullPtr());
    Assert(push_forward != nullptr, ExcNullPtr());
    return shared_ptr<self_t>(new self_t(ref_space,push_forward,id));
}



template <class RefSpace_, class PushForward_>
auto
PhysicalSpace<RefSpace_,PushForward_>::
begin() const -> ElementIterator
{
    return ElementIterator(this->shared_from_this(), 0);
}



template <class RefSpace_, class PushForward_>
auto
PhysicalSpace<RefSpace_,PushForward_>::
last() const -> ElementIterator
{
    return ElementIterator(this->shared_from_this(),
                           ref_space_->get_grid()->get_num_elements() - 1);
}



template <class RefSpace_, class PushForward_>
auto
PhysicalSpace<RefSpace_,PushForward_>::
end() const -> ElementIterator
{
    return ElementIterator(this->shared_from_this(),
                           IteratorState::pass_the_end);
}

template <class RefSpace_, class PushForward_>
auto
PhysicalSpace<RefSpace_,PushForward_>::
get_element(const Index elem_flat_id) const -> ElementAccessor
{
    Assert(elem_flat_id >= 0 && elem_flat_id < ref_space_->get_grid()->get_num_elements(),
           ExcIndexRange(elem_flat_id,0,ref_space_->get_grid()->get_num_elements()));

    auto elem = this->begin();
    for (int i = 0 ; i < elem_flat_id ; ++i)
        ++elem;

    return *elem;
}



template <class RefSpace_, class PushForward_>
Index
PhysicalSpace<RefSpace_,PushForward_>::
get_num_basis() const
{
    return ref_space_->get_num_basis();
}



template <class RefSpace_, class PushForward_>
int
PhysicalSpace<RefSpace_,PushForward_>::
get_num_basis_per_element() const
{
    return ref_space_->get_num_basis_per_element();
}



template <class RefSpace_, class PushForward_>
auto
PhysicalSpace<RefSpace_,PushForward_>::
get_push_forward() const -> shared_ptr<const PushForwardType>
{
    return shared_ptr<const PushForwardType>(push_forward_);
}



template <class RefSpace_, class PushForward_>
auto
PhysicalSpace<RefSpace_,PushForward_>::
get_reference_space() const -> shared_ptr<const RefSpace>
{
    return shared_ptr<const RefSpace>(ref_space_);
}


template <class RefSpace_, class PushForward_>
auto
PhysicalSpace<RefSpace_,PushForward_>::
get_face_space(const Index face_id,
               vector<Index> &face_to_element_dofs) const -> shared_ptr<FaceSpace>
{
    auto elem_map = std::make_shared<std::map<int,int> >();
    auto face_ref_sp = ref_space_->get_ref_face_space(face_id, face_to_element_dofs, *elem_map);
    auto map  = push_forward_->get_mapping();

    auto fmap = MappingSlice<FaceSpace::PushForwardType::dim, FaceSpace::PushForwardType::codim>::
    create(map, face_id, face_ref_sp->get_grid(), elem_map);
    auto fpf = FaceSpace::PushForwardType::create(fmap);
    auto face_space = FaceSpace::create(face_ref_sp,fpf);

    return face_space;
}

#if 0
template <class RefSpace_, class PushForward_>
auto
PhysicalSpace<RefSpace_,PushForward_>::
get_degree() const -> const ComponentTable<TensorIndex<dim>> &
{
    return ref_space_->get_degree();
}


template <class RefSpace_, class PushForward_>
auto
PhysicalSpace<RefSpace_,PushForward_>::
get_element_global_dofs() const -> const std::vector<std::vector<Index>> &
{
    return ref_space_->get_element_global_dofs();
}
#endif

template <class RefSpace_, class PushForward_>
Index
PhysicalSpace<RefSpace_,PushForward_>::
get_id() const
{
    return id_;
}

template <class RefSpace_, class PushForward_>
void
PhysicalSpace<RefSpace_,PushForward_>::
print_info(LogStream &out) const
{
    using std::endl;
    out << "PhysicalSpace info" << endl;

    out.push("\t");
    out << "Reference space:" << endl;
    ref_space_->print_info(out);
    out << endl;

    out << "Push-forward:" << endl;
    push_forward_->print_info(out);
    out << endl;

    out.pop();
}

template <class RefSpace_, class PushForward_>
std::shared_ptr<DofsManager>
PhysicalSpace<RefSpace_,PushForward_>::
get_dofs_manager() const
{
	return this->get_reference_space()->get_dofs_manager();
}


template <class RefSpace_, class PushForward_>
void
PhysicalSpace<RefSpace_,PushForward_>::
print_memory_info(LogStream &out) const
{
    using std::endl;
    out << "PHYSICAL SPACE memory info" << endl;
    out << "this address = " << this << endl;

    out.push("\t");
    out << "ref_space_ memory address = " << ref_space_ << endl;
    out << endl;

    out << "push_forward_ memory address = " << push_forward_ << endl;
    push_forward_->print_memory_info(out);
    out << endl;

    out.pop();
}

IGA_NAMESPACE_CLOSE

#include <igatools/basis_functions/physical_space.inst>


