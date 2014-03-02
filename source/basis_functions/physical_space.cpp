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

using std::array;
using std::shared_ptr;
using std::const_pointer_cast;

IGA_NAMESPACE_OPEN

template <class RefSpace_, class PushForward_>
PhysicalSpace<RefSpace_,PushForward_>::
PhysicalSpace(
    shared_ptr<RefSpace> ref_space,
    shared_ptr<PushForwardType> push_forward)
    :
    BaseSpace(ref_space->get_grid()),
    ref_space_(ref_space),
    push_forward_(push_forward)
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
            shared_ptr<PushForwardType>(new PushForwardType(*push_forward_))));
};



template <class RefSpace_, class PushForward_>
auto
PhysicalSpace<RefSpace_,PushForward_>::
create(
    shared_ptr<RefSpace> ref_space,
    shared_ptr<PushForwardType> push_forward) -> shared_ptr<self_t>
{
    Assert(ref_space != nullptr, ExcNullPtr());
    Assert(push_forward != nullptr, ExcNullPtr());
    return shared_ptr<self_t>(new self_t(ref_space,push_forward));
}



template <class RefSpace_, class PushForward_>
auto
PhysicalSpace<RefSpace_,PushForward_>::
begin() const -> ElementIterator
{
    return ElementIterator(const_cast<self_t &>(*this), 0);
}



template <class RefSpace_, class PushForward_>
auto
PhysicalSpace<RefSpace_,PushForward_>::
last() const -> ElementIterator
{
    return ElementIterator(
               const_cast<self_t &>(*this),
               ref_space_->get_grid()->get_num_elements() - 1);
}



template <class RefSpace_, class PushForward_>
auto
PhysicalSpace<RefSpace_,PushForward_>::
end() const -> ElementIterator
{
    return ElementIterator(const_cast<self_t & >(*this),
                           IteratorState::pass_the_end);
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


