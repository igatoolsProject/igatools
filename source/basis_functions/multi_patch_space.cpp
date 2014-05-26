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


#include <igatools/basis_functions/multi_patch_space.h>
#include <igatools/basis_functions/physical_space.h>
#include <igatools/utils/vector_tools.h>

#include <igatools/base/logstream.h>

using std::string;
using std::vector;
using std::shared_ptr;
using std::unique_ptr;


IGA_NAMESPACE_OPEN



template <class PhysicalSpace>
void
MultiPatchSpace<PhysicalSpace>::
arrangement_open()
{
    is_arrangement_open_ = true;
}

template <class PhysicalSpace>
void
MultiPatchSpace<PhysicalSpace>::
arrangement_close()
{
    is_arrangement_open_ = false;


    //------------------------------------------------------------------------
    // check that each reference space is used in only one physical space -- begin
    vector<shared_ptr<const RefSpace>> ref_spaces;
    for (const auto & phys_space : patches_)
        ref_spaces.push_back(phys_space->get_reference_space());


    vector<shared_ptr<const RefSpace>> ref_spaces_no_duplicates;
    vector<int> ref_spaces_multiplicities;
    vector_tools::count_and_remove_duplicates(
        ref_spaces,ref_spaces_no_duplicates,ref_spaces_multiplicities) ;

    for (const int mult : ref_spaces_multiplicities)
        AssertThrow(mult == 1,ExcMessage("At least one reference space is used to define multiple physical spaces."));
    // check that a reference space is used in only one physical space -- end
    //------------------------------------------------------------------------


    //------------------------------------------------------------------------
    // check that a each mapping is used in only one reference space -- begin
    vector<shared_ptr<const Map>> maps;
    for (const auto & phys_space : patches_)
        maps.push_back(phys_space->get_push_forward()->get_mapping());

    vector<shared_ptr<const Map>> maps_no_duplicates;
    vector<int> maps_multiplicities;
    vector_tools::count_and_remove_duplicates(
        maps,maps_no_duplicates,maps_multiplicities) ;

    for (const int mult : maps_multiplicities)
        AssertThrow(mult == 1,ExcMessage("At least one mapping is used to define multiple physical spaces."));
    // check that a mapping is used in only one reference space -- end
    //------------------------------------------------------------------------


    //------------------------------------------------------------------------
    // Renumber the dofs in the reference spaces in order to avoid same dof ids between different spaces -- begin
    this->perform_ref_spaces_dofs_renumbering();
    // Renumber the dofs in the reference spaces in order to avoid same dof ids between different spaces -- end
    //------------------------------------------------------------------------
}


template <class PhysicalSpace>
void
MultiPatchSpace<PhysicalSpace>::
perform_ref_spaces_dofs_renumbering()
{
    Index dofs_offset = 0;
    for (const auto & phys_space : patches_)
    {
        shared_ptr<RefSpace> ref_space = std::const_pointer_cast<RefSpace>(phys_space->get_reference_space());
        ref_space->add_dofs_offset(dofs_offset);

        dofs_offset += ref_space->get_num_basis();
    }
}

template <class PhysicalSpace>
void
MultiPatchSpace<PhysicalSpace>::
add_patch(std::shared_ptr<const PhysicalSpace> patch)
{
    Assert(is_arrangement_open_,ExcInvalidState());

    //------------------------------------------------------------------------
    // check if the patch is already present in the vector of patches -- begin
    Assert(std::count(patches_.begin(),patches_.end(),patch) == 0,
           ExcMessage("Patch (physical space) already added."))
    // check if the patch is already present in the vector of patches -- end
    //------------------------------------------------------------------------

    patches_.push_back(patch);
}


template <class PhysicalSpace>
int
MultiPatchSpace<PhysicalSpace>::
get_num_patches() const
{
    return patches_.size();
}

template <class PhysicalSpace>
int
MultiPatchSpace<PhysicalSpace>::
get_num_interfaces() const
{
    return interfaces_.size();
}

template <class PhysicalSpace>
void
MultiPatchSpace<PhysicalSpace>::
print_info(LogStream &out) const
{
    using std::endl;
    string tab = "   ";

    out<< "MultiPatchSpace infos:" << endl;
    out.push(tab);

    out << "Num. patches = " << this->get_num_patches() << endl;

    out.push(tab);
    int patch_id = 0 ;
    for (const auto & patch : patches_)
    {
        out << "Patch id = " << patch_id++ << endl;
        patch->print_info(out);
        out.push(tab);
    }



    out.pop();

    out << "Num. interfaces = " << this->get_num_interfaces() << endl;
    int interface_id = 0 ;
    for (const auto & interface : interfaces_)
    {
        out << "Interface id = " << interface_id++ << endl;
        interface->print_info(out);
        out.push(tab);
    }


    out.pop();
}

template <class PhysicalSpace>
void
MultiPatchSpace<PhysicalSpace>::
add_interface(const InterfaceType &type,
              Patch patch_0,const int side_id_patch_0,
              Patch patch_1,const int side_id_patch_1)
{
    unique_ptr<Interface> interface_to_be_added(
        new Interface(type,patch_0,side_id_patch_0,patch_1,side_id_patch_1));

#ifndef NDEBUG
    for (const auto & interface : interfaces_)
        Assert(*interface_to_be_added != *interface, ExcMessage("Interface already added."));
#endif

    interfaces_.push_back(std::move(interface_to_be_added));
}


template <class PhysicalSpace>
MultiPatchSpace<PhysicalSpace>::
Interface::
Interface(const InterfaceType &type, Patch patch_0,const int side_id_patch_0,Patch patch_1,const int side_id_patch_1)
    :
    type_(type)
{
    Assert(patch_0 != patch_1,ExcMessage("Impossible to use the same patch to define an interface."));
    Assert(side_id_patch_0 >= 0 && side_id_patch_0 < (UnitElement<dim>::faces_per_element),
           ExcIndexRange(side_id_patch_0,0,UnitElement<dim>::faces_per_element));
    Assert(side_id_patch_1 >= 0 && side_id_patch_1 < (UnitElement<dim>::faces_per_element),
           ExcIndexRange(side_id_patch_1,0,UnitElement<dim>::faces_per_element));


    patch_[0] = patch_0;
    side_id_[0] = side_id_patch_0;

    patch_[1] = patch_1;
    side_id_[1] = side_id_patch_1;
}




template <class PhysicalSpace>
bool
MultiPatchSpace<PhysicalSpace>::
Interface::
operator==(const Interface &interface_to_compare) const
{
    return (type_ == interface_to_compare.type_ &&
            patch_[0] == interface_to_compare.patch_[0] &&
            patch_[1] == interface_to_compare.patch_[1] &&
            side_id_[0] == interface_to_compare.side_id_[0] &&
            side_id_[1] == interface_to_compare.side_id_[1]);
}

template <class PhysicalSpace>
bool
MultiPatchSpace<PhysicalSpace>::
Interface::
operator!=(const Interface &interface_to_compare) const
{
    return !(*this == interface_to_compare);
}


template <class PhysicalSpace>
void
MultiPatchSpace<PhysicalSpace>::
Interface::
print_info(LogStream &out) const
{
    using std::endl;
    string tab = "   ";

    out << "Interface type = " << static_cast<int>(type_) << endl;
    out.push(tab);

    out << "Patch 0 infos:" << endl;
    out.push(tab);
    out << "shared_ptr = " << patch_[0] << endl;
    out << "Side id = " << side_id_[0] << endl;
    out.pop();

    out << "Patch 1 infos:" << endl;
    out.push(tab);
    out << "shared_ptr = " << patch_[1] << endl;
    out << "Side id = " << side_id_[1] << endl;
    out.pop();




    out.pop();

    out.pop();

}



IGA_NAMESPACE_CLOSE


#include <igatools/basis_functions/multi_patch_space.inst>
