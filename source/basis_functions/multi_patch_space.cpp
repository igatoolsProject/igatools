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
#include <igatools/basis_functions/dofs_manager.h>
#include <igatools/utils/vector_tools.h>

#include <igatools/utils/dynamic_multi_array.h>

#include <igatools/base/logstream.h>

using std::string;
using std::vector;
using std::shared_ptr;
using std::unique_ptr;


IGA_NAMESPACE_OPEN

template <class PhysicalSpace>
MultiPatchSpace<PhysicalSpace>::
MultiPatchSpace(shared_ptr<DofsManager> dofs_manager)
    :
    dofs_manager_(dofs_manager)
{
    Assert(dofs_manager != nullptr,ExcNullPtr());
}


template <class PhysicalSpace>
void
MultiPatchSpace<PhysicalSpace>::
patch_insertion_open()
{
    is_patch_insertion_open_ = true;

    dofs_manager_->space_insertion_open();
}

template <class PhysicalSpace>
void
MultiPatchSpace<PhysicalSpace>::
patch_insertion_close(const bool automatic_dofs_renumbering)
{
    //------------------------------------------------------------------------
    // check that each reference space is used in only one physical space -- begin
    vector<shared_ptr<const RefSpace>> ref_spaces;
    for (const auto &phys_space : patches_)
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
    for (const auto &phys_space : patches_)
        maps.push_back(phys_space->get_push_forward()->get_mapping());

    vector<shared_ptr<const Map>> maps_no_duplicates;
    vector<int> maps_multiplicities;
    vector_tools::count_and_remove_duplicates(
        maps,maps_no_duplicates,maps_multiplicities) ;

    for (const int mult : maps_multiplicities)
        AssertThrow(mult == 1,ExcMessage("At least one mapping is used to define multiple physical spaces."));
    // check that a mapping is used in only one reference space -- end
    //------------------------------------------------------------------------


    dofs_manager_->space_insertion_close(automatic_dofs_renumbering);

    is_patch_insertion_open_ = false;
}

template <class PhysicalSpace>
void
MultiPatchSpace<PhysicalSpace>::
build_graph()
{
    Assert(is_patch_insertion_open_ == false,ExcInvalidState());
    Assert(is_interface_insertion_open_ == false,ExcInvalidState());


    //------------------------------------------------------------------------
    // creating the graph representing the multipatch structure -- begin
    // the vertices of the graph represents the patches
    std::map<PatchPtr,Index> map_patch_id;
    Index id = 0;
    for (const auto &patch : patches_)
    {
        boost::add_vertex(patch, multipatch_graph_);

        map_patch_id.emplace(patch,id++);
    }


    // the edges of the graph represents the interfaces
    for (const auto &interfaces_same_type : interfaces_)
    {
        for (const auto &interface : interfaces_same_type.second)
        {
            const Index id_patch_0 = map_patch_id.at(interface->patch_and_side_[0].first);
            const Index id_patch_1 = map_patch_id.at(interface->patch_and_side_[1].first);

            boost::add_edge(
                boost::vertex(id_patch_0,multipatch_graph_),
                boost::vertex(id_patch_1,multipatch_graph_),
                interface,
                multipatch_graph_);
        }
    }
    Assert(patches_.size() == boost::num_vertices(multipatch_graph_),
           ExcDimensionMismatch(patches_.size(),boost::num_vertices(multipatch_graph_)));
    Assert(get_num_interfaces() == boost::num_edges(multipatch_graph_),
           ExcDimensionMismatch(get_num_interfaces(),boost::num_edges(multipatch_graph_)));

    // creating the graph representing the multipatch structure -- end
    //------------------------------------------------------------------------

    is_graph_built_ = true;
}

template <class PhysicalSpace>
void
MultiPatchSpace<PhysicalSpace>::
compute_constraints()
{
    Assert(is_patch_insertion_open_ == false,ExcInvalidState());
    Assert(is_interface_insertion_open_ == false,ExcInvalidState());
    Assert(are_constraints_computed_ == false,ExcInvalidState());

    //---------------------------------------------------------------------------
    if (is_graph_built_ == false)
        this->build_graph();
    //---------------------------------------------------------------------------



    Assert(false,ExcNotImplemented());
    AssertThrow(false,ExcNotImplemented());

    are_constraints_computed_ = true;
}


template <class PhysicalSpace>
void
MultiPatchSpace<PhysicalSpace>::
perform_ref_spaces_add_dofs_offset()
{
    Index dofs_offset = 0;
    for (const auto &phys_space : patches_)
    {
        shared_ptr<RefSpace> ref_space = std::const_pointer_cast<RefSpace>(phys_space->get_reference_space());
        ref_space->add_dofs_offset(dofs_offset);

        dofs_offset += ref_space->get_num_basis();
    }
}

template <class PhysicalSpace>
void
MultiPatchSpace<PhysicalSpace>::
add_patch(PatchPtr patch)
{
    Assert(is_patch_insertion_open_ == true,ExcInvalidState());

    //------------------------------------------------------------------------
    // check if the patch is already present in the vector of patches -- begin
    Assert(std::count(patches_.begin(),patches_.end(),patch) == 0,
           ExcMessage("Patch (physical space) already added."))
    // check if the patch is already present in the vector of patches -- end
    //------------------------------------------------------------------------

    patches_.push_back(patch);


    dofs_manager_->add_space(patch);
}


template <class PhysicalSpace>
auto
MultiPatchSpace<PhysicalSpace>::
get_patches() const -> vector<PatchPtr>
{
    return patches_;
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
get_num_interfaces(const InterfaceType interface_type) const
{
    return const_cast<std::map<InterfaceType,std::set<InterfacePtr>> &>(interfaces_)[interface_type].size();
}


template <class PhysicalSpace>
int
MultiPatchSpace<PhysicalSpace>::
get_num_interfaces() const
{
    Index n_interfaces = 0;
    for (const auto &interfaces_same_type : interfaces_)
    {
        n_interfaces += interfaces_same_type.second.size();
    }

    return n_interfaces;
}

template <class PhysicalSpace>
auto
MultiPatchSpace<PhysicalSpace>::
get_dofs_manager() -> shared_ptr<DofsManager>
{
    /*
    Assert(is_patch_insertion_open_ == false,ExcInvalidState());
    Assert(is_interface_insertion_open_ == false,ExcInvalidState());
    Assert(is_graph_built_ == true,ExcInvalidState());
    Assert(are_constraints_computed_ == true,ExcInvalidState());
    //*/
    return dofs_manager_;
}

template <class PhysicalSpace>
auto
MultiPatchSpace<PhysicalSpace>::
get_dofs_manager() const -> shared_ptr<const DofsManager>
{
    /*
    Assert(is_patch_insertion_open_ == false,ExcInvalidState());
    Assert(is_interface_insertion_open_ == false,ExcInvalidState());
    Assert(is_graph_built_ == true,ExcInvalidState());
    Assert(are_constraints_computed_ == true,ExcInvalidState());
    //*/
    return dofs_manager_;
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

    Assert(is_patch_insertion_open_ == false,ExcInvalidState());
    out << "Num. patches = " << this->get_num_patches() << endl;

    out.push(tab);
    for (const auto &patch : patches_)
    {
        out << "Patch id = " << patch->get_id() << endl;
//        patch->print_info(out);
        out.push(tab);
    }



    out.pop();

    Assert(is_interface_insertion_open_ == false,ExcInvalidState());
    out << "Num. interfaces = " << this->get_num_interfaces() << endl;
    int interface_id = 0 ;
    for (const auto &interfaces_same_type : interfaces_)
    {
        out << "Interfaces of type " << to_integral(interfaces_same_type.first)
            << " --- num=" << interfaces_same_type.second.size() << endl;

        out.push(tab);
        for (const auto &interface : interfaces_same_type.second)
        {
            out << "Interface id = " << interface_id++ << endl;
            interface->print_info(out);
            out.push(tab);
        }
        out.pop();
    }




    //---------------------------------------------------------------------------
    Assert(is_graph_built_ == true,ExcInvalidState());
    out << "Patches in the graph:" << endl;
    using vertex_iterator = typename boost::graph_traits<Graph>::vertex_iterator;
    vertex_iterator vertex;
    vertex_iterator vertex_end;
    boost::tie(vertex, vertex_end) = boost::vertices(multipatch_graph_);

    out.push(tab);
    for (; vertex != vertex_end ; ++vertex)
    {
        //printing the id of the patch represented by the vertex
        out << "Patch ID = " << multipatch_graph_[*vertex]->get_id() << endl;
    }
    out.pop();
    //---------------------------------------------------------------------------



    //---------------------------------------------------------------------------
    out << "Interfaces in the graph:" << endl;
    using edge_iterator = typename boost::graph_traits<Graph>::edge_iterator;
    edge_iterator edge;
    edge_iterator edge_end;
    boost::tie(edge, edge_end) = boost::edges(multipatch_graph_);

    out.push(tab);
    for (; edge != edge_end ; ++edge)
    {
        //printing the information about the interface represented by the edge
        multipatch_graph_[*edge]->print_info(out);
        out << endl;
    }
    out.pop();
    //---------------------------------------------------------------------------



    //---------------------------------------------------------------------------
    out.push(tab);
    out << "DOFs manager:" << endl;
    out.push(tab);
    dofs_manager_->print_info(out);
    out << endl;
    out.pop();
    //---------------------------------------------------------------------------


    out.pop();
}

template <class PhysicalSpace>
void
MultiPatchSpace<PhysicalSpace>::
interface_insertion_open()
{
    Assert(is_patch_insertion_open_ == false,ExcInvalidState());
    is_interface_insertion_open_ = true;
}

template <class PhysicalSpace>
void
MultiPatchSpace<PhysicalSpace>::
interface_insertion_close()
{
    is_interface_insertion_open_ = false;
}

template <class PhysicalSpace>
void
MultiPatchSpace<PhysicalSpace>::
add_interface(const InterfaceType &type,
              PatchPtr patch_0,const int side_id_patch_0,
              PatchPtr patch_1,const int side_id_patch_1)
{
    Assert(is_interface_insertion_open_ == true,ExcInvalidState());

    //------------------------------------------------------------------------
    // Verify that patch 0 is present in the vector of patches -- begin
    Assert(std::count(patches_.begin(),patches_.end(),patch_0) == 1,
           ExcMessage("Patch 0 is not present in the vector of patches."))
    // Verify that patch 0 is present in the vector of patches -- end
    //------------------------------------------------------------------------


    //------------------------------------------------------------------------
    // Verify that patch 1 is present in the vector of patches -- begin
    Assert(std::count(patches_.begin(),patches_.end(),patch_1) == 1,
           ExcMessage("Patch 1 is not present in the vector of patches."))
    // Verify that patch 1 is present in the vector of patches -- end
    //------------------------------------------------------------------------


    InterfacePtr interface_to_be_added(
        new Interface(type,patch_0,side_id_patch_0,patch_1,side_id_patch_1));

#ifndef NDEBUG
    for (const auto &interface : interfaces_[type])
        Assert(*interface_to_be_added != *interface, ExcMessage("Interface already added."));
#endif

    interfaces_[type].insert(interface_to_be_added);
}


template <class PhysicalSpace>
MultiPatchSpace<PhysicalSpace>::
Interface::
Interface(const InterfaceType &type, PatchPtr patch_0,const int side_id_patch_0,PatchPtr patch_1,const int side_id_patch_1)
    :
    type_(type)
{
    Assert(patch_0 != patch_1,ExcMessage("Impossible to use the same patch to define an interface."));
    Assert(side_id_patch_0 >= 0 && side_id_patch_0 < (UnitElement<dim>::faces_per_element),
           ExcIndexRange(side_id_patch_0,0,UnitElement<dim>::faces_per_element));
    Assert(side_id_patch_1 >= 0 && side_id_patch_1 < (UnitElement<dim>::faces_per_element),
           ExcIndexRange(side_id_patch_1,0,UnitElement<dim>::faces_per_element));


    patch_and_side_[0] = std::make_pair(patch_0, side_id_patch_0);
    patch_and_side_[1] = std::make_pair(patch_1, side_id_patch_1);
}




template <class PhysicalSpace>
bool
MultiPatchSpace<PhysicalSpace>::
Interface::
operator==(const Interface &interface_to_compare) const
{
    return (type_ == interface_to_compare.type_ &&
            patch_and_side_[0] == interface_to_compare.patch_and_side_[0] &&
            patch_and_side_[1] == interface_to_compare.patch_and_side_[1]);
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
process_interfaces_C0_strong()
{
    Assert(is_patch_insertion_open_ == false,ExcInvalidState());
    Assert(is_interface_insertion_open_ == false,ExcInvalidState());
    Assert(is_graph_built_ == true,ExcInvalidState());

    Assert(false,ExcNotImplemented());
    AssertThrow(false,ExcNotImplemented());
}

template <class PhysicalSpace>
void
MultiPatchSpace<PhysicalSpace>::
process_interfaces_C0_strong_renumbering()
{
    Assert(is_patch_insertion_open_ == false,ExcInvalidState());
    Assert(is_interface_insertion_open_ == false,ExcInvalidState());
    Assert(is_graph_built_ == true,ExcInvalidState());

    Assert(false,ExcNotImplemented());
    AssertThrow(false,ExcNotImplemented());
}

template <class PhysicalSpace>
void
MultiPatchSpace<PhysicalSpace>::
process_interfaces_mortar()
{
    Assert(is_patch_insertion_open_ == false,ExcInvalidState());
    Assert(is_interface_insertion_open_ == false,ExcInvalidState());
    Assert(is_graph_built_ == true,ExcInvalidState());

    this->process_interfaces_C0_strong();
    this->process_interfaces_C0_strong_renumbering();
    this->process_interfaces_mortar();

    Assert(false,ExcNotImplemented());
    AssertThrow(false,ExcNotImplemented());
}


template <class PhysicalSpace>
void
MultiPatchSpace<PhysicalSpace>::
process_interfaces()
{
    Assert(false,ExcNotImplemented());
    AssertThrow(false,ExcNotImplemented());
}


template <class PhysicalSpace>
auto
MultiPatchSpace<PhysicalSpace>::
get_interfaces_same_type(const InterfaceType interface_type) -> std::set<InterfacePtr>
{
    return interfaces_[interface_type];
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
    out << "Patch id = " << patch_and_side_[0].first->get_id() << endl;
    out << "Side  id = " << patch_and_side_[0].second << endl;
    out.pop();

    out << "Patch 1 infos:" << endl;
    out.push(tab);
    out << "Patch id = " << patch_and_side_[1].first->get_id() << endl;
    out << "Side id = " << patch_and_side_[1].second << endl;
    out.pop();



    out.pop();

    out.pop();

}



IGA_NAMESPACE_CLOSE


#include <igatools/basis_functions/multi_patch_space.inst>
