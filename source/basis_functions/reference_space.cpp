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


#include <igatools/basis_functions/reference_space.h>
#include <igatools/basis_functions/reference_element.h>
#include <igatools/base/array_utils.h>

#include <igatools/basis_functions/bspline_space.h>
#include <igatools/basis_functions/nurbs_space.h>
#include <igatools/basis_functions/space_manager.h>


using std::shared_ptr;
using std::make_shared;
using std::const_pointer_cast;

IGA_NAMESPACE_OPEN


template<int dim_, int range_, int rank_>
ReferenceSpace<dim_, range_, rank_>::
ReferenceSpace(
    const std::shared_ptr<CartesianGrid<dim_>> grid,
    const std::shared_ptr<DofDistribution<dim_,range_,rank_>> dof_distribution)
    :
    GridSpace(grid),
    dof_distribution_(dof_distribution)
{
    Assert(this->get_grid() != nullptr,ExcNullPtr());
    Assert(dof_distribution_ != nullptr,ExcNullPtr());

    //------------------------------------------------------------------------------
    /*
    using DofDistrib = DofDistribution<dim,range,rank>;
    dof_distribution_ = shared_ptr<DofDistrib>(new DofDistrib(
                                                   space_data_->get_num_basis_table(),
                                                   space_data_->get_degree(),
                                                   space_data_->get_periodic_table()));
    //*/
    //------------------------------------------------------------------------------
}



template<int dim, int range, int rank>
template<int k>
auto
ReferenceSpace<dim, range, rank>::
get_ref_sub_space(const int sub_elem_id,
                  InterSpaceMap<k> &dof_map,
                  std::shared_ptr<CartesianGrid<k>> sub_grid) const
-> std::shared_ptr< SubRefSpace<k> >
{
    std::shared_ptr< SubRefSpace<k> > sub_ref_space;
    if (this->is_bspline())
    {
        const auto bsp_space = dynamic_cast<const BSplineSpace<dim,range,rank> *>(this);
        Assert(bsp_space != nullptr,ExcNullPtr());
        sub_ref_space = bsp_space->get_ref_sub_space(sub_elem_id,dof_map,sub_grid);
    }
    else
    {
#ifdef NURBS
        //TODO (MM, Dec 22, 2014): implement NURBSSpace::get_ref_sub_space()
        const auto nrb_space = dynamic_cast<const NURBSSpace<dim,range,rank> *>(this);
        Assert(nrb_space != nullptr,ExcNullPtr());
        Assert(false,ExcNotImplemented());
#else
        Assert(false,ExcMessage("NURBS support disabled from configuration cmake parameters."));
        AssertThrow(false,ExcMessage("NURBS support disabled from configuration cmake parameters."));
#endif
    }

    Assert(sub_ref_space != nullptr, ExcNullPtr());
    return sub_ref_space;
}



template<int dim, int range, int rank>
template<int k>
auto
ReferenceSpace<dim, range, rank>::
get_sub_space(const int s_id, InterSpaceMap<k> &dof_map,
              std::shared_ptr<CartesianGrid<k>> sub_grid,
              std::shared_ptr<InterGridMap<k>> elem_map) const
-> std::shared_ptr<SubSpace<k> >
{
    std::shared_ptr<SubSpace<k> > sub_space;
    if (this->is_bspline())
    {
        const auto bsp_space =
        dynamic_cast<const BSplineSpace<dim,range,rank> *>(this);
        Assert(bsp_space != nullptr, ExcNullPtr());
        sub_space = bsp_space->get_sub_space(s_id,dof_map,sub_grid, elem_map);
    }
    else
    {
#ifdef NURBS

        Assert(false,ExcNotImplemented());
#else
        Assert(false,ExcMessage("NURBS support disabled from configuration cmake parameters."));
        AssertThrow(false,ExcMessage("NURBS support disabled from configuration cmake parameters."));
#endif
    }

    Assert(sub_space != nullptr, ExcNullPtr());
    return sub_space;
}


template<int dim, int range, int rank>
auto
ReferenceSpace<dim, range, rank>::
begin(const std::string &element_property) const -> ElementIterator
{
    return ElementIterator(
               this->create_element(
                   this->get_grid()->get_first_element_id_same_property(element_property)),
               element_property);
}



template<int dim, int range, int rank>
auto
ReferenceSpace<dim, range, rank>::
last(const std::string &element_property) const -> ElementIterator
{
    return ElementIterator(
               this->create_element(
                   this->get_grid()->get_last_element_id_same_property(element_property)),
               element_property);
}



template<int dim, int range, int rank>
auto
ReferenceSpace<dim, range, rank>::
end(const std::string &element_property) const -> ElementIterator
{
    return ElementIterator(this->create_element(IteratorState::pass_the_end),element_property);
}



template<int dim, int range, int rank>
auto
ReferenceSpace<dim, range, rank>::
get_space_manager() -> shared_ptr<SpaceManager>
{
    auto space_manager = make_shared<SpaceManager>(SpaceManager());


    shared_ptr<ReferenceSpace<dim,range,rank> > this_space;
    if (this->is_bspline())
    {
        using BSpSpace = BSplineSpace<dim,range,rank>;
        this_space = dynamic_cast<BSpSpace &>(*this).shared_from_this();
    }
    else
    {
#ifdef NURBS
        using NrbSpace = NURBSSpace<dim,range,rank>;
        this_space = dynamic_cast<NrbSpace &>(*this).shared_from_this();
#else
        Assert(false,ExcMessage("NURBS support disabled from configuration cmake parameters."));
        AssertThrow(false,ExcMessage("NURBS support disabled from configuration cmake parameters."));
#endif
    }
    Assert(this_space != nullptr,ExcNullPtr());

    space_manager->spaces_insertion_open();
    space_manager->add_space(this_space);
    space_manager->spaces_insertion_close();


    space_manager->spaces_connectivity_open();
    space_manager->add_spaces_connection(this_space);
    space_manager->spaces_connectivity_close();

    return space_manager;
}



template<int dim, int range, int rank>
auto
ReferenceSpace<dim, range, rank>::
get_space_manager() const -> std::shared_ptr<const SpaceManager>
{
    return const_cast<ReferenceSpace<dim,range,rank> &>(*this).get_space_manager();
}






template<int dim, int range, int rank>
vector<Index>
ReferenceSpace<dim, range, rank>::
get_loc_to_global(const CartesianGridElement<dim> &element) const
{
    return this->get_element_dofs(element);
}



template<int dim, int range, int rank>
vector<Index>
ReferenceSpace<dim, range, rank>::
get_loc_to_patch(const CartesianGridElement<dim> &element) const
{
    const auto elem_dofs_global = this->get_loc_to_global(element);
    vector<Index> elem_dofs_local;

    for (const auto dof_global : elem_dofs_global)
        elem_dofs_local.push_back(
            dof_distribution_->global_to_patch_local(dof_global));

    return elem_dofs_local;
}


template<int dim, int range, int rank>
void
ReferenceSpace<dim, range, rank>::
add_dofs_offset(const Index offset)
{
    dof_distribution_->add_dofs_offset(offset);
}


template<int dim, int range, int rank>
Index
ReferenceSpace<dim, range, rank>::
get_global_dof_id(const TensorIndex<dim> &tensor_index,
                  const Index comp) const
{
    return dof_distribution_->get_index_table()[comp](tensor_index);
}


template<int dim, int range, int rank>
auto
ReferenceSpace<dim, range, rank>::
get_dof_distribution() const -> shared_ptr<const DofDistribution<dim,range,rank> >
{
    return dof_distribution_;
}

template<int dim, int range, int rank>
auto
ReferenceSpace<dim, range, rank>::
get_dof_distribution() -> shared_ptr<DofDistribution<dim,range,rank> >
{
    return dof_distribution_;
}


IGA_NAMESPACE_CLOSE

#include <igatools/basis_functions/reference_space.inst>

