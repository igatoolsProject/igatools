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
#include <igatools/utils/multi_array_utils.h>


using std::shared_ptr;
using std::make_shared;
using std::const_pointer_cast;

IGA_NAMESPACE_OPEN


template<int dim, int range, int rank>
ReferenceSpace<dim, range, rank>::
ReferenceSpace(const std::shared_ptr<SpaceData> space_data)
    :
    GridSpace(std::const_pointer_cast<CartesianGrid<dim>>(space_data->get_grid())),
    space_data_(space_data)
{
    Assert(this->get_grid() != nullptr,ExcNullPtr());
    Assert(space_data_ != nullptr,ExcNullPtr());
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
get_space_data() const -> std::shared_ptr<SpaceData>
{
    Assert(space_data_ != nullptr,ExcNullPtr());
    return space_data_;
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
get_element_dofs(
    const CartesianGridElement<dim> &element,
    const DofDistribution<dim, range, rank> &dofs_distribution) const
{
    const auto &accum_mult = this->space_data_->accumulated_interior_multiplicities();
    const auto &index_table = dofs_distribution.get_index_table();

    const auto &degree_table = this->space_data_->get_degree();

    vector<Index> element_dofs;
    const auto &elem_tensor_id = element.get_tensor_index();

    using Topology = UnitElement<dim>;

    for (int comp = 0 ; comp < SpaceData::n_components ; ++comp)
    {
        //-----------------------------------------------------------------
        // building the lookup table for the local dof id on the current component of the element --- begin
        // TODO (MM, March 06, 2015): this can be put on the SplineSpace constructor for optimization
        const auto &degree_comp = degree_table[comp];

        TensorSize<dim> dofs_t_size_elem_comp;
        for (const auto dir : Topology::active_directions)
            dofs_t_size_elem_comp[dir] = degree_comp[dir] + 1;

        const auto dofs_f_size_elem_comp = dofs_t_size_elem_comp.flat_size();

        vector<Index> elem_comp_dof_f_id(dofs_f_size_elem_comp);
        std::iota(elem_comp_dof_f_id.begin(),elem_comp_dof_f_id.end(),0);

        vector<TensorIndex<dim>> elem_comp_dof_t_id;
        const auto w_dofs_elem_comp = MultiArrayUtils<dim>::compute_weight(dofs_t_size_elem_comp);
        for (const auto dof_f_id : elem_comp_dof_f_id)
            elem_comp_dof_t_id.emplace_back(MultiArrayUtils<dim>::flat_to_tensor_index(dof_f_id,w_dofs_elem_comp));
        // building the lookup table for the local dof id on the current component of the element --- end
        //-----------------------------------------------------------------



        //-----------------------------------------------------------------
        const auto &index_table_comp = index_table[comp];

        const auto dof_t_origin = accum_mult[comp].cartesian_product(elem_tensor_id);
        for (const auto loc_dof_t_id : elem_comp_dof_t_id)
        {
            const auto dof_t_id = dof_t_origin + loc_dof_t_id;
            element_dofs.emplace_back(index_table_comp(dof_t_id));
        }
        //-----------------------------------------------------------------

    } // end comp loop

    return element_dofs;
}



template<int dim, int range, int rank>
vector<Index>
ReferenceSpace<dim, range, rank>::
get_loc_to_global(const CartesianGridElement<dim> &element) const
{
    return this->get_element_dofs(
               element,
               this->get_dof_distribution_global());
}



template<int dim, int range, int rank>
vector<Index>
ReferenceSpace<dim, range, rank>::
get_loc_to_patch(const CartesianGridElement<dim> &element) const
{
    const auto elem_dofs_global = this->get_loc_to_global(element);
    vector<Index> elem_dofs_local;

    const auto &dof_distribution = this->get_dof_distribution_global();
    for (const auto dof_global : elem_dofs_global)
        elem_dofs_local.push_back(
            dof_distribution.global_to_patch_local(dof_global));

    return elem_dofs_local;
}


template<int dim, int range, int rank>
void
ReferenceSpace<dim, range, rank>::
add_dofs_offset(const Index offset)
{
    this->get_dof_distribution_global().add_dofs_offset(offset);
}


template<int dim, int range, int rank>
Index
ReferenceSpace<dim, range, rank>::
get_global_dof_id(const TensorIndex<dim> &tensor_index,
                  const Index comp) const
{
    return this->get_dof_distribution_global().get_index_table()[comp](tensor_index);
}


template<int dim, int range, int rank>
auto
ReferenceSpace<dim, range, rank>::
get_dof_distribution_global() const -> const DofDistribution<dim, range, rank> &
{
    return space_data_->get_dof_distribution_global();
}



template<int dim, int range, int rank>
auto
ReferenceSpace<dim, range, rank>::
get_dof_distribution_global() -> DofDistribution<dim, range, rank> &
{
    return space_data_->get_dof_distribution_global();
}


IGA_NAMESPACE_CLOSE

#include <igatools/basis_functions/reference_space.inst>

