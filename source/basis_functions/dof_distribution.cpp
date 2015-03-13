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


#include <igatools/basis_functions/dof_distribution.h>
#include <igatools/utils/multi_array_utils.h>


using std::map;
using std::shared_ptr;
using std::make_shared;

IGA_NAMESPACE_OPEN

template<int dim, int range, int rank>
DofDistribution<dim, range, rank>::
DofDistribution(const TensorSizeTable &n_basis,
                const DegreeTable &degree_table,
                const PeriodicTable &periodic,
                DistributionPolicy pol)
    :
    num_dofs_table_(n_basis),
    index_table_size_(n_basis),
    policy_(pol)
{
    Assert(pol == DistributionPolicy::standard, ExcNotImplemented());

    using Topology = UnitElement<dim>;

    //-----------------------------------------------------------------------
    for (const auto comp : Space::components)
        for (const auto dir : Topology::active_directions)
            if (periodic[comp][dir])
                index_table_size_[comp][dir] += degree_table[comp][dir] + 1;
    //-----------------------------------------------------------------------


    //-----------------------------------------------------------------------
    // fills the standard distribution, sorted by component and
    // by direction x moves faster
    Index comp_offset = 0;
    for (const auto comp : Space::components)
    {
        const auto &index_table_size_comp = index_table_size_[comp];
        auto &index_table_comp = index_table_[comp];
        index_table_comp.resize(index_table_size_comp);

        const auto &n_dofs_comp = num_dofs_table_[comp];
        const auto w_dofs_comp = MultiArrayUtils<dim>::compute_weight(n_dofs_comp);


        const auto n_indices_comp = index_table_size_.get_component_size(comp);
        for (int i = 0 ; i < n_indices_comp ; ++i)
        {
            auto t_ind = index_table_comp.flat_to_tensor(i);
            for (const auto dir : Topology::active_directions)
                t_ind[dir] %= n_dofs_comp[dir];

            const auto f_ind = MultiArrayUtils<dim>::tensor_to_flat_index(t_ind,w_dofs_comp);
            index_table_comp[i] = comp_offset + f_ind;

        }
        comp_offset += num_dofs_table_.get_component_size(comp);//n_basis_comp_size;
    }
    //-----------------------------------------------------------------------


    //-----------------------------------------------------------------------
    // creating the dofs view from the dofs components views -- begin
    vector<DofsComponentView> components_views;
    for (auto &index_table_comp : index_table_)
        components_views.emplace_back(index_table_comp.get_flat_view());

    dofs_view_ = DofsView(
                     DofsIterator(components_views,0),
                     DofsIterator(components_views,IteratorState::pass_the_end));
    // creating the dofs view from the dofs components views -- end
    //-----------------------------------------------------------------------
}



template<int dim, int range, int rank>
Index
DofDistribution<dim, range, rank>::
get_min_dof_id() const
{
    return *std::min_element(dofs_view_.cbegin(),dofs_view_.cend());
}



template<int dim, int range, int rank>
Index
DofDistribution<dim, range, rank>::
get_max_dof_id() const
{
    return *std::max_element(dofs_view_.cbegin(),dofs_view_.cend());
}



template<int dim, int range, int rank>
bool
DofDistribution<dim, range, rank>::
find_dof_id(const Index dof_id, int &comp_id, TensorIndex<dim> &tensor_index) const
{
    bool dof_is_found = false;

    for (const auto comp: Space::components)
    {
        const auto &index_table_comp = index_table_[comp];

        const auto &index_table_comp_data = index_table_comp.get_data();

        const auto dofs_begin = index_table_comp_data.begin();
        const auto dofs_end   = index_table_comp_data.end();
        const auto it = std::find(dofs_begin, dofs_end, dof_id);

        if (it != dofs_end)
        {
            dof_is_found = true;
            comp_id = comp;
            tensor_index = index_table_comp.flat_to_tensor(it-dofs_begin);

            break;
        }
    }

    return dof_is_found;
}



template<int dim, int range, int rank>
void
DofDistribution<dim, range, rank>::
add_dofs_offset(const Index offset)
{
    for (auto &dofs_component : index_table_)
        for (auto &dof_id : dofs_component)
            dof_id += offset;
}



template<int dim, int range, int rank>
auto
DofDistribution<dim, range, rank>::
get_index_table() const -> const IndexDistributionTable &
{
    return index_table_;
}



template<int dim, int range, int rank>
auto
DofDistribution<dim, range, rank>::
get_dofs_view() -> DofsView &
{
    return dofs_view_;
}

template<int dim, int range, int rank>
auto
DofDistribution<dim, range, rank>::
get_dofs_view() const -> const DofsView &
{
    return dofs_view_;
}


template<int dim, int range, int rank>
Index
DofDistribution<dim, range, rank>::
global_to_patch_local(const Index global_dof_id) const
{
    int comp_id;
    TensorIndex<dim> tensor_index;
#ifndef NDEBUG
    bool global_dof_is_found = this->find_dof_id(global_dof_id,comp_id,tensor_index);
    Assert(global_dof_is_found,
           ExcMessage("The global dof id " + std::to_string(global_dof_id) +
                      " is not present in the DofDistribution."));
#else
    this->find_dof_id(global_dof_id,comp_id,tensor_index);
#endif

    Index offset = 0;
    for (int comp = 0 ; comp < comp_id ; ++comp)
        offset += index_table_[comp].flat_size();

    const Index local_dof_id = offset + index_table_[comp_id].tensor_to_flat(tensor_index);

    return local_dof_id;
}

template<int dim, int range, int rank>
auto
DofDistribution<dim, range, rank>::
get_dofs_offset() const -> OffsetTable
{
    /*
    OffsetTable offset;
    offset[0] = 0;
    for (int comp = 1; comp < Space::n_components; ++comp)
        offset[comp] = offset[comp-1] + num_dofs_table_.get_component_size(comp-1);
    //*/
    return num_dofs_table_.get_offset();
}

template<int dim, int range, int rank>
void
DofDistribution<dim, range, rank>::
print_info(LogStream &out) const
{
    using std::endl;

    out.begin_item("Dof indices:");
    for (const auto &index_table_comp : index_table_)
    {
        index_table_comp.print_info(out);
        out << endl;
    }
    out.end_item();


    if (!properties_dofs_.empty())
    {
        out.begin_item("Dof properties:");
        properties_dofs_.print_info(out);
        out.end_item();
    }

}


template<int dim, int range, int rank>
bool
DofDistribution<dim, range, rank>::
test_if_dof_has_property(const Index dof_id, const std::string &property) const
{
    return properties_dofs_.test_id_for_property(dof_id, property);
}


template<int dim, int range, int rank>
void
DofDistribution<dim, range, rank>::
add_dofs_property(const std::string &property)
{
    properties_dofs_.add_property(property);
}


template<int dim, int range, int rank>
std::set<Index> &
DofDistribution<dim, range, rank>::
get_dofs_id_same_property(const std::string &property)
{
    return properties_dofs_.get_ids_same_property(property);
}

template<int dim, int range, int rank>
const std::set<Index> &
DofDistribution<dim, range, rank>::
get_dofs_id_same_property(const std::string &property) const
{
    return properties_dofs_.get_ids_same_property(property);
}


template<int dim, int range, int rank>
void
DofDistribution<dim, range, rank>::
set_dof_property_status(const std::string &property, const Index dof_id, const bool status)
{
    properties_dofs_.set_id_property_status(property,dof_id,status);
}


template<int dim, int range, int rank>
void
DofDistribution<dim, range, rank>::
set_all_dofs_property_status(const std::string &property, const bool status)
{
    properties_dofs_.set_ids_property_status(
        property,
        std::set<Index>(dofs_view_.cbegin(),dofs_view_.cend()),
        status);
}


template<int dim, int range, int rank>
auto
DofDistribution<dim, range, rank>::
get_num_dofs_table() const -> const TensorSizeTable &
{
    return num_dofs_table_;
}


IGA_NAMESPACE_CLOSE

#include <igatools/basis_functions/dof_distribution.inst>

