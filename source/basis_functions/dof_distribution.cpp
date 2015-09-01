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
#include <igatools/base/properties.h>
#include <igatools/utils/tensor_range.h>

using std::map;
using std::shared_ptr;
using std::make_shared;

IGA_NAMESPACE_OPEN

template<int dim, int range, int rank>
DofDistribution<dim, range, rank>::
DofDistribution(const TensorSizeTable &n_basis,
                const DegreeTable &degree_table,
                const PeriodicityTable &periodic)
  :
  num_dofs_table_(n_basis),
  index_table_size_(n_basis)
{
  //-----------------------------------------------------------------------
  for (const auto comp : Space::components)
    for (const auto dir : UnitElement<dim>::active_directions)
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
      for (const auto dir : UnitElement<dim>::active_directions)
        t_ind[dir] %= n_dofs_comp[dir];

      const auto f_ind = MultiArrayUtils<dim>::tensor_to_flat_index(t_ind,w_dofs_comp);
      index_table_comp[i] = comp_offset + f_ind;

    }
    comp_offset += num_dofs_table_.get_component_size(comp);//n_basis_comp_size;
  }
  //-----------------------------------------------------------------------


  //-----------------------------------------------------------------------
  properties_dofs_.add_property(DofProperties::active);
  this->set_all_dofs_property_status(DofProperties::active,true);
  //-----------------------------------------------------------------------
}



template<int dim, int range, int rank>
Index
DofDistribution<dim, range, rank>::
get_min_dof_id() const
{
  const auto dofs_view = this->get_dofs_const_view();
  return *std::min_element(dofs_view.cbegin(),dofs_view.cend());
}



template<int dim, int range, int rank>
Index
DofDistribution<dim, range, rank>::
get_max_dof_id() const
{
  const auto dofs_view = this->get_dofs_const_view();
  return *std::max_element(dofs_view.cbegin(),dofs_view.cend());
}



template<int dim, int range, int rank>
bool
DofDistribution<dim, range, rank>::
find_dof_id(const Index dof_id, int &comp_id, Index &dof_id_in_component) const
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

      const auto dofs_offset = this->get_dofs_offset();
      dof_id_in_component = dof_id - dofs_offset[comp_id];
//            tensor_index = index_table_comp.flat_to_tensor(it-dofs_begin);

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
  for (auto &index_table_comp : index_table_)
    for (auto &dof : index_table_comp.get_flat_view())
      dof += offset;


  for (auto &property_dofs : properties_dofs_)
    for (auto &dof : property_dofs.second)
      const_cast<Index &>(dof) += offset;
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
get_dofs_const_view() const -> DofsConstView
{
  // creating the dofs view from the dofs components views
  SafeSTLVector<DofsComponentConstView> components_views;
  for (auto &index_table_comp : index_table_)
    components_views.emplace_back(index_table_comp.get_flat_view());

  return DofsConstView(DofsConstIterator(components_views,0),
                       DofsConstIterator(components_views,IteratorState::pass_the_end));
}

template<int dim, int range, int rank>
auto
DofDistribution<dim, range, rank>::
get_dofs_view() -> DofsView
{
  // creating the dofs view from the dofs components views
  SafeSTLVector<DofsComponentView> components_views;
  for (auto &index_table_comp : index_table_)
    components_views.emplace_back(index_table_comp.get_flat_view());

  return DofsView(DofsIterator(components_views,0),
  DofsIterator(components_views,IteratorState::pass_the_end));
}

template<int dim, int range, int rank>
Index
DofDistribution<dim, range, rank>::
global_to_patch_local(const Index global_dof_id) const
{
  int comp;
  Index dof_id_comp;
  this->global_to_comp_local(global_dof_id,comp,dof_id_comp);

  //TODO (martinelli, 19Mar2015): this is wrong for periodic dofs
  const Index offset = index_table_size_.get_offset()[comp];

  const Index local_dof_id = offset + dof_id_comp;

  return local_dof_id;
}

template<int dim, int range, int rank>
void
DofDistribution<dim, range, rank>::
global_to_comp_local(const Index global_dof_id,int &comp,int &dof_id_comp) const
{
#ifndef NDEBUG
  bool global_dof_is_found = this->find_dof_id(global_dof_id,comp,dof_id_comp);
  Assert(global_dof_is_found,
         ExcMessage("The global dof id " + std::to_string(global_dof_id) +
                    " is not present in the DofDistribution."));
#else
  this->find_dof_id(global_dof_id,comp,dof_id_comp);
#endif
}


template<int dim, int range, int rank>
auto
DofDistribution<dim, range, rank>::
get_dofs_offset() const -> OffsetTable
{
  return num_dofs_table_.get_offset();
}

template<int dim, int range, int rank>
void
DofDistribution<dim, range, rank>::
print_info(LogStream &out) const
{
  using std::endl;

  out.begin_item("Num dofs table:");
  num_dofs_table_.print_info(out);
  out.end_item();

  out.begin_item("Index table size:");
  index_table_size_.print_info(out);
  out.end_item();

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
is_property_defined(const std::string &property) const
{
  return properties_dofs_.is_property_defined(property);
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
Size
DofDistribution<dim, range, rank>::
get_num_dofs(const std::string &property) const
{
  return properties_dofs_[property].size();
}




template<int dim, int range, int rank>
std::set<Index> &
DofDistribution<dim, range, rank>::
get_dofs_id_same_property(const std::string &property)
{
  return properties_dofs_[property];
}



template<int dim, int range, int rank>
const std::set<Index> &
DofDistribution<dim, range, rank>::
get_dofs_id_same_property(const std::string &property) const
{
  return properties_dofs_[property];
}


template<int dim, int range, int rank>
void
DofDistribution<dim, range, rank>::
set_dof_property_status(const std::string &property, const Index dof_id, const bool status)
{
  if (status)
    properties_dofs_[property].insert(dof_id);
  else
    properties_dofs_[property].erase(dof_id);

}



template<int dim, int range, int rank>
void
DofDistribution<dim, range, rank>::
set_dof_property_status(const std::string &property,
                        const std::set<Index> ids,
                        const bool status)
{
  if (status)
    properties_dofs_[property].insert(ids.begin(),ids.end());
  else
    properties_dofs_[property].erase(ids.begin(),ids.end());
}



template<int dim, int range, int rank>
void
DofDistribution<dim, range, rank>::
set_all_dofs_property_status(const std::string &property, const bool status)
{
  const auto dofs_view = this->get_dofs_const_view();
  if (status)
    properties_dofs_[property].insert(dofs_view.cbegin(),dofs_view.cend());
  else
  {
    for (const auto &dof : dofs_view)
      properties_dofs_[property].erase(dof);
  }
}


template<int dim, int range, int rank>
auto
DofDistribution<dim, range, rank>::
get_num_dofs_table() const -> const TensorSizeTable &
{
  return num_dofs_table_;
}




template<int dim, int range, int rank>
auto
DofDistribution<dim, range, rank>::
get_interior_dofs() const -> std::set<Index>
{
  Assert(num_dofs_table_ == index_table_size_,ExcNotImplemented());
  /*
  #ifndef NDEBUG
  for (int comp : end_b_.get_active_components_id())
      for (int j=0; j<dim; ++j)
          Assert(end_b_[comp][j] == BasisEndBehaviour::interpolatory,
          ExcNotImplemented());
  #endif
  //*/

  std::set<Index> dofs;

  TensorIndex<dim> first(1);
  TensorIndex<dim> last;

  for (int comp = 0 ; comp < IndexDistributionTable::n_entries ; ++comp)
  {
    for (int j = 0; j < dim ; ++j)
    {
//            first[j] = 1;
      last[j] = num_dofs_table_[comp][j]-1;
    }

    auto tensor_ind = el_tensor_range<dim>(first, last);
    const auto &elem_global_indices = index_table_[comp];

    for (auto &tensor_index : tensor_ind)
      dofs.insert(elem_global_indices(tensor_index));
  }

  return dofs;
}

#ifdef SERIALIZATION
template<int dim, int range, int rank>
template<class Archive>
void
DofDistribution<dim, range, rank>::
serialize(Archive &ar, const unsigned int version)
{
  ar &boost::serialization::make_nvp("index_table_",index_table_);

  ar &boost::serialization::make_nvp("num_dofs_table_",num_dofs_table_);

  ar &boost::serialization::make_nvp("index_table_size_",index_table_size_);

  ar &boost::serialization::make_nvp("properties_dofs_",properties_dofs_);
}
#endif // SERIALIZATION


IGA_NAMESPACE_CLOSE

#include <igatools/basis_functions/dof_distribution.inst>

