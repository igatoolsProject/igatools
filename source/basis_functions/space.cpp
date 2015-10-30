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

#include <igatools/basis_functions/space.h>
#include <igatools/utils/unique_id_generator.h>
#include <igatools/basis_functions/space_element.h>
#include <igatools/functions/ig_function.h>
#include <igatools/functions/identity_function.h>
#include <igatools/basis_functions/dof_distribution.h>
#include <igatools/basis_functions/nurbs_space.h>



using std::shared_ptr;
using std::unique_ptr;

IGA_NAMESPACE_OPEN







template <int dim_,int codim_,int range_,int rank_>
Space<dim_,codim_,range_,rank_>::
Space()
  :
  object_id_(UniqueIdGenerator::get_unique_id())
{}



template <int dim_,int codim_,int range_,int rank_>
Index
Space<dim_,codim_,range_,rank_>::
get_object_id() const
{
  return object_id_;
}


template <int dim_,int codim_,int range_,int rank_>
const std::string &
Space<dim_,codim_,range_,rank_>::
get_name() const
{
  return name_;
}

template <int dim_,int codim_,int range_,int rank_>
void
Space<dim_,codim_,range_,rank_>::
set_name(const std::string &name)
{
  name_ = name;
}

template <int dim_,int codim_,int range_,int rank_>
auto
Space<dim_,codim_,range_,rank_>::
begin(const PropId &prop) const -> ElementIterator
{
  return this->cbegin(prop);
}



template <int dim_,int codim_,int range_,int rank_>
auto
Space<dim_,codim_,range_,rank_>::
end(const PropId &prop) const -> ElementIterator
{
  return this->cend(prop);
}

template <int dim_,int codim_,int range_,int rank_>
auto
Space<dim_,codim_,range_,rank_>::
cbegin(const PropId &prop) const -> ElementIterator
{
  return ElementIterator(
           this->create_element(
             this->get_grid()->get_elements_with_property(prop).begin(),prop));
}



template <int dim_,int codim_,int range_,int rank_>
auto
Space<dim_,codim_,range_,rank_>::
cend(const PropId &prop) const -> ElementIterator
{
  return ElementIterator(
           this->create_element(
             this->get_grid()->get_elements_with_property(prop).end(),prop));
}


template <int dim_,int codim_,int range_,int rank_>
auto
Space<dim_,codim_,range_,rank_>::
get_num_basis() const -> Size
{
  return this->get_ptr_const_dof_distribution()->get_num_dofs_table().total_dimension();
}


template <int dim_,int codim_,int range_,int rank_>
auto
Space<dim_,codim_,range_,rank_>::
get_num_basis(const int comp) const -> Size
{
  return this->get_ptr_const_dof_distribution()->get_num_dofs_table().get_component_size(comp);
}

template <int dim_,int codim_,int range_,int rank_>
auto
Space<dim_,codim_,range_,rank_>::
get_num_basis(const int comp, const int dir) const -> Size
{
  return this->get_ptr_const_dof_distribution()->get_num_dofs_table()[comp][dir];
}


template <int dim_,int codim_,int range_,int rank_>
auto
Space<dim_,codim_,range_,rank_>::
get_elem_num_basis() const -> Size
{
  return this->get_ptr_const_dof_distribution()->get_num_dofs_table().total_dimension();
}

template <int dim_,int codim_,int range_,int rank_>
auto
Space<dim_,codim_,range_,rank_>::
get_global_dof_id(const TensorIndex<dim> &tensor_index,
                  const Index comp) const -> Index
{
  return this->get_ptr_const_dof_distribution()->get_index_table()[comp](tensor_index);
}

template <int dim_,int codim_,int range_,int rank_>
auto
Space<dim_,codim_,range_,rank_>::
get_interior_dofs() const -> std::set<Index>
{
  return this->get_ptr_const_dof_distribution()->get_interior_dofs();
}

template <int dim_,int codim_,int range_,int rank_>
auto
Space<dim_,codim_,range_,rank_>::
get_boundary_dofs(const int s_id, const topology_variant &topology) const -> std::set<Index>
{
  return this->get_ptr_const_dof_distribution()->get_boundary_dofs(s_id,topology);
}

#ifdef MESH_REFINEMENT

#if 0
template <int dim_,int codim_,int range_,int rank_>
void
Space<dim_,codim_,range_,rank_>::
refine_h(const Size n_subdivisions)
{
  this->get_ptr_grid()->refine(n_subdivisions);
}


template <int dim_,int codim_,int range_,int rank_>
void
Space<dim_,codim_,range_,rank_>::
create_connection_for_insert_knots(const std::shared_ptr<self_t> &space)
{
  Assert(space != nullptr, ExcNullPtr());
  Assert(&(*space) == &(*this), ExcMessage("Different objects."));

  auto func_to_connect =
    std::bind(&self_t::rebuild_after_insert_knots,
              space.get(),
              std::placeholders::_1,
              std::placeholders::_2);

  using SlotType = typename Grid<dim>::SignalInsertKnotsSlot;
  this->get_ptr_grid()->connect_insert_knots(SlotType(func_to_connect).track_foreign(space));
}
#endif

#endif // MESH_REFINEMENT




#ifdef SERIALIZATION

template <int dim_,int codim_,int range_,int rank_>
template<class Archive>
void
Space<dim_,codim_,range_,rank_>::
serialize(Archive &ar)
{
  ar &make_nvp("object_id_",object_id_);

  ar &make_nvp("name_",name_);
}
#endif // SERIALIZATION

IGA_NAMESPACE_CLOSE


#include <igatools/basis_functions/space.inst>
