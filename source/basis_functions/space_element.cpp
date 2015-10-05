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



#include <igatools/basis_functions/space_element.h>
#include <igatools/basis_functions/reference_element.h>
#include <igatools/basis_functions/physical_space_element.h>


IGA_NAMESPACE_OPEN


template<int dim_,int codim_,int range_,int rank_,Transformation type_>
SpaceElement<dim_,codim_,range_,rank_,type_>::
SpaceElement(const std::shared_ptr<Sp> space,
             const ListIt &index,
             const PropId &prop)
  :
  Element(prop),
  space_(space)
{
//  grid_elem_ = space_->get_ptr_const_grid()->create_element(index,prop);
}



#if 0
template<int dim_,int codim_,int range_,int rank_,Transformation type_>
auto
SpaceElement<dim_,codim_,range_,rank_,type_>::
get_grid_element() -> GridElem &
{
  return *grid_elem_;
}

template<int dim_,int codim_,int range_,int rank_,Transformation type_>
auto
SpaceElement<dim_,codim_,range_,rank_,type_>::
get_grid_element() const -> const  GridElem &
{
  return *grid_elem_;
}
#endif

template<int dim_,int codim_,int range_,int rank_,Transformation type_>
void
SpaceElement<dim_,codim_,range_,rank_,type_>::
print_info(LogStream &out) const
{
  this->get_grid_element().print_info(out);

  out.begin_item("Element global connectivity (property=\"" + DofProperties::active + "\"):");
  const auto glob_dofs = this->get_local_to_global(DofProperties::active);
  glob_dofs.print_info(out);
  out.end_item();
}


template<int dim_,int codim_,int range_,int rank_,Transformation type_>
void
SpaceElement<dim_,codim_,range_,rank_,type_>::
print_cache_info(LogStream &out) const
{
  out.begin_item("GridElement<" + std::to_string(dim) + "> cache:");
  this->get_grid_element().print_cache_info(out);
//  grid_elem_->print_cache_info(out);
  out.end_item();

  //    Assert(all_sub_elems_cache_ != nullptr,ExcNullPtr());
  all_sub_elems_cache_.print_info(out);
}



template<int dim_,int codim_,int range_,int rank_,Transformation type_>
auto
SpaceElement<dim_,codim_,range_,rank_,type_>::
get_index() const -> IndexType
{
  return this->get_grid_element().get_index();
}

template<int dim_,int codim_,int range_,int rank_,Transformation type_>
std::shared_ptr<const Grid<dim_> >
SpaceElement<dim_,codim_,range_,rank_,type_>::
get_grid() const
{
  return this->get_grid_element().get_grid();
}




template<int dim_,int codim_,int range_,int rank_,Transformation type_>
SafeSTLVector<Index>
SpaceElement<dim_,codim_,range_,rank_,type_>::
get_local_to_global(const std::string &dofs_property) const
{
  SafeSTLVector<Index> dofs_global;
  SafeSTLVector<Index> dofs_loc_to_patch;
  SafeSTLVector<Index> dofs_loc_to_elem;
  this->space_->get_element_dofs(
    this->get_index(),
    dofs_global,
    dofs_loc_to_patch,
    dofs_loc_to_elem,
    dofs_property);

  return dofs_global;
}

template<int dim_,int codim_,int range_,int rank_,Transformation type_>
SafeSTLVector<Index>
SpaceElement<dim_,codim_,range_,rank_,type_>::
get_local_to_patch(const std::string &dofs_property) const
{
  SafeSTLVector<Index> dofs_global;
  SafeSTLVector<Index> dofs_loc_to_patch;
  SafeSTLVector<Index> dofs_loc_to_elem;
  this->space_->get_element_dofs(
    this->get_index(),
    dofs_global,
    dofs_loc_to_patch,
    dofs_loc_to_elem,
    dofs_property);

  return dofs_loc_to_patch;
}

template<int dim_,int codim_,int range_,int rank_,Transformation type_>
SafeSTLVector<Index>
SpaceElement<dim_,codim_,range_,rank_,type_>::
get_local_dofs(const std::string &dofs_property) const
{
  SafeSTLVector<Index> dofs_global;
  SafeSTLVector<Index> dofs_loc_to_patch;
  SafeSTLVector<Index> dofs_loc_to_elem;
  this->space_->get_element_dofs(
    this->get_index(),
    dofs_global,
    dofs_loc_to_patch,
    dofs_loc_to_elem,
    dofs_property);

  return dofs_loc_to_elem;
}

template<int dim_,int codim_,int range_,int rank_,Transformation type_>
Size
SpaceElement<dim_,codim_,range_,rank_,type_>::
get_num_basis(const std::string &dofs_property) const
{
  const auto dofs_global = this->get_local_to_global(dofs_property);
  return dofs_global.size();
}


template<int dim_,int codim_,int range_,int rank_,Transformation type_>
bool
SpaceElement<dim_,codim_,range_,rank_,type_>::
operator==(const self_t &a) const
{
  Assert(this->is_comparable_with(a),
         ExcMessage("Comparison between elements defined on different spaces"));
  return this->get_grid_element() == a.get_grid_element();
}

template<int dim_,int codim_,int range_,int rank_,Transformation type_>
bool
SpaceElement<dim_,codim_,range_,rank_,type_>::
operator!=(const self_t &a) const
{
  Assert(this->is_comparable_with(a),
         ExcMessage("Comparison between elements defined on different spaces"));
  return this->get_grid_element() != a.get_grid_element();
}

template<int dim_,int codim_,int range_,int rank_,Transformation type_>
bool
SpaceElement<dim_,codim_,range_,rank_,type_>::
operator<(const self_t &a) const
{
  Assert(this->is_comparable_with(a),
         ExcMessage("Comparison between elements defined on different spaces"));
//  return *grid_elem_ < *a.grid_elem_;
  return this->get_grid_element() < a.get_grid_element();
}

template<int dim_,int codim_,int range_,int rank_,Transformation type_>
bool
SpaceElement<dim_,codim_,range_,rank_,type_>::
operator>(const self_t &a) const
{
  Assert(this->is_comparable_with(a),
         ExcMessage("Comparison between elements defined on different spaces"));
//  return *grid_elem_ > *a.grid_elem_;
  return this->get_grid_element() > a.get_grid_element();
}





template<int dim_,int codim_,int range_,int rank_,Transformation type_>
template <int k>
ValueVector<Real>
SpaceElement<dim_,codim_,range_,rank_,type_>::
get_w_measures(const int j) const
{
  ValueVector<Real> w_measures;

  using RefElem = const ReferenceElement<dim_,range_,rank_>;
  RefElem *as_ref_elem = dynamic_cast<RefElem *>(this);
  if (as_ref_elem)
    w_measures = as_ref_elem->template get_w_measures<k>(j);

  using PhysElem = const PhysicalSpaceElement<dim_,range_,rank_,codim_>;
  PhysElem *as_phys_elem = dynamic_cast<PhysElem *>(this);
  if (as_phys_elem)
    w_measures = as_phys_elem->template get_w_measures<k>(j);

  return w_measures;
}

template<int dim_,int codim_,int range_,int rank_,Transformation type_>
auto
SpaceElement<dim_,codim_,range_,rank_,type_>::
get_space() const -> std::shared_ptr<Sp>
{
  return space_;
}


template<int dim_,int codim_,int range_,int rank_,Transformation type_>
bool
SpaceElement<dim_,codim_,range_,rank_,type_>::
is_comparable_with(const self_t &elem) const
{
	return (space_ == elem.space_);
}

IGA_NAMESPACE_CLOSE


#include <igatools/basis_functions/space_element.inst>


