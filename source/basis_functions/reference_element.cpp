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


#include <igatools/basis_functions/reference_element.h>
#include <igatools/basis_functions/reference_space.h>
#include <igatools/basis_functions/space_element.h>
#include <igatools/basis_functions/reference_element_handler.h>

IGA_NAMESPACE_OPEN



template <int dim, int range, int rank>
ReferenceElement<dim, range, rank>::
ReferenceElement(const std::shared_ptr<ConstSpace> space,
                 const ListIt &index,
                 const PropId &prop)
  :
  parent_t(space,prop),
  space_(space)
{
//    Assert(this->get_space() != nullptr,ExcNullPtr());

  //-------------------------------------------------
  const auto &degree_table = space->get_degree_table();
  TensorSizeTable n_basis(degree_table.get_comp_map());
  for (auto comp : degree_table.get_active_components_id())
    n_basis[comp] = TensorSize<dim>(degree_table[comp]+1);

  n_basis_direction_ = n_basis;
  //-------------------------------------------------


  //----------------------------------------------------------------
  comp_offset_[0] = 0;
  for (int comp = 1; comp <= Space::n_components; ++comp)
    comp_offset_[comp] = comp_offset_[comp-1] +
                         this->n_basis_direction_.get_component_size(comp-1);
  //----------------------------------------------------------------


  //----------------------------------------------------------------
  for (int comp : basis_functions_indexer_.get_active_components_id())
  {
    // creating the objects for fast conversion from flat-to-tensor indexing
    // (in practice it is an hash-table from flat to tensor indices)
    basis_functions_indexer_[comp] =
      std::shared_ptr<Indexer>(new Indexer(this->n_basis_direction_[comp]));
  }
  //----------------------------------------------------------------
};




template <int dim, int range, int rank>
int
ReferenceElement<dim, range, rank>::
get_num_basis_comp(const int i) const
{
  return this->n_basis_direction_[i].flat_size();
}



template <int dim, int range, int rank>
auto
ReferenceElement<dim, range, rank>::
get_basis_offset() const -> OffsetTable
{
  return this->comp_offset_;
}




template <int dim, int range, int rank>
void
ReferenceElement<dim, range, rank>::
print_info(LogStream &out) const
{
  parent_t::print_info(out);


  out.begin_item("Number of element basis: ");
  n_basis_direction_.print_info(out);
  out.end_item();
}





template <int dim, int range, int rank>
auto
ReferenceElement<dim, range, rank>::
get_element_w_measures() const -> ValueVector<Real>
{
  return this->template get_w_measures<dim>(0);
}
#if 0
template <int dim, int range, int rank>
auto
ReferenceElement<dim, range, rank>::
get_ig_space() const -> std::shared_ptr<const Space>
{
  return space_;
}
#endif

IGA_NAMESPACE_CLOSE

#include <igatools/basis_functions/reference_element.inst>


