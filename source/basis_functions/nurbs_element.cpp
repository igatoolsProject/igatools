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


#include <igatools/basis_functions/nurbs_element.h>
#include <igatools/basis_functions/nurbs_space.h>
#include <igatools/basis_functions/nurbs_element_handler.h>

#ifdef NURBS

IGA_NAMESPACE_OPEN


template <int dim, int range, int rank>
NURBSElement<dim, range, rank>::
NURBSElement(const std::shared_ptr<ContainerType> space,
             const ListIt &index,
             const PropId &prop)
  :
  parent_t(space,index,prop),
  bspline_elem_(space->get_spline_space(),index,prop),
  weight_elem_(space->weight_func_->create_element(index,prop))
{
//    weight_elem_ =
//        std::shared_ptr<WeightElem>(new WeightElem(space->weight_func_,index));
}




#if 0
template <int dim, int range, int rank>
bool
NURBSElement<dim, range, rank>::
jump(const TensorIndex<dim> &increment)
{
  const bool    grid_elem_active =     parent_t::jump(increment);
  const bool bspline_elem_active = bspline_elem_.jump(increment);

  bool  weight_elems_active = true;
  for (const auto &comp_id : weight_elem_table_.get_active_components_id())
    weight_elems_active = weight_elems_active && weight_elem_table_[comp_id]->jump(increment);

  return grid_elem_active && bspline_elem_active && weight_elems_active;
}
#endif

#if 0
template <int dim, int range, int rank>
void
NURBSElement<dim, range, rank>::
move_to(const Index flat_index)
{
  parent_t::move_to(flat_index);
  bspline_elem_.move_to(flat_index);

  weight_elem_.move_to(flat_index);
}
#endif

template <int dim, int range, int rank>
auto
NURBSElement<dim, range, rank>::
get_nurbs_space() const -> std::shared_ptr<const Space>
{
  const auto nrb_space = std::dynamic_pointer_cast<const Space>(this->space_);
  Assert(nrb_space != nullptr,ExcNullPtr());
  return nrb_space;
}


template <int dim, int range, int rank>
void
NURBSElement<dim, range, rank>::
operator++()
{
  ++bspline_elem_;
  ++(*weight_elem_);
}

template <int dim, int range, int rank>
void
NURBSElement<dim, range, rank>::
move_to(const IndexType &elem_id)
{
  bspline_elem_.move_to(elem_id);
  weight_elem_->move_to(elem_id);
}


template <int dim, int range, int rank>
auto
NURBSElement<dim, range, rank>::
get_grid_element() -> GridElem &
{
  return bspline_elem_.get_grid_element();
}

template <int dim, int range, int rank>
auto
NURBSElement<dim, range, rank>::
get_grid_element() const -> const GridElem &
{
  return bspline_elem_.get_grid_element();
}













IGA_NAMESPACE_CLOSE

#include <igatools/basis_functions/nurbs_element.inst>

#endif // #ifdef NURBS



