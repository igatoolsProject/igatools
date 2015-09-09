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

#include <igatools/geometry/grid_function_element.h>

IGA_NAMESPACE_OPEN

template<int dim_, int space_dim_, class ContainerType_>
GridFunctionElementBase<dim_, space_dim_, ContainerType_>::
GridFunctionElementBase(std::shared_ptr<ContainerType_> grid_function,
                        const ListIt &index,
                        const PropId &prop)
  :
  grid_function_(grid_function),
  grid_elem_(grid_function_->get_grid()->create_element(index,prop))
{}



template<int dim_, int space_dim_, class ContainerType_>
bool
GridFunctionElementBase<dim_, space_dim_, ContainerType_>::
operator ==(const self_t &elem) const
{
  Assert(grid_function_ == elem.grid_function_,
         ExcMessage("Cannot compare elements on different grid."));
  return (*grid_elem_ == *(elem.grid_elem_));
}



template<int dim_, int space_dim_, class ContainerType_>
bool
GridFunctionElementBase<dim_, space_dim_, ContainerType_>::
operator !=(const self_t &elem) const
{
  Assert(grid_function_ == elem.grid_function_,
         ExcMessage("Cannot compare elements on different grid."));
  return (*grid_elem_ != *(elem.grid_elem_));
}



template<int dim_, int space_dim_, class ContainerType_>
bool
GridFunctionElementBase<dim_, space_dim_, ContainerType_>::
operator <(const self_t &elem) const
{
  Assert(grid_function_ == elem.grid_function_,
         ExcMessage("Cannot compare elements on different grid."));
  return (*grid_elem_ < *(elem.grid_elem_));
}



template<int dim_, int space_dim_, class ContainerType_>
bool
GridFunctionElementBase<dim_, space_dim_, ContainerType_>::
operator >(const self_t &elem) const
{
  Assert(grid_function_ == elem.grid_function_,
         ExcMessage("Cannot compare elements on different grid."));
  return (*grid_elem_ > *(elem.grid_elem_));
}


//template<int dim_, int space_dim_, class ContainerType_>
//template<int sdim>
//auto
//GridFunctionElementBase<dim_, space_dim_, ContainerType_>::
//get_w_measures(const int s_id) const -> ValueVector<Real>
//{
//  const auto &meas = get_values_from_cache<_Measure, sdim>(s_id);
//  const auto &w = grid_elem_->template get_weights<sdim>(s_id);
//  auto w_meas = meas;
//  auto it_w = w.begin();
//  for (auto &w_m : w_meas)
//    w_m *= *(it_w);
//  return w_meas;
//}

#if 0
template<int dim_, int space_dim_, class ContainerType_>
auto
GridFunctionElementBase<dim_, space_dim_, ContainerType_>::
get_exterior_normals() const -> ValueVector<SafeSTLArray<Value, space_dim_> >
{
  const int sdim = dim_;
  const int s_id = 0;
  Assert(space_dim_ == 1, ExcNotImplemented());
  ValueVector<SafeSTLArray<Value, space_dim_>> res;
  const auto &DF = this->template get_values_from_cache<_Gradient, sdim>(s_id);
  const auto n_points = DF.get_num_points();
  res.resize(n_points);

  for (int pt = 0; pt < n_points; ++pt)
  {
    res[0][pt] = cross_product<dim_, space_dim_>(DF[pt]);
    res[0][pt] /= res[0][pt].norm();
  }

  return res;
}
#endif

IGA_NAMESPACE_CLOSE

#include <igatools/geometry/grid_function_element.inst>

