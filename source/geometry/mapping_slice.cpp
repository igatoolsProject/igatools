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
#if 0
#include <igatools/base/sub_function.h>
#include <igatools/base/quadrature_lib.h>
#include <igatools/base/exceptions.h>



using std::shared_ptr;
using std::endl;

IGA_NAMESPACE_OPEN

template<int dim_, int codim_>
MappingSlice<dim_, codim_>::
MappingSlice(const std::shared_ptr<const SupMap> map,
             const int face_id,
             const std::shared_ptr<GridType > grid,
             const std::shared_ptr<typename SupMap::GridType::FaceGridMap> elem_map)
  :
  base_t::Mapping(grid),
  map_(map),
  direction_(UnitElement<dim + 1>::face_constant_direction[face_id]),
  value_(UnitElement<dim + 1>::face_side[face_id]),
  element(map_->begin()),
  elem_map_ {elem_map}
{}



template<int dim_, int codim_>
MappingSlice<dim_, codim_>::
MappingSlice(const self_t &map_slice)
  :
  base_t::Mapping(map_slice),
  map_(map_slice.map_),
  direction_(map_slice.direction_),
  value_(map_slice.value_),
  element(map_slice.element)
{}



template<int dim_, int codim_>
auto
MappingSlice<dim_, codim_>::
create(const std::shared_ptr<const SupMap> map,
       const int face_id,
       const std::shared_ptr<GridType > grid,
       const std::shared_ptr<typename SupMap::GridType::FaceGridMap> elem_map) -> shared_ptr<base_t>
{
  return shared_ptr<base_t>(new self_t(map, face_id, grid, elem_map));
}



template<int dim_, int codim_>
auto
MappingSlice<dim_, codim_>::
build_extended_quadrature(const Quadrature<dim> &quad) const -> Quadrature<dim+1>
{
  const auto points  = quad.get_points();
  const auto weights = quad.get_weights();

  auto ext_quad = Quadrature<dim+1>(
                    insert(points, direction_,SafeSTLVector<Real>(1,value_)),
                    insert(weights,direction_,SafeSTLVector<Real>(1,1.0))) ;

  return ext_quad;
}



template<int dim_, int codim_>
void
MappingSlice<dim_, codim_>::
evaluate(ValueVector<Value> &values) const
{
  values = element->get_map_values();
}



template<int dim_, int codim_>
void
MappingSlice<dim_, codim_>::
evaluate_gradients(ValueVector<Gradient> &gradients) const
{
  auto grad = element->get_map_gradients();

  const auto active_dir = UnitElement<dim+1>::active_directions[direction_];

  const int num_points = grad.size();
  for (int p = 0 ; p < num_points; p++)
    for (int i = 0; i < dim; i++)
      gradients[p][i] = grad[p][active_dir[i]];

}



template<int dim_, int codim_>
void
MappingSlice<dim_, codim_>::
init_element(const ValueFlags flag, const Quadrature<dim> &quad) const
{
  element->init_cache(flag, build_extended_quadrature(quad));
}



template<int dim_, int codim_>
void
MappingSlice<dim_, codim_>::
set_element(const GridIterator &elem) const
{
  typename Grid<dim_>::ElementIterator el_tmp(this->get_grid(),elem.get_flat_index());
  element->move_to((*elem_map_)[el_tmp]->get_flat_index());
  element->fill_cache();
}


template<int dim_, int codim_>
void
MappingSlice<dim_,codim_>::
set_face_element(const Index face_id,
                 const GridIterator &elem) const
{
  Assert(false, ExcNotImplemented());
  AssertThrow(false, ExcNotImplemented());
}


template<int dim_, int codim_>
auto
MappingSlice<dim_,codim_>::
inject_points(const ValueVector<Point> &points) const -> ValueVector<typename SupMap::Point>
{
  Assert(!points.empty(),ExcEmptyObject());

  const Size n_points = points.size();

  using SupMapPoint = typename SupMap::Point;
  ValueVector<SupMapPoint> sup_map_points(n_points);
  if (dim_ == 1)
  {
    if (direction_ == 0)
    {
      for (Index ipt = 0 ; ipt < n_points ; ++ipt)
      {
        // constant along direction 0
        sup_map_points[ipt][0] = value_;
        sup_map_points[ipt][1] = points[ipt][0];
      }
    }
    else if (direction_ == 1)
    {
      for (Index ipt = 0 ; ipt < n_points ; ++ipt)
      {
        // constant along direction 1
        sup_map_points[ipt][0] = points[ipt][0];
        sup_map_points[ipt][1] = value_;
      }
    }
  }
  else if (dim_ == 2)
  {
    if (direction_ == 0)
    {
      for (Index ipt = 0 ; ipt < n_points ; ++ipt)
      {
        // constant along direction 0
        sup_map_points[ipt][0] = value_;
        sup_map_points[ipt][1] = points[ipt][0];
        sup_map_points[ipt][2] = points[ipt][1];
      }
    }
    else if (direction_ == 1)
    {
      for (Index ipt = 0 ; ipt < n_points ; ++ipt)
      {
        // constant along direction 1
        sup_map_points[ipt][0] = points[ipt][0];
        sup_map_points[ipt][1] = value_;
        sup_map_points[ipt][2] = points[ipt][1];
      }
    }
    else if (direction_ == 2)
    {
      for (Index ipt = 0 ; ipt < n_points ; ++ipt)
      {
        // constant along direction 2
        sup_map_points[ipt][0] = points[ipt][0];
        sup_map_points[ipt][1] = points[ipt][1];
        sup_map_points[ipt][2] = value_;
      }
    }
  }
  else
  {
    Assert(false,ExcNotImplemented());
  }

  return sup_map_points;
}

template<int dim_, int codim_>
void
MappingSlice<dim_,codim_>::
evaluate_at_points(const ValueVector<Point> &points, ValueVector<Value> &values) const
{
  Assert(!points.empty(),ExcEmptyObject());
  Assert(values.size() == points.size(),ExcDimensionMismatch(values.size(),points.size()))

  const auto sup_map_points = this->inject_points(points);
  map_->evaluate_at_points(sup_map_points,values);
}




template<int dim_, int codim_>
void
MappingSlice<dim_,codim_>::
evaluate_gradients_at_points(const ValueVector<Point> &points, ValueVector<Gradient> &gradients) const
{
  Assert(!points.empty(),ExcEmptyObject());
  Assert(gradients.size() == points.size(),ExcDimensionMismatch(gradients.size(),points.size()))

  const Size n_points = points.size();

  const auto sup_map_points = this->inject_points(points);

  using SupMapGrad  = typename SupMap::Gradient;
  ValueVector<SupMapGrad> sup_map_gradients(n_points);
  map_->evaluate_gradients_at_points(sup_map_points,sup_map_gradients);


  using UnitElem = UnitElement<dim_+1>;
  for (int current_dir = 0 ; current_dir < dim_ ; ++current_dir)
  {
    // getting the active direction id form the current one
    const int active_direction_id = UnitElem::active_directions[direction_][current_dir];

    for (Index ipt = 0 ; ipt < n_points ; ++ipt)
      for (int i = 0 ; i <= dim_ ; ++i)
        gradients[ipt][current_dir][i] = sup_map_gradients[ipt][active_direction_id][i];
  }
  Assert(dim_ == 1 || dim_ == 2,ExcNotImplemented());
}

template<int dim_, int codim_>
void
MappingSlice<dim_,codim_>::
evaluate_hessians_at_points(const ValueVector<Point> &points, ValueVector<Hessian> &hessians) const
{
  Assert(false,ExcNotImplemented());
}



//TODO(pauletti, Jun 20, 2014): simplify the output
template<int dim_, int codim_>
void
MappingSlice<dim_, codim_>::
print_info(LogStream &out) const
{
  out << "Type = MappingSlice<" << dim_ << "," << dim_+codim_ << ">" << endl;

  out.push("\t");
  out << "Direction = " << direction_ << endl ;
  out << "    ValueType = " << value_ << endl ;

  out << "Sliced Map:" << endl ;
  out.push("\t");
  map_->print_info(out) ;
  out << endl;

  out.pop();
}

IGA_NAMESPACE_CLOSE

#include <igatools/geometry/mapping_slice.inst>
#endif
