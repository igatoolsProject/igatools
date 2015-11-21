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

#include <igatools/geometry/grid_handler.h>

using std::shared_ptr;

IGA_NAMESPACE_OPEN


template <int dim>
GridHandler<dim>::
GridHandler(const shared_ptr<GridType> &grid)
  :
  grid_(grid),
  flags_(Flags::none)
{}


template <int dim>
template<int sdim>
void
GridHandler<dim>::
set_flags(const Flags &flag)
{
  Flags grid_flag = Flags::none;
  SafeSTLVector<Flags> all_flags = {Flags::point, Flags::weight};
  for (auto &fl : all_flags)
    if (contains(flag, fl))
    {
      grid_flag |= grid_element::activate::grid[fl];
    }
  flags_[sdim] |= grid_flag;
}



template <int dim>
void
GridHandler<dim>::
set_flags(const topology_variant &sdim,
          const Flags &flag)
{
  auto disp = SetFlagsDispatcher(flag, *this);
  boost::apply_visitor(disp, sdim);
}


template <int dim>
void
GridHandler<dim>::
set_element_flags(const Flags &flag)
{
  this->template set_flags<dim>(flag);
}




template <int dim>
template <int sdim>
void
GridHandler<dim>::
init_cache(ElementAccessor &elem,
           std::shared_ptr<const Quadrature<sdim>> quad) const
{
  Assert(quad != nullptr, ExcNullPtr());

  elem.quad_list_.template get_quad<sdim>() = quad;

  auto &cache = elem.all_sub_elems_cache_;
  /*
  if (cache == nullptr)
  {
    using Cache = typename ElementAccessor::CacheType;
    cache = std::make_unique<Cache>();
  }
  //*/

  for (auto &s_id: UnitElement<dim>::template elems_ids<sdim>())
  {
    auto &s_cache = cache.template get_sub_elem_cache<sdim>(s_id);
    s_cache.resize(flags_[sdim], quad->get_num_points());

  }
}

template <int dim>
void
GridHandler<dim>::
init_cache(ElementAccessor &elem,
           const eval_pts_variant &quad) const
{
  auto disp = InitCacheDispatcher(*this, elem);
  boost::apply_visitor(disp, quad);
}


template <int dim>
void
GridHandler<dim>::
init_element_cache(ElementIterator &elem,
                   std::shared_ptr<const Quadrature<dim>> quad) const
{
  init_cache<dim>(*elem, quad);
}

template <int dim>
void
GridHandler<dim>::
init_face_cache(ElementIterator &elem,
                std::shared_ptr<const Quadrature<(dim > 0) ? dim-1 : 0>> quad) const
{
  Assert(dim > 0,ExcMessage("No face defined for element with topological dimension 0."));
  init_cache<(dim > 0) ? dim-1 : 0>(*elem, quad);
}


template <int dim>
void
GridHandler<dim>::
fill_cache(const topology_variant &sdim,
           ElementAccessor &elem,
           const int s_id) const
{
  auto disp = FillCacheDispatcher(*this, elem, s_id);
  boost::apply_visitor(disp, sdim);
}


template <int dim>
template <int sdim>
void
GridHandler<dim>::
fill_cache(ElementAccessor &elem, const int s_id) const
{
  using _Point = typename ElementAccessor::_Point;
  using _Weight = typename ElementAccessor::_Weight;
//  Assert(elem.all_sub_elems_cache_ != nullptr, ExcNullPtr());
  auto &cache = elem.all_sub_elems_cache_.template get_sub_elem_cache<sdim>(s_id);

  const auto &s_quad = elem.quad_list_.template get_quad<sdim>();

  if (cache.template status_fill<_Point>())
  {
    auto quad = extend_sub_elem_quad<sdim,dim>(*s_quad, s_id);

    const auto translate = elem.vertex(0);
    const auto dilate    = elem.template get_side_lengths<dim>(0);
    quad.dilate(dilate);
    quad.translate(translate);
    cache.template get_data<_Point>().fill(quad.get_points());
//    cache.template set_status_filled<_Point>(true);
  }

  if (cache.template status_fill<_Weight>())
  {
    cache.template get_data<_Weight>().fill(
      elem.template get_measure<sdim>(s_id) * s_quad->get_weights());
//    cache.template set_status_filled<_Weight>(true);
  }

  cache.set_filled(true);
}


template <int dim>
void
GridHandler<dim>::
fill_element_cache(ElementIterator &elem) const
{
  fill_cache<dim>(*elem,0);
}


template <int dim>
void
GridHandler<dim>::
fill_face_cache(ElementIterator &elem, const int s_id) const
{
  Assert(dim > 0,ExcMessage("No face defined for element with topological dimension 0."));
  fill_cache<(dim > 0) ? dim-1 : 0>(*elem,s_id);
}


template <int dim>
auto
GridHandler<dim>::
get_grid() const -> std::shared_ptr<const GridType>
{
  return grid_;
}



template <int dim>
void
GridHandler<dim>::
print_info(LogStream &out) const
{
  out.begin_item("Flags for each subdimension");
  // flags_.print_info(out);
  out.end_item();

}


IGA_NAMESPACE_CLOSE

#include <igatools/geometry/grid_handler.inst>
