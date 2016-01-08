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

#include <igatools/geometry/grid_function_handler.h>
#include <igatools/geometry/grid_function.h>
#include <igatools/geometry/grid_function_element.h>

IGA_NAMESPACE_OPEN

template<int dim_, int range_>
GridFunctionHandler<dim_, range_>::
GridFunctionHandler(const std::shared_ptr<GridFunctionType> &grid_function)
  :
  grid_function_(grid_function),
  grid_handler_(grid_function->get_grid()->create_cache_handler()),
  flags_(Flags::none)
{
  Assert(grid_function_ != nullptr, ExcNullPtr());
}




template<int dim_, int range_>
auto
GridFunctionHandler<dim_, range_>::
get_grid_function() const -> std::shared_ptr<GridFunctionType>
{
  return grid_function_;
}

template<int dim_, int range_>
auto
GridFunctionHandler<dim_, range_>::
get_grid_handler() const -> const GridHandler &
{
  return *grid_handler_;
}

template<int dim_, int range_>
auto
GridFunctionHandler<dim_, range_>::
get_grid_handler() -> GridHandler &
{
  return *grid_handler_;
}


template<int dim_, int range_>
auto
GridFunctionHandler<dim_, range_>::
get_element_cache(ElementAccessor &elem) const
-> typename ElementAccessor::CacheType &
{
  return  elem.local_cache_;
}


template<int dim_, int range_>
auto
GridFunctionHandler<dim_, range_>::
set_flags(const topology_variant &sdim,
          const Flags &flag) -> void
{
  using GridFlags = typename GridType::ElementHandler::Flags;
  GridFlags  grid_flag = GridFlags::none;
  Flags dom_flag = Flags::none;

  for (auto &fl :  grid_function_element::all_flags)
    if (contains(flag, fl))
    {
      grid_flag |= grid_function_element::activate::grid[fl];
      dom_flag  |= grid_function_element::activate::grid_function[fl];
    }

  grid_handler_->set_flags(sdim, grid_flag);

  auto disp = SetFlagsDispatcher(dom_flag, flags_);
  boost::apply_visitor(disp, sdim);
}


template<int dim_, int range_>
template <int sdim>
void
GridFunctionHandler<dim_, range_>::
set_flags(const Flags &flag)
{
  this->set_flags(Topology<sdim>(), flag);
}


template<int dim_, int range_>
void
GridFunctionHandler<dim_, range_>::
set_element_flags(const Flags &flag)
{
  this->set_flags(Topology<dim_>(), flag);
}



template<int dim_, int range_>
void
GridFunctionHandler<dim_, range_>::
init_cache(ElementAccessor &elem,
           const eval_pts_variant &quad) const
{
  grid_handler_->init_cache(*(elem.grid_elem_), quad);

  auto disp = InitCacheDispatcher(*this, elem, flags_);
  boost::apply_visitor(disp, quad);
}

template<int dim_, int range_>
void
GridFunctionHandler<dim_, range_>::
init_cache(ElementIterator &elem,
           const eval_pts_variant &quad) const
{
  this->init_cache(*elem, quad);
}

template<int dim_, int range_>
void
GridFunctionHandler<dim_, range_>::
init_element_cache(ElementAccessor &elem,
                   const std::shared_ptr<const Quadrature<dim_>> &quad) const
{
  this->init_cache(elem,quad);
}

template<int dim_, int range_>
void
GridFunctionHandler<dim_, range_>::
init_element_cache(ElementIterator &elem,
                   const std::shared_ptr<const Quadrature<dim_>> &quad) const
{
  this->init_cache(elem,quad);
}


#if 0
template<int dim_, int range_>
void
GridFunctionHandler<dim_, range_>::
fill_cache(const topology_variant &sdim,
           ElementAccessor &elem,
           const int s_id) const
{
  grid_handler_->fill_cache(sdim, *(elem.grid_elem_), s_id);
}
#endif

template<int dim_, int range_>
void
GridFunctionHandler<dim_, range_>::
fill_cache(const topology_variant &sdim,
           ElementIterator &elem,
           const int s_id) const
{
  this->fill_cache(sdim, *elem, s_id);
}

template<int dim_, int range_>
template <int sdim>
void
GridFunctionHandler<dim_, range_>::
fill_cache(ElementIterator &elem,
           const int s_id) const
{
  this->fill_cache(Topology<sdim>(), elem, s_id);
}

template<int dim_, int range_>
template <int sdim>
void
GridFunctionHandler<dim_, range_>::
fill_cache(ElementAccessor &elem,
           const int s_id) const
{
  this->fill_cache(Topology<sdim>(), elem, s_id);
}


template<int dim_, int range_>
void
GridFunctionHandler<dim_, range_>::
fill_element_cache(ElementAccessor &elem) const
{
  this->fill_cache(Topology<dim_>(), elem,0);
}

template<int dim_, int range_>
void
GridFunctionHandler<dim_, range_>::
fill_element_cache(ElementIterator &elem) const
{
  this->fill_cache(Topology<dim_>(), elem,0);
}


template<int dim_, int range_>
GridFunctionHandler<dim_, range_>::
SetFlagsDispatcher::
SetFlagsDispatcher(const Flags flag, FlagsArray &flags)
  :
  flag_(flag),
  flags_(flags)
{}


template<int dim_, int range_>
template<int sdim>
void
GridFunctionHandler<dim_, range_>::
SetFlagsDispatcher::
operator()(const Topology<sdim> &)
{
  flags_[sdim] |= flag_;
}


template<int dim_, int range_>
GridFunctionHandler<dim_, range_>::
InitCacheDispatcher::
InitCacheDispatcher(const self_t &grid_function_handler,
                    ElementAccessor &elem,
                    const FlagsArray &flags)
  :
  grid_function_handler_(grid_function_handler),
  elem_(elem),
  flags_(flags)
{}

template<int dim_, int range_>
template<int sdim>
void
GridFunctionHandler<dim_, range_>::
InitCacheDispatcher::
operator()(const std::shared_ptr<const Quadrature<sdim>> &quad)
{
  auto &cache = grid_function_handler_.get_element_cache(elem_);

  const auto n_points = elem_.get_grid_element().template get_quad<sdim>()
                        ->get_num_points();
  for (auto &s_id: UnitElement<dim_>::template elems_ids<sdim>())
  {
    auto &s_cache = cache.template get_sub_elem_cache<sdim>(s_id);
    s_cache.resize(flags_[sdim], n_points);
  }
}


IGA_NAMESPACE_CLOSE

#include <igatools/geometry/grid_function_handler.inst>
