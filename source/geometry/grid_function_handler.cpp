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

template<int dim_, int space_dim_>
GridFunctionHandler<dim_, space_dim_>::
GridFunctionHandler(std::shared_ptr<GridFunctionType> grid_function)
  :
  grid_function_(grid_function),
  grid_handler_(grid_function->get_grid()->create_cache_handler()),
  flags_(CacheFlags::none)
{
  Assert(grid_function_ != nullptr, ExcNullPtr());
}



template<int dim_, int space_dim_>
GridFunctionHandler<dim_, space_dim_>::
~GridFunctionHandler()
{}



template<int dim_, int space_dim_>
auto
GridFunctionHandler<dim_, space_dim_>::
set_flags(const topology_variant &sdim,
          const Flags &flag) -> void
{
  using GridFlags = typename GridType::ElementHandler::Flags;
  GridFlags  grid_flag = GridFlags::none;
  CacheFlags dom_flag = CacheFlags::none;

  SafeSTLVector<Flags> all_flags ={Flags::point, Flags::measure, Flags::w_measure};
  for (auto &fl : all_flags)
    if (contains(flag, fl))
    {
      grid_flag |= grid_function_element::activate::grid[fl];
      dom_flag  |= grid_function_element::activate::grid_function[fl];
    }

  grid_handler_->set_flags(sdim, grid_flag);

  auto disp = SetFlagsDispatcher(dom_flag, flags_);
  boost::apply_visitor(disp, sdim);
}



template<int dim_, int space_dim_>
void
GridFunctionHandler<dim_, space_dim_>::
init_cache(ConstElementAccessor &elem,
           const eval_pts_variant &quad) const
{
  grid_handler_->init_cache(*(elem.grid_elem_), quad);

  auto &cache = elem.local_cache_;
  if (cache == nullptr)
  {
    using Cache = typename ElementAccessor::CacheType;
    cache = std::make_shared<Cache>();
  }

  auto disp = InitCacheDispatcher(this, elem, flags_);
  boost::apply_visitor(disp, quad);

}



template<int dim_, int space_dim_>
auto
GridFunctionHandler<dim_, space_dim_>::
fill_cache(const topology_variant &sdim,
           ConstElementAccessor &elem,
           const int s_id) const-> void
{
  grid_handler_->fill_cache(sdim, *(elem.grid_elem_), s_id);

  auto disp = FillCacheDispatcher(elem, s_id);
  boost::apply_visitor(disp, sdim);

}



IGA_NAMESPACE_CLOSE

#include <igatools/geometry/grid_function_handler.inst>

