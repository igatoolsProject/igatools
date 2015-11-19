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

#include <igatools/geometry/sub_grid_function_handler.h>
#include <igatools/geometry/sub_grid_function.h>
#include <igatools/geometry/sub_grid_function_element.h>

IGA_NAMESPACE_OPEN

template<int sdim,int dim, int space_dim>
SubGridFunctionHandler<sdim,dim,space_dim>::
SubGridFunctionHandler(const std::shared_ptr<GridFunctionType> &grid_function)
  :
  parent_t(grid_function),
  sup_grid_func_handler_(grid_function->get_sup_grid_function()->create_cache_handler())
{}


template<int sdim,int dim, int space_dim>
void
SubGridFunctionHandler<sdim,dim,space_dim>::
set_flags(const topology_variant &topology,
          const Flags &flag)
{
  parent_t::set_flags(topology,flag);

  sup_grid_func_handler_->set_flags(topology,flag);
}

template<int sdim,int dim, int space_dim>
void
SubGridFunctionHandler<sdim,dim,space_dim>::
init_cache(GridFunctionElement<sdim,space_dim> &sub_grid_func_elem,
           const eval_pts_variant &quad) const
{
  parent_t::init_cache(sub_grid_func_elem,quad);

  auto &as_sub_grid_func_elem =
    dynamic_cast<SubGridFunctionElement<sdim,dim,space_dim> &>(sub_grid_func_elem);

  this->sup_grid_func_handler_->init_cache(as_sub_grid_func_elem.get_sup_grid_function_element(),quad);
}


template<int sdim,int dim, int space_dim>
void
SubGridFunctionHandler<sdim,dim,space_dim>::
fill_cache(const topology_variant &topology,
           GridFunctionElement<sdim,space_dim> &sub_grid_func_elem,
           const int s_id) const
{
  this->grid_handler_->fill_cache(topology, sub_grid_func_elem.get_grid_element(), s_id);

  auto &as_sub_grid_func_elem =
    dynamic_cast<SubGridFunctionElement<sdim,dim,space_dim> &>(sub_grid_func_elem);

  this->sup_grid_func_handler_->fill_cache(topology,as_sub_grid_func_elem.get_sup_grid_function_element(),s_id);

  auto fill_dispatcher = FillCacheDispatcher(
                           *this,
                           as_sub_grid_func_elem,
                           s_id);
  boost::apply_visitor(fill_dispatcher, topology);
}


template<int sdim,int dim, int space_dim>
SubGridFunctionHandler<sdim,dim,space_dim>::
FillCacheDispatcher::
FillCacheDispatcher(const SubGridFunctionHandler<sdim,dim,space_dim> &sub_grid_func_handler,
                    SubGridFunctionElement<sdim,dim,space_dim> &sub_grid_func_elem,
                    const int s_id)
  :
  sub_grid_func_handler_(sub_grid_func_handler),
  sub_grid_func_elem_(sub_grid_func_elem),
  s_id_(s_id)
{}


IGA_NAMESPACE_CLOSE

#include <igatools/geometry/sub_grid_function_handler.inst>
