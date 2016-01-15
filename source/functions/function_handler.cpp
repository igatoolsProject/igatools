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

#include <igatools/functions/function.h>
#include <igatools/functions/function_element.h>
//#include <igatools/utils/unique_id_generator.h>
#include <igatools/geometry/domain.h>
#include <igatools/geometry/domain_element.h>
#include <igatools/geometry/domain_handler.h>

IGA_NAMESPACE_OPEN

template<int dim_, int codim_, int range_, int rank_ >
FunctionHandler<dim_, codim_, range_, rank_ >::
FunctionHandler(std::shared_ptr<FuncType> func)
  :
  func_(func),
  domain_handler_(func_->get_domain()->create_cache_handler())
{}



template<int dim_, int codim_, int range_, int rank_>
void
FunctionHandler<dim_, codim_, range_, rank_ >::
set_flags(const topology_variant &sdim,
          const Flags &flag)
{
  using DomainFlags = typename DomainType::Handler::Flags;
  DomainFlags  dom_flag = DomainFlags::none;
  Flags func_flag = Flags::none;

  for (auto &fl : function_element::all_flags)
    if (contains(flag, fl))
      func_flag  |= function_element::activate::function[fl];

  for (auto &fl : function_element::all_flags)
    if (contains(flag, fl))
      dom_flag  |= function_element::activate::domain[fl];

  domain_handler_->set_flags(sdim, dom_flag);

  auto disp = SetFlagsDispatcher(func_flag, flags_);
  boost::apply_visitor(disp, sdim);
}


template<int dim_, int codim_, int range_, int rank_>
void
FunctionHandler<dim_, codim_, range_, rank_ >::
set_element_flags(const Flags &flag)
{
  this->set_flags(Topology<dim_>(), flag);
}


template<int dim_, int codim_, int range_, int rank_>
void
FunctionHandler<dim_, codim_, range_, rank_ >::
init_cache(ElementAccessor &elem,
           const eval_pts_variant &quad) const
{
  domain_handler_->init_cache(*(elem.domain_elem_), quad);

  auto disp = InitCacheDispatcher(elem, flags_);
  boost::apply_visitor(disp, quad);
}

template<int dim_, int codim_, int range_, int rank_>
void
FunctionHandler<dim_, codim_, range_, rank_ >::
init_cache(ElementIterator &elem,
           const eval_pts_variant &quad) const
{
  this->init_cache(*elem, quad);
}



template<int dim_, int codim_, int range_, int rank_>
void
FunctionHandler<dim_, codim_, range_, rank_ >::
fill_cache(const topology_variant &sdim,
           ElementAccessor &elem,
           const int s_id) const
{
  domain_handler_->fill_cache(sdim, *(elem.domain_elem_), s_id);
}


template<int dim_, int codim_, int range_, int rank_>
void
FunctionHandler<dim_, codim_, range_, rank_ >::
fill_cache(const topology_variant &sdim,
           ElementIterator &elem,
           const int s_id) const
{
  this->fill_cache(sdim, *elem, s_id);
}


template<int dim_, int codim_, int range_, int rank_>
void
FunctionHandler<dim_, codim_, range_, rank_ >::
fill_element_cache(ElementAccessor &elem)
{
  this->fill_cache(Topology<dim_>(),elem,0);
}

template<int dim_, int codim_, int range_, int rank_>
void
FunctionHandler<dim_, codim_, range_, rank_ >::
fill_element_cache(ElementIterator &elem)
{
  this->fill_cache(Topology<dim_>(),*elem,0);
}

IGA_NAMESPACE_CLOSE

#include <igatools/functions/function_handler.inst>

