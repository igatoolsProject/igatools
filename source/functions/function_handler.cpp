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

IGA_NAMESPACE_OPEN

template<int dim_, int codim_, int range_, int rank_ >
FunctionHandler<dim_, codim_, range_, rank_ >::
FunctionHandler(std::shared_ptr<FuncType> func)
  :
  func_(func),
  domain_handler_(func_->get_domain()->create_cache_handler())
  // object_id_(UniqueIdGenerator::get_unique_id())
{}



template<int dim_, int codim_, int range_, int rank_>
void
FunctionHandler<dim_, codim_, range_, rank_ >::
set_flags(const topology_variant &sdim,
          const Flags &flag)
{
  using DomainFlags = typename DomainType::ElementHandler::Flags;
  DomainFlags  dom_flag = DomainFlags::none;
  CacheFlags func_flag = CacheFlags::none;

  SafeSTLVector<Flags> all_flags = {Flags::value, Flags::gradient, Flags::D2};
  for (auto &fl : all_flags)
    if (contains(flag, fl))
    {
      dom_flag  |= function_element::activate::domain[fl];
      func_flag |= function_element::activate::function[fl];
    }

  domain_handler_->set_flags(sdim, dom_flag);

  auto disp = SetFlagsDispatcher(func_flag, flags_);
  boost::apply_visitor(disp, sdim);
}



template<int dim_, int codim_, int range_, int rank_>
void
FunctionHandler<dim_, codim_, range_, rank_ >::
init_cache(ConstElementAccessor &elem,
           const eval_pts_variant &quad) const
{
  domain_handler_->init_cache(*(elem.domain_elem_), quad);

  auto &cache = elem.local_cache_;
  if (cache == nullptr)
  {
    using Cache = typename ElementAccessor::CacheType;
    cache = std::make_shared<Cache>();
  }

  // auto disp = InitCacheDispatcher(this, elem, flags_);
  auto disp = InitCacheDispatcher(elem, flags_);
  boost::apply_visitor(disp, quad);
}



template<int dim_, int codim_, int range_, int rank_>
void
FunctionHandler<dim_, codim_, range_, rank_ >::
fill_cache(const topology_variant &sdim,
           ConstElementAccessor &elem,
           const int s_id) const
{
  domain_handler_->fill_cache(sdim, *(elem.domain_elem_), s_id);
}



#ifdef SERIALIZATION
template<int dim_, int codim_, int range_, int rank_>
Index
FunctionHandler<dim_, codim_, range_, rank_ >::
get_object_id() const
{
  return object_id_;
}


template<int dim_, int codim_, int range_, int rank_>
const std::string &
FunctionHandler<dim_, codim_, range_, rank_ >::
get_name() const
{
  return name_;
}

template<int dim_, int codim_, int range_, int rank_>
void
FunctionHandler<dim_, codim_, range_, rank_ >::
set_name(const std::string &name)
{
  name_ = name;
}


template<class Archive>
void
FunctionHandler<dim_, codim_, range_, rank_ >::
serialize(Archive &ar, const unsigned int version)
{
  ar &boost::serialization::make_nvp("grid_elem_handler_",
                                     boost::serialization::base_object<GridHandler<dim_>>(*this));

  ar &boost::serialization::make_nvp("object_id_",object_id_);
  ar &boost::serialization::make_nvp("name_",name_);

  ar &boost::serialization::make_nvp("flags_",flags_);

  ar &boost::serialization::make_nvp("grid_",grid_);

#ifdef MESH_REFINEMENT
  ar &boost::serialization::make_nvp("function_previous_refinement_",function_previous_refinement_);
#endif // MESH_REFINEMENT
}
#endif // SERIALIZATION


IGA_NAMESPACE_CLOSE

#include <igatools/functions/function_handler.inst>

