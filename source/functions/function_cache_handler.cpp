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
#include <igatools/geometry/physical_domain.h>
#include <igatools/geometry/physical_domain_element.h>

IGA_NAMESPACE_OPEN

template<int dim_, int codim_, int range_, int rank_ >
FunctionElementHandler<dim_, codim_, range_, rank_ >::
FunctionElementHandler(std::shared_ptr<FuncType> func)
  :
  func_(func),
  domain_handler_(func_->get_physical_domain()->create_cache_handler())
  // object_id_(UniqueIdGenerator::get_unique_id())
{
// Assert(phys_dom != nullptr,ExcNullPtr());
}



template<int dim_, int codim_, int range_, int rank_>
void
FunctionElementHandler<dim_, codim_, range_, rank_ >::
set_flags(const topology_variant &sdim,
          const Flags &flag)
{
  auto disp = SetFlagsDispatcher(flag, flags_);
  boost::apply_visitor(disp, sdim);
}



template<int dim_, int codim_, int range_, int rank_>
void
FunctionElementHandler<dim_, codim_, range_, rank_ >::
init_cache(ElementAccessor &elem,
           const eval_pts_variant &quad) const
{
  auto disp = InitCacheDispatcher(flags_, elem);
  boost::apply_visitor(disp, quad);
}



template<int dim_, int codim_, int range_, int rank_>
void
FunctionElementHandler<dim_, codim_, range_, rank_ >::
fill_cache(const topology_variant &sdim,
           ElementAccessor &elem,
           const int s_id) const
{
#if 0
  auto fill_dispatcher = FillCacheDispatcher(s_id, *this, elem);
  boost::apply_visitor(fill_dispatcher, sdim);
#endif
}



#if 0
template<int dim_, int codim_, int range_, int rank_>
void
FunctionElementHandler<dim_, codim_, range_, rank_ >::
init_cache(ElementIterator &elem, const topology_variant &sdim) const
{
  init_cache(*elem, k);
}


template<int dim_, int codim_, int range_, int rank_>
void
FunctionElementHandler<dim_, codim_, range_, rank_ >::
init_element_cache(ElementAccessor &elem) const
{
  this->init_cache(elem, Topology<dim_>());
}


template<int dim_, int codim_, int range_, int rank_>
void
FunctionElementHandler<dim_, codim_, range_, rank_ >::
init_element_cache(ElementIterator &elem) const
{
  this->init_cache(*elem, Topology<dim_>());
}
#endif



#if 0
template<int dim_, int codim_, int range_, int rank_>
void
FunctionElementHandler<dim_, codim_, range_, rank_ >::
fill_cache(ElementIterator &elem, const topology_variant &sdim, const int j) const
{
  this->fill_cache(*elem, k, j);
}


template<int dim_, int codim_, int range_, int rank_>
void
FunctionElementHandler<dim_, codim_, range_, rank_ >::
fill_element_cache(ElementAccessor &elem) const
{
  this->fill_cache(elem, Topology<dim_>(),0);
}


template<int dim_, int codim_, int range_, int rank_>
void
FunctionElementHandler<dim_, codim_, range_, rank_ >::
fill_element_cache(ElementIterator &elem) const
{
  this->fill_cache(*elem, Topology<dim_>(),0);
}
#endif



//template<int dim_, int codim_, int range_, int rank_>
//auto
//FunctionElementHandler<dim_, codim_, range_, rank_ >::
//get_cache(ElementAccessor &elem)
//-> std::shared_ptr<typename ElementAccessor::CacheType> &
//{
//    Assert(elem.all_sub_elems_cache_ != nullptr,ExcNullPtr());
//    return elem.all_sub_elems_cache_;
//}






#ifdef SERIALIZATION
template<int dim_, int codim_, int range_, int rank_>
Index
FunctionElementHandler<dim_, codim_, range_, rank_ >::
get_object_id() const
{
  return object_id_;
}


template<int dim_, int codim_, int range_, int rank_>
const std::string &
FunctionElementHandler<dim_, codim_, range_, rank_ >::
get_name() const
{
  return name_;
}

template<int dim_, int codim_, int range_, int rank_>
void
FunctionElementHandler<dim_, codim_, range_, rank_ >::
set_name(const std::string &name)
{
  name_ = name;
}


template<class Archive>
void
FunctionElementHandler<dim_, codim_, range_, rank_ >::
serialize(Archive &ar, const unsigned int version)
{
  ar &boost::serialization::make_nvp("grid_elem_handler_",
                                     boost::serialization::base_object<GridElementHandler<dim_>>(*this));

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

#include <igatools/functions/function_cache_handler.inst>

