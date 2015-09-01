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
Function<dim_, codim_, range_, rank_ >::
Function(std::shared_ptr<const DomainType> phys_dom)
  :
  phys_domain_(phys_dom)
  // object_id_(UniqueIdGenerator::get_unique_id())
{
  Assert(phys_dom != nullptr,ExcNullPtr());
}



template<int dim_, int codim_, int range_, int rank_>
Function<dim_, codim_, range_, rank_ >::
Function(const self_t &func)
  :
  phys_domain_(func.phys_domain_)
{}






template<int dim_, int codim_, int range_, int rank_>
auto
Function<dim_, codim_, range_, rank_ >::
create_cache_handler() const
-> std::shared_ptr<ElementHandler>
{
  return std::make_shared<ElementHandler>(this->shared_from_this());
}



template<int dim_, int codim_, int range_, int rank_>
auto
Function<dim_, codim_, range_, rank_ >::
create_element(const ListIt &index, const PropId &prop) const
-> std::shared_ptr<ConstElementAccessor>
{
  using Elem = ConstElementAccessor;
  auto elem = std::make_shared<Elem>(this->shared_from_this(), index, prop);
  Assert(elem != nullptr,ExcNullPtr());

  return elem;
}

template<int dim_, int codim_, int range_, int rank_>
auto
Function<dim_, codim_, range_, rank_ >::
create_element(const ListIt &index, const PropId &prop)
-> std::shared_ptr<ElementAccessor>
{
  using Elem = ElementAccessor;
  auto elem = std::make_shared<Elem>(this->shared_from_this(), index, prop);
  Assert(elem != nullptr,ExcNullPtr());

  return elem;
}



template<int dim_, int codim_, int range_, int rank_>
auto
Function<dim_, codim_, range_, rank_ >::
begin(const PropId &prop) -> ElementIterator
{
  return ElementIterator(this->shared_from_this(),
  phys_domain_->get_grid()->get_element_property(prop).begin(),
  prop);
}



template<int dim_, int codim_, int range_, int rank_>
auto
Function<dim_, codim_, range_, rank_ >::
end(const PropId &prop) -> ElementIterator
{
  return ElementIterator(this->shared_from_this(),
  phys_domain_->get_grid()->get_element_property(prop).end(),
  prop);
}



template<int dim_, int codim_, int range_, int rank_>
auto
Function<dim_, codim_, range_, rank_ >::
begin(const PropId &prop) const -> ElementConstIterator
{
  return this->cbegin(prop);
}

template<int dim_, int codim_, int range_, int rank_>
auto
Function<dim_, codim_, range_, rank_ >::
end(const PropId &prop) const -> ElementConstIterator
{
  return this->cend(prop);
}

template<int dim_, int codim_, int range_, int rank_>
auto
Function<dim_, codim_, range_, rank_ >::
cbegin(const PropId &prop) const -> ElementConstIterator
{
  return ElementConstIterator(this->shared_from_this(),
                              phys_domain_->get_grid()->get_element_property(prop).end(),
                              prop);
}


template<int dim_, int codim_, int range_, int rank_>
auto
Function<dim_, codim_, range_, rank_ >::
cend(const PropId &prop) const -> ElementConstIterator
{
  return ElementConstIterator(this->shared_from_this(),
                              phys_domain_->get_grid()->get_element_property(prop).end(),
                              prop);
}



#ifdef SERIALIZATION
template<int dim_, int codim_, int range_, int rank_>
Index
Function<dim_, codim_, range_, rank_ >::
get_object_id() const
{
  return object_id_;
}


template<int dim_, int codim_, int range_, int rank_>
const std::string &
Function<dim_, codim_, range_, rank_ >::
get_name() const
{
  return name_;
}

template<int dim_, int codim_, int range_, int rank_>
void
Function<dim_, codim_, range_, rank_ >::
set_name(const std::string &name)
{
  name_ = name;
}


template<class Archive>
void
Function<dim_, codim_, range_, rank_ >::
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

#include <igatools/functions/function.inst>

