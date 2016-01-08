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
#include <igatools/utils/unique_id_generator.h>
//#include <igatools/geometry/domain.h>
//#include <igatools/geometry/domain_element.h>

IGA_NAMESPACE_OPEN

using std::to_string;

template<int dim_, int codim_, int range_, int rank_ >
Function<dim_, codim_, range_, rank_ >::
Function(const SharedPtrConstnessHandler<DomainType> &domain)
  :
  domain_(domain),
  object_id_(UniqueIdGenerator::get_unique_id())
{

}



//template<int dim_, int codim_, int range_, int rank_>
//Function<dim_, codim_, range_, rank_ >::
//Function(const self_t &func)
//  :
//  domain_(func.domain_)
//{}



template<int dim_, int codim_, int range_, int rank_>
auto
Function<dim_, codim_, range_, rank_ >::
create_cache_handler() const
-> std::unique_ptr<ElementHandler>
{
  return std::unique_ptr<ElementHandler>(new ElementHandler(this->shared_from_this()));
}


#if 0
template<int dim_, int codim_, int range_, int rank_>
auto
Function<dim_, codim_, range_, rank_ >::
create_element(const ListIt &index, const PropId &prop) const
-> std::unique_ptr<ElementAccessor>
{
  using Elem = ElementAccessor;
  return std::unique_ptr<Elem>(new Elem(this->shared_from_this(), index, prop));
}
#endif

template<int dim_, int codim_, int range_, int rank_>
auto
Function<dim_, codim_, range_, rank_ >::
create_element_begin(const PropId &prop) const
-> std::unique_ptr<ElementAccessor>
{
  using Elem = ElementAccessor;
  return std::unique_ptr<Elem>(new Elem(
    this->shared_from_this(),
    domain_->create_element_begin(prop)));
}


template<int dim_, int codim_, int range_, int rank_>
auto
Function<dim_, codim_, range_, rank_ >::
create_element_end(const PropId &prop) const
-> std::unique_ptr<ElementAccessor>
{
  using Elem = ElementAccessor;
  return std::unique_ptr<Elem>(new Elem(
    this->shared_from_this(),
    domain_->create_element_end(prop)));
}


template<int dim_, int codim_, int range_, int rank_>
auto
Function<dim_, codim_, range_, rank_ >::
begin(const PropId &prop) const -> ElementIterator
{
  return this->cbegin(prop);
}



template<int dim_, int codim_, int range_, int rank_>
auto
Function<dim_, codim_, range_, rank_ >::
end(const PropId &prop) const -> ElementIterator
{
  return this->cend(prop);
}



template<int dim_, int codim_, int range_, int rank_>
auto
Function<dim_, codim_, range_, rank_ >::
cbegin(const PropId &prop) const -> ElementIterator
{
  return ElementIterator(this->create_element_begin(prop));
}



template<int dim_, int codim_, int range_, int rank_>
auto
Function<dim_, codim_, range_, rank_ >::
cend(const PropId &prop) const -> ElementIterator
{
  return ElementIterator(this->create_element_end(prop));
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


template<int dim_, int codim_, int range_, int rank_>
Index
Function<dim_, codim_, range_, rank_ >::
get_object_id() const
{
  return object_id_;
}




#ifdef MESH_REFINEMENT
template<int dim_, int codim_, int range_, int rank_>
auto
Function<dim_, codim_, range_, rank_ >::
get_function_previous_refinement() const -> const std::shared_ptr<const self_t> &
{
  return function_previous_refinement_;
}


template<int dim_, int codim_, int range_, int rank_>
boost::signals2::connection
Function<dim_, codim_, range_, rank_ >::
connect_insert_knots(const typename Grid<dim_>::SignalInsertKnotsSlot &subscriber)
{
  return domain_.data_is_const() ?
         std::const_pointer_cast<DomainType>(domain_.get_ptr_const_data())
         ->connect_insert_knots(subscriber) :
         domain_.get_ptr_data()->connect_insert_knots(subscriber);
}

template<int dim_, int codim_, int range_, int rank_>
void
Function<dim_, codim_, range_, rank_ >::
create_connection_for_insert_knots(const std::shared_ptr<self_t> &function)
{
  Assert(function != nullptr, ExcNullPtr());
  Assert(&(*function) == &(*this), ExcMessage("Different objects."));

  auto func_to_connect =
    std::bind(&self_t::rebuild_after_insert_knots,
              function.get(),
              std::placeholders::_1,
              std::placeholders::_2);

  using SlotType = typename Grid<dim_>::SignalInsertKnotsSlot;
  this->connect_insert_knots(SlotType(func_to_connect).track_foreign(function));
}


#endif // MESH_REFINEMENT



#ifdef SERIALIZATION
template<int dim_, int codim_, int range_, int rank_>
template<class Archive>
void
Function<dim_, codim_, range_, rank_ >::
serialize(Archive &ar)
{
  ar &make_nvp("domain_",domain_);
  ar &make_nvp("name_",name_);
  ar &make_nvp("object_id_",object_id_);

#ifdef MESH_REFINEMENT
  auto tmp = std::const_pointer_cast<self_t>(function_previous_refinement_);
  ar &make_nvp("function_previous_refinement_",tmp);
  function_previous_refinement_ = tmp;
#endif // MESH_REFINEMENT
}
///@}
#endif // SERIALIZATION


IGA_NAMESPACE_CLOSE

#include <igatools/functions/function.inst>

