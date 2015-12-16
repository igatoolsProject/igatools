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

#include <igatools/geometry/domain.h>
#include <igatools/geometry/domain_element.h>
#include <igatools/utils/unique_id_generator.h>

using std::shared_ptr;
using std::to_string;

IGA_NAMESPACE_OPEN

template<int dim_, int codim_>
Domain<dim_, codim_>::
Domain(const SharedPtrConstnessHandler<GridFuncType> &func,
       const std::string &name)
  :
  grid_func_(func),
  name_(name),
  object_id_(UniqueIdGenerator::get_unique_id())
{}



template<int dim_, int codim_>
auto
Domain<dim_, codim_>::
create(const std::shared_ptr<GridFuncType> &func,
       const std::string &name)
-> std::shared_ptr<self_t>
{
  auto domain = std::shared_ptr<self_t>(
    new self_t(SharedPtrConstnessHandler<GridFuncType>(func),name));

#ifdef MESH_REFINEMENT
  domain->create_connection_for_insert_knots(domain);
#endif

  return domain;
}


template<int dim_, int codim_>
auto
Domain<dim_, codim_>::
const_create(const std::shared_ptr<const GridFuncType> &func,
             const std::string &name)
-> std::shared_ptr<const self_t>
{
  return std::shared_ptr<self_t>(
    new self_t(SharedPtrConstnessHandler<GridFuncType>(func),name));
}


#ifdef MESH_REFINEMENT
template<int dim_, int codim_>
void
Domain<dim_, codim_>::
create_connection_for_insert_knots(const std::shared_ptr<self_t> &domain)
{
  Assert(domain != nullptr, ExcNullPtr());
  Assert(&(*domain) == &(*this), ExcMessage("Different objects."));

  auto func_to_connect =
    std::bind(&self_t::rebuild_after_insert_knots,
              domain.get(),
              std::placeholders::_1,
              std::placeholders::_2);

  using SlotType = typename Grid<dim_>::SignalInsertKnotsSlot;
//  grid_func_.get_ptr_data()->get_grid()
//  ->connect_insert_knots(SlotType(func_to_connect).track_foreign(domain));
  connect_insert_knots(SlotType(func_to_connect).track_foreign(domain));
}


template<int dim_, int codim_>
void
Domain<dim_, codim_>::
rebuild_after_insert_knots(
  const SafeSTLArray<SafeSTLVector<Real>,dim_> &knots_to_insert,
  const Grid<dim_> &old_grid)
{
  this->domain_previous_refinement_ =
    Domain<dim_,codim_>::const_create(
      this->grid_func_->get_grid_function_previous_refinement());
}

template<int dim_, int codim_>
boost::signals2::connection
Domain<dim_, codim_>::
connect_insert_knots(const typename Grid<dim_>::SignalInsertKnotsSlot &subscriber)
{
//  grid_func_.get_ptr_data()->get_grid()->connect_insert_knots(subscriber);
  return grid_func_.get_ptr_data()->connect_insert_knots(subscriber);
}


#endif


template<int dim_, int codim_>
auto
Domain<dim_, codim_>::
get_grid_function() const -> std::shared_ptr<const GridFuncType>
{
  return grid_func_.get_ptr_const_data();
}


template<int dim_, int codim_>
const std::string &
Domain<dim_, codim_>::
get_name() const
{
  return name_;
}


template<int dim_, int codim_>
void
Domain<dim_, codim_>::
set_name(const std::string &name)
{
  name_ = name;
}

template<int dim_, int codim_>
int
Domain<dim_, codim_>::
get_object_id() const
{
  return object_id_;
}


template<int dim_, int codim_>
auto
Domain<dim_, codim_>::
create_cache_handler() const
-> std::unique_ptr<ElementHandler>
{
  return std::unique_ptr<ElementHandler>(new ElementHandler(this->shared_from_this()));
}


#if 0
template<int dim_, int codim_>
auto
Domain<dim_, codim_>::
create_element(const ListIt &index, const PropId &prop) const
-> std::unique_ptr<ElementAccessor>
{
  std::unique_ptr<ElementAccessor> elem;

  const auto &elements_with_property =
  this->get_grid_function()->get_grid()->get_elements_with_property(prop);

  if (&(*index) != &(*elements_with_property.end()))
  {
    elem = this->create_element_begin(prop);
    elem->move_to(*index);
  }
  else
  {
    elem = this->create_element_end(prop);
  }

  return std::move(elem);
}
#endif

template<int dim_, int codim_>
auto
Domain<dim_, codim_>::
create_element_begin(const PropId &prop) const
-> std::unique_ptr<ElementAccessor>
{
  using Elem = DomainElement<dim_,codim_>;
  return std::make_unique<Elem>(
    this->shared_from_this(),
    grid_func_->create_element_begin(prop));
}

template<int dim_, int codim_>
auto
Domain<dim_, codim_>::
create_element_end(const PropId &prop) const
-> std::unique_ptr<ElementAccessor>
{
  using Elem = DomainElement<dim_,codim_>;
  return std::make_unique<Elem>(
    this->shared_from_this(),
    grid_func_->create_element_end(prop));
}

template<int dim_, int codim_>
auto
Domain<dim_, codim_>::
begin(const PropId &prop) const -> ElementIterator
{
  return this->cbegin(prop);
}



template<int dim_, int codim_>
auto
Domain<dim_, codim_>::
end(const PropId &prop) const -> ElementIterator
{
  return this->cend(prop);
}



template<int dim_, int codim_>
auto
Domain<dim_, codim_>::
cbegin(const PropId &prop) const -> ElementIterator
{
  return ElementIterator(this->create_element_begin(prop));
}



template<int dim_, int codim_>
auto
Domain<dim_, codim_>::
cend(const PropId &prop) const -> ElementIterator
{
  return ElementIterator(this->create_element_end(prop));
}


template<int dim_, int codim_>
void
Domain<dim_, codim_>::
print_info(LogStream &out) const
{
  using std::to_string;
  out.begin_item("GridFunction<"
                 + to_string(dim)
                 + ","
                 + to_string(space_dim)
                 +">");
  grid_func_->print_info(out);
  out.end_item();

  if (this->name_.size() > 0)
      out << "Name: " << this->name_ << std::endl;

//    AssertThrow(false,ExcNotImplemented());
}


IGA_NAMESPACE_CLOSE

#include <igatools/geometry/domain.inst>

