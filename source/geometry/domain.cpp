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

using std::shared_ptr;

IGA_NAMESPACE_OPEN

template<int dim_, int codim_>
Domain<dim_, codim_>::
Domain(std::shared_ptr<GridType> grid)
  :
  grid_(grid)
{
  Assert(grid_ != nullptr, ExcNullPtr());
}



template<int dim_, int codim_>
Domain<dim_, codim_>::
~Domain()
{}



template<int dim_, int codim_>
auto
Domain<dim_, codim_>::
get_grid() const -> std::shared_ptr<GridType>
{
  return grid_;
}



template<int dim_, int codim_>
auto
Domain<dim_, codim_>::
create_cache_handler() const
-> std::shared_ptr<ElementHandler>
{
  return std::make_shared<ElementHandler>(this->shared_from_this());
}



template<int dim_, int codim_>
auto
Domain<dim_, codim_>::
create_element(const ListIt &index, const PropId &prop) const
-> std::unique_ptr<ConstElementAccessor>
{
  using Elem = ConstElementAccessor;
  auto elem = std::make_unique<Elem>(this->shared_from_this(), index, prop);
  Assert(elem != nullptr, ExcNullPtr());

  return elem;
}



template<int dim_, int codim_>
auto
Domain<dim_, codim_>::
create_element(const ListIt &index, const PropId &prop)
-> std::unique_ptr<ElementAccessor>
{
  using Elem = ElementAccessor;
  auto elem = std::make_unique<Elem>(this->shared_from_this(), index, prop);
  Assert(elem != nullptr, ExcNullPtr());

  return elem;
}



template<int dim_, int codim_>
auto
Domain<dim_, codim_>::
begin(const PropId &prop) -> ElementIterator
{
  return ElementIterator(this->shared_from_this(),
  grid_->get_elements_with_property(prop).begin(), prop);
}



template<int dim_, int codim_>
auto
Domain<dim_, codim_>::
end(const PropId &prop) -> ElementIterator
{
  return ElementIterator(this->shared_from_this(),
  grid_->get_elements_with_property(prop).end(),
  prop);
}



template<int dim_, int codim_>
auto
Domain<dim_, codim_>::
begin(const PropId &prop) const -> ElementConstIterator
{
  return this->cbegin(prop);
}



template<int dim_, int codim_>
auto
Domain<dim_, codim_>::
end(const PropId &prop) const -> ElementConstIterator
{
  return this->cend(prop);
}



template<int dim_, int codim_>
auto
Domain<dim_, codim_>::
cbegin(const PropId &prop) const -> ElementConstIterator
{
  return ElementConstIterator(this->shared_from_this(),
                              grid_->get_elements_with_property(prop).begin(),
                              prop);
}



template<int dim_, int codim_>
auto
Domain<dim_, codim_>::
cend(const PropId &prop) const -> ElementConstIterator
{
  return ElementConstIterator(this->shared_from_this(),
                              grid_->get_elements_with_property(prop).end(),
                              prop);
}

IGA_NAMESPACE_CLOSE

#include <igatools/geometry/domain.inst>

