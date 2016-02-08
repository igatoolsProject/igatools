//-+--------------------------------------------------------------------
// Igatools a general purpose Isogeometric analysis library.
// Copyright (C) 2012-2016  by the igatools authors (see authors.txt).
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

#include <igatools/base/properties_id_container.h>
#include <igatools/base/exceptions.h>

#include <algorithm>

IGA_NAMESPACE_OPEN

template <typename IndexType,template <class T> class STLContainer>
bool
PropertiesIdContainer<IndexType,STLContainer>::
is_property_defined(const PropId &property) const
{
  return (properties_id_.count(property) > 0) ? true : false;
}



template <typename IndexType,template <class T> class STLContainer>
bool
PropertiesIdContainer<IndexType,STLContainer>::
test_id_for_property(const IndexType id, const PropId &property) const
{
  const auto &ids_same_property = (*this)[property];
  return std::binary_search(ids_same_property.begin(),ids_same_property.end(),id);
}



template <typename IndexType,template <class T> class STLContainer>
void
PropertiesIdContainer<IndexType,STLContainer>::
add_property(const PropId &property)
{
  Assert(!is_property_defined(property), ExcPropAlreadyDefined(property));
  properties_id_[property] = List();
}



template <typename IndexType,template <class T> class STLContainer>
auto
PropertiesIdContainer<IndexType,STLContainer>::
operator[](const PropId &property) -> List &
{
  Assert(is_property_defined(property), ExcPropNotDefined(property));
  return properties_id_.at(property);
}



template <typename IndexType,template <class T> class STLContainer>
auto
PropertiesIdContainer<IndexType,STLContainer>::
operator[](const PropId &property) const
-> const List &
{
  Assert(is_property_defined(property), ExcPropNotDefined(property));
  return properties_id_.at(property);
}






template <typename IndexType,template <class T> class STLContainer>
auto
PropertiesIdContainer<IndexType,STLContainer>::
begin() -> iterator
{
  return properties_id_.begin();
}



template <typename IndexType,template <class T> class STLContainer>
auto
PropertiesIdContainer<IndexType,STLContainer>::
end() -> iterator
{
  return properties_id_.end();
}



template <typename IndexType,template <class T> class STLContainer>
auto
PropertiesIdContainer<IndexType,STLContainer>::
begin() const -> const_iterator
{
  return properties_id_.begin();
}



template <typename IndexType,template <class T> class STLContainer>
auto
PropertiesIdContainer<IndexType,STLContainer>::
end() const -> const_iterator
{
  return properties_id_.end();
}



template <typename IndexType,template <class T> class STLContainer>
SafeSTLVector<PropId>
PropertiesIdContainer<IndexType,STLContainer>::
get_properties() const
{
  SafeSTLVector<PropId> properties;

  for (const auto &ids_same_property : properties_id_)
    properties.emplace_back(ids_same_property.first);

  return properties;
}




//template <typename IndexType,template <class T> class STLContainer>
//void
//PropertiesIdContainer<IndexType,STLContainer>::
//add_offset(const IndexType offset)
//{
//    for (auto &property_id : properties_id_)
//    {
//        const std::set<IndexType> &old_dofs = property_id.second;
//        std::set<IndexType> new_dofs;
//        for (const auto &dof : old_dofs)
//            new_dofs.insert(dof + offset);
//
//        property_id.second = std::move(new_dofs);
////      property_id.second = new_dofs;
//    }
//}



template <typename IndexType,template <class T> class STLContainer>
void
PropertiesIdContainer<IndexType,STLContainer>::
print_info(LogStream &out) const
{
  for (const auto &entry : properties_id_)
  {
    out.begin_item("IDs with property \"" + entry.first + "\":");
    entry.second.print_info(out);
    out.end_item();
  }
}



template <typename IndexType,template <class T> class STLContainer>
bool
PropertiesIdContainer<IndexType,STLContainer>::
empty() const noexcept
{
  return properties_id_.empty();
}




void
PropertiesDofs::
set_property_status_for_id(const PropId &property,
                           const int id,
                           const bool status)
{
  auto &list = (*this)[property];
  if (status)
  {
    list.insert(id);
  }
  else
  {
    Assert(!list.empty(),ExcEmptyObject());
    list.erase(id);
  }
}



void
PropertiesDofs::
set_property_status_for_ids(const PropId &property,
                            const List &ids,
                            const bool status)
{
  auto &list = (*this)[property];
  if (status)
  {
    list.insert(ids.begin(),ids.end());
  }
  else
  {
    Assert(!list.empty(),ExcEmptyObject());
    for (int id : ids)
      list.erase(id);
  }
}



template <int dim>
void
PropertiesElementID<dim>::
set_property_status_for_id(const PropId &property,
                           const ElementIndex<dim> &elem_id,
                           const bool status)
{
  auto &list = (*this)[property];
  if (status)
  {
    list.emplace_back(elem_id);
    std::sort(list.begin(),list.end());
  }
  else
  {
    Assert(!list.empty(),ExcEmptyObject());
    list.erase(std::find(list.cbegin(),list.cend(),elem_id));
  }
}



template <int dim>
void
PropertiesElementID<dim>::
set_property_status_for_ids(const PropId &property,
                            const List &ids,
                            const bool status)
{
  auto &list = (*this)[property];
  if (status)
  {
    for (const auto &elem_id : ids)
      list.emplace_back(elem_id);

    std::sort(list.begin(),list.end());
  }
  else
  {
    Assert(!list.empty(),ExcEmptyObject());
    for (const auto &elem_id : ids)
      list.erase(std::find(list.cbegin(),list.cend(),elem_id));
  }
}

IGA_NAMESPACE_CLOSE

#include <igatools/base/properties_id_container.inst>
