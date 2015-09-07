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

#include <igatools/base/properties_id_container.h>
#include <igatools/base/exceptions.h>

#include <algorithm>

IGA_NAMESPACE_OPEN

template <typename IndexType>
bool
PropertiesIdContainer<IndexType>::
is_property_defined(const PropId &property) const
{
  return (properties_id_.count(property) > 0) ? true : false;
}



template <typename IndexType>
bool
PropertiesIdContainer<IndexType>::
test_id_for_property(const IndexType id, const PropId &property) const
{
  const auto &ids_same_property = (*this)[property];
  return std::binary_search(ids_same_property.begin(),ids_same_property.end(),id);
}



template <typename IndexType>
void
PropertiesIdContainer<IndexType>::
add_property(const PropId &property)
{
  Assert(!is_property_defined(property), ExcPropAlreadyDefined(property));
  properties_id_[property] = List();
}



template <typename IndexType>
auto
PropertiesIdContainer<IndexType>::
operator[](const PropId &property) -> List &
{
  Assert(is_property_defined(property), ExcPropNotDefined(property));
  return properties_id_.at(property);
}



template <typename IndexType>
auto
PropertiesIdContainer<IndexType>::
operator[](const PropId &property) const
-> const List &
{
  Assert(is_property_defined(property), ExcPropNotDefined(property));
  return properties_id_.at(property);
}



template <typename IndexType>
void
PropertiesIdContainer<IndexType>::
set_id_property_status(const PropId &property,
                       const IndexType id,
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



template <typename IndexType>
void
PropertiesIdContainer<IndexType>::
set_ids_property_status(const PropId &property,
                        const List ids,
                        const bool status)
{
  for (const auto id : ids)
    set_id_property_status(property,id,status);
}



template <typename IndexType>
auto
PropertiesIdContainer<IndexType>::
begin() -> iterator
{
  return properties_id_.begin();
}



template <typename IndexType>
auto
PropertiesIdContainer<IndexType>::
end() -> iterator
{
  return properties_id_.end();
}



template <typename IndexType>
auto
PropertiesIdContainer<IndexType>::
begin() const -> const_iterator
{
  return properties_id_.begin();
}



template <typename IndexType>
auto
PropertiesIdContainer<IndexType>::
end() const -> const_iterator
{
  return properties_id_.end();
}



template <typename IndexType>
SafeSTLVector<PropId>
PropertiesIdContainer<IndexType>::
get_properties() const
{
  SafeSTLVector<PropId> properties;

  for (const auto &ids_same_property : properties_id_)
    properties.emplace_back(ids_same_property.first);

  return properties;
}




//template <typename IndexType>
//void
//PropertiesIdContainer<IndexType>::
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



template <typename IndexType>
void
PropertiesIdContainer<IndexType>::
print_info(LogStream &out) const
{
  for (const auto &entry : properties_id_)
  {
    out.begin_item("IDs with property \"" + entry.first + "\":");
    entry.second.print_info(out);
    out.end_item();
  }
}



template <typename IndexType>
bool
PropertiesIdContainer<IndexType>::
empty() const noexcept
{
  return properties_id_.empty();
}

#if 0
#ifdef SERIALIZATION

template <typename IndexType>
template<class Archive>
void
PropertiesIdContainer<IndexType>::
serialize(Archive &ar, const unsigned int version)
{
  ar &make_nvp("properties_id_",properties_id_);
}

#endif // SERIALIZATION
#endif
IGA_NAMESPACE_CLOSE

#include <igatools/base/properties_id_container.inst>
