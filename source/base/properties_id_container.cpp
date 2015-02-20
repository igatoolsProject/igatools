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



const std::string PropertiesIdContainer::property_none = "";

bool
PropertiesIdContainer::
test_id_for_property(const Index id, const std::string &property) const
{
    if (property == property_none)
    {
        return true; // an element can always be considered without any property
    }
    else
    {
        const auto &ids_same_property = this->get_ids_same_property(property);
        return std::binary_search(ids_same_property.begin(),ids_same_property.end(),id);
    }
}


void
PropertiesIdContainer::
add_property(const std::string &property)
{
    Assert(properties_id_.count(property) == 0,
           ExcMessage("The property \"" + property + "\" is already defined."));
    properties_id_[property] = std::set<int>();
}


std::set<Index> &
PropertiesIdContainer::
get_ids_same_property(const std::string &property)
{
    Assert(properties_id_.count(property) > 0,
           ExcMessage("The property \"" + property + "\" is not defined."));
    return properties_id_.at(property);
}



const std::set<Index> &
PropertiesIdContainer::
get_ids_same_property(const std::string &property) const
{
    Assert(properties_id_.count(property) > 0,
           ExcMessage("The property \"" + property + "\" is not defined."));
    return properties_id_.at(property);
}




void
PropertiesIdContainer::
set_id_property_status(const std::string &property,
                       const Index id,
                       const bool status)
{
    auto &ids_same_property = get_ids_same_property(property);
    if (status)
    {
        ids_same_property.insert(id);
    }
    else
    {
        Assert(!ids_same_property.empty(),ExcEmptyObject());
        ids_same_property.erase(id);
    }
}


void
PropertiesIdContainer::
set_ids_property_status(const std::string &property,
                        const std::set<Index> ids,
                        const bool status)
{
    for (const auto id : ids)
        set_id_property_status(property,id,status);
}


auto
PropertiesIdContainer::
begin() -> iterator
{
    return properties_id_.begin();
}


auto
PropertiesIdContainer::
end() -> iterator
{
    return properties_id_.end();
}

auto
PropertiesIdContainer::
begin() const -> const_iterator
{
    return properties_id_.begin();
}


auto
PropertiesIdContainer::
end() const -> const_iterator
{
    return properties_id_.end();
}

IGA_NAMESPACE_CLOSE
