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

#ifndef PROPERTIES_ID_CONTAINER_H_
#define PROPERTIES_ID_CONTAINER_H_

#include <igatools/base/config.h>

#include <set>
#include <map>
#include <string>


IGA_NAMESPACE_OPEN

class PropertiesIdContainer
{
public:
	static const std::string property_none;


    bool test_id_for_property(const Index id, const std::string &property) const;


    /**
     * Adds a new <tt>property</tt> definition.
     *
     * @note If the <tt>property</tt> is already present, an assertion will be raised (in Debug mode).
     */
    void add_property(const std::string &property);

    /**
     * Returns the the set of IDs having a certain @p property (non-const version).
     */
    std::set<Index> &get_ids_same_property(const std::string &property);

    /**
     * Returns the flat id of IDs having a certain @p property (const version).
     */
    const std::set<Index> &get_ids_same_property(const std::string &property) const;

    /**
     * Sets the <tt>status</tt> of the given <tt>property</tt> for the given <tt>id</tt>.
     */
    void set_id_property_status(const std::string &property,
                                     const Index id,
                                     const bool status);

private:
    using ContainerType = std::map<std::string,std::set<Index>>;
    using iterator = typename ContainerType::iterator;

public:
    iterator begin();
    iterator end();

private:
	std::map<std::string,std::set<Index>> properties_id_;
};



IGA_NAMESPACE_CLOSE

#endif // #ifndef PROPERTIES_ID_CONTAINER_H_


