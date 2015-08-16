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
#include <igatools/base/logstream.h>
#include <igatools/base/properties.h>
#include <igatools/utils/safe_stl_vector.h>


#include <set>
#include <map>
#include <string>

IGA_NAMESPACE_OPEN
//TODO (pauletti, Aug 15, 2015): document this class
/**
 * The property name is the key of the std::map.
 *
 * @ingroup serializable
 */
class PropertiesIdContainer
{
public:
    /**
     * Returns TRUE if the @p property is defined.
     *
     * @warning Returns true also if the @p property is defined but no ids are associated to the @p property.
     */
    bool is_property_defined(const std::string &property) const;

    /**
     * Returns TRUE if the @p id has the given @p property.
     */
    bool test_id_for_property(const Index id, const std::string &property) const;


    /**
     * Adds a new <tt>property</tt> definition.
     *
     * @note If the <tt>property</tt> is already present,
     * an assertion will be raised (in Debug mode).
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

    /**
     * Sets the <tt>status</tt> of the given <tt>property</tt> for the given <tt>ids</tt>.
     */
    void set_ids_property_status(const std::string &property,
                                 const std::set<Index> ids,
                                 const bool status);


    /**
     * Prints the contents of the class. Its use is intended for testing and debugging purposes.
     */
    void print_info(LogStream &out) const;


    /**
     * Returns whether the container is empty (i.e. whether no properties are defined).
     * @note This function does not modify the container in any way.
     */
    bool empty() const noexcept;


    /**
     * Returns the properties defined.
     */
    SafeSTLVector<std::string> get_properties() const;


//    /**
//     * Add the @p offset value to the ids contained in the object.
//     * @param offset
//     */
//    void add_offset(const Index offset);

private:
    using ContainerType = std::map<std::string,std::set<Index>>;
    using iterator = typename ContainerType::iterator;
    using const_iterator = typename ContainerType::const_iterator;

public:
    iterator begin();
    iterator end();

    const_iterator begin() const;
    const_iterator end() const;

private:
    /** Property table */
    ContainerType properties_id_;


#ifdef SERIALIZATION
    /**
     * @name Functions needed for boost::serialization
     * @see <a href="http://www.boost.org/doc/libs/release/libs/serialization/">boost::serialization</a>
     */
    ///@{
    friend class boost::serialization::access;

    template<class Archive>
    void
    serialize(Archive &ar, const unsigned int version);
    ///@}
#endif //SERIALIZATION

};



IGA_NAMESPACE_CLOSE

#endif // #ifndef PROPERTIES_ID_CONTAINER_H_


