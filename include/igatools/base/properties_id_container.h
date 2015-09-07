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

#ifndef __PROPERTIES_ID_CONTAINER_H_
#define __PROPERTIES_ID_CONTAINER_H_

#include <igatools/base/config.h>
#include <igatools/base/properties.h>
#include <igatools/base/logstream.h>
#include <igatools/utils/safe_stl_set.h>
#include <igatools/utils/safe_stl_vector.h>

#include <map>

IGA_NAMESPACE_OPEN
//TODO (pauletti, Aug 15, 2015): document this class
/**
 * Container indexed by a property, storing in each entry
 * a list of of IndexType.
 *
 * @author martinelli 2014,2015
 * @author pauletti 2015
 *
 * @ingroup serializable
 */
template <typename IndexType>
class PropertiesIdContainer
{
public:
  using List = SafeSTLSet<IndexType>;
private:
  using ContainerType = std::map<PropId, List>;
  using iterator = typename ContainerType::iterator;
  using const_iterator = typename ContainerType::const_iterator;

public:
  /**
   * Adds a new <tt>property</tt> definition.
   *
   * @note If the <tt>property</tt> is already present,
   * an assertion will be raised (in Debug mode).
   */
  void add_property(const PropId &property);

  /**
   * Returns TRUE if the @p property is defined.
   *
   * @warning Returns true also if the @p property is defined but no
   * ids are associated to the @p property.
   */
  bool is_property_defined(const PropId &property) const;

  /**
   * Returns TRUE if the @p id has the given @p property.
   */
  bool test_id_for_property(const IndexType id, const PropId &property) const;

  /**
   * Returns the the set of IDs having a certain @p property (non-const version).
   */
  List &operator[](const PropId &property);

  /**
   * Returns the flat id of IDs having a certain @p property (const version).
   */
  const List &operator[](const PropId &property) const;

private:
  /**
   * Sets the <tt>status</tt> of the given <tt>property</tt> for the given <tt>id</tt>.
   */
  void set_id_property_status(const PropId &property,
                              const IndexType id,
                              const bool status);

  /**
   * Sets the <tt>status</tt> of the given <tt>property</tt> for the given <tt>ids</tt>.
   */
  void set_ids_property_status(const PropId &property,
                               const List ids,
                               const bool status);
public:

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
  SafeSTLVector<PropId> get_properties() const;


//    /**
//     * Add the @p offset value to the ids contained in the object.
//     * @param offset
//     */
//    void add_offset(const IndexType offset);


public:
  iterator begin();
  iterator end();

  const_iterator begin() const;
  const_iterator end() const;

private:
  /** Property table */
  ContainerType properties_id_;

  DeclException1(ExcPropNotDefined, PropId,
                 << "Property \"" << arg1 << "\" is not present.");
  DeclException1(ExcPropAlreadyDefined, PropId,
                 << "Property \"" << arg1 << "\" is already present.");

#ifdef SERIALIZATION
  /**
   * @name Functions needed for boost::serialization
   * @see <a href="http://www.boost.org/doc/libs/release/libs/serialization/">boost::serialization</a>
   */
  ///@{
  friend class boost::serialization::access;
  friend class serialization_access;

  template<class Archive>
  void
  serialize(Archive &ar, const unsigned int version)
  {
    ar &make_nvp("properties_id_",properties_id_);
  }

  ///@}
#endif //SERIALIZATION

};



IGA_NAMESPACE_CLOSE

#endif // #ifndef PROPERTIES_ID_CONTAINER_H_


