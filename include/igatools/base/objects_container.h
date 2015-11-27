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

#ifndef __OBJECTS_CONTAINER_H_
#define __OBJECTS_CONTAINER_H_

#include <igatools/base/config.h>

#include <igatools/base/instantiated_types.inst>
#include <igatools/utils/safe_stl_vector.h>

#include <boost/fusion/container/vector.hpp>
#include <boost/fusion/container/map.hpp>

#include <boost/fusion/algorithm/transformation/transform.hpp>
#include <boost/fusion/include/transform.hpp>
#include <boost/mpl/copy.hpp>

IGA_NAMESPACE_OPEN

class LogStream;

/**
 * @brief Container class for igatools objects.
 *
 * This class is a container for high-level igatools objects.
 *
 * The objects are stored inside wrapped into shared pointers, that is,
 * the container is actually storing shared pointers.
 * And they can be stored as shared pointer to constant o non constant
 * objects.
 *
 * The objects that can be stored inside are
 * - @ref Grid
 * - @ref SplineSpace
 * - @ref ReferenceSpaceBasis
 * - @ref GridFunction
 * - @ref Domain
 * - @ref PhysicalSpaceBasis
 * - @ref Function

 * All the instantiated dimensions of the above listed types can be
 * stored. In @ref InstantiatedTypes a list of the valid types is defined.
 *
 * To every single type (const or non constant) a @p vector of shared
 * pointer to (constant o non constant) object of the type is associated.
 *
 * The class provides methods for inserting and getting objects,
 * and also for querying if an object of a given type and given object
 * unique @p Id number is present.
 *
 * This container can be create void and filled by using the
 * @p insert_object method, or it can be parsed directly from an XML file
 * by using the class @ref ObjectsContainerParser.
 *
 * This class is not intended to provide multi-patch support, dofs
 * management, or similar features. It is just a container of pointers
 * for different types.
 *
 * @see InstantiatedTypes
 *
 * @author P. Antolin
 * @date 2015
 */
class ObjectsContainer
{
private:

  /** Type for current class. */
  using self_t = ObjectsContainer;

  template< class T >
  struct as_fusion_vector_shared_ptr
  {
      /**
       * This functor transform a <tt>boost::mpl::vector</tt> of type into a
       *  <tt>boost::fusion::vector</tt> of <tt>shared_ptr</tt> of the types.
       */

      typedef typename boost::fusion::result_of::as_vector<
        typename boost::mpl::transform<T, std::shared_ptr<boost::mpl::_1>>::type>::type type;
  };

  /** Alias for all instantiated grids. */
  using Grids        = typename InstantiatedTypes::Grids;

  /** Alias for all instantiated spline spaces. */
  using SpSpaces = typename InstantiatedTypes::SplineSpaces;

  /** Alias for all instantiated spline spaces. */
  using RefSpaces = typename InstantiatedTypes::RefSpaceBases;

  /** Alias for all instantiated grids functions. */
  using GridFunc     = typename InstantiatedTypes::GridFunctions;

  /** Alias for all instantiated domains. */
  using Domains      = typename InstantiatedTypes::Domains;

  /** Alias for all instantiated physical space basis. */
  using PhysSpaces   = typename InstantiatedTypes::PhysSpaces;

  /** Alias for all instantiated  functions.*/
  using Functions    = typename InstantiatedTypes::Functions;

public:

  /** <tt>fusion::vector</tt> of <tt>shared_ptr</tt> of all instantiated grids. */
  using GridPtrs      = as_fusion_vector_shared_ptr<Grids>::type;

  /** <tt>fusion::vector</tt> of <tt>shared_ptr</tt> of all instantiated spline spaces. */
  using SpSpacePtrs   = as_fusion_vector_shared_ptr<SpSpaces>::type;

  /** <tt>fusion::vector</tt> of <tt>shared_ptr</tt> of all instantiated reference spaces. */
  using RefSpacePtrs   = as_fusion_vector_shared_ptr<RefSpaces>::type;

  /** <tt>fusion::vector</tt> of <tt>shared_ptr</tt> of all instantiated grid functions. */
  using GridFuncPtrs  = as_fusion_vector_shared_ptr<GridFunc>::type;

  /** <tt>fusion::vector</tt> of <tt>shared_ptr</tt> of all instantiated domains. */
  using DomainPtrs    = as_fusion_vector_shared_ptr<Domains>::type;

  /** <tt>fusion::vector</tt> of <tt>shared_ptr</tt> of all instantiated physical spaces. */
  using PhysSpacePtrs = as_fusion_vector_shared_ptr<PhysSpaces>::type;

  /** <tt>fusion::vector</tt> of <tt>shared_ptr</tt> of all instantiated functions. */
  using FunctionPtrs  = as_fusion_vector_shared_ptr<Functions>::type;


private:
  /** <tt>mpl::vector</tt> containing all the instantiated types together. */
  typedef boost::mpl::copy<Grids,
          boost::mpl::back_inserter<boost::mpl::copy<SpSpaces,
          boost::mpl::back_inserter<boost::mpl::copy<RefSpaces,
          boost::mpl::back_inserter<boost::mpl::copy<GridFunc,
          boost::mpl::back_inserter<boost::mpl::copy<Domains,
          boost::mpl::back_inserter<boost::mpl::copy<PhysSpaces,
          boost::mpl::back_inserter<Functions>
      >::type> >::type> >::type> >::type> >::type> >::type JointTypes_;

  /** @p JointTypes_ converted to constant types. */
  typedef boost::mpl::transform<JointTypes_,
          boost::add_const<boost::mpl::_1> >::type JointConstTypes_;

  /** @p All types, const and not const, together in single <tt>mpl::vector</tt>. */
  typedef boost::mpl::copy<JointTypes_,
          boost::mpl::back_inserter<JointConstTypes_>>::type AllTypes_;

  template <class T>
  struct as_fusion_map
  {
      /**
       * This functor transform a sequence of types @p T into a
       * <tt>boost::fusion::map</tt> composed of <tt>pair</tt>s of the form
       * <tt>pair<T, SafeSTLVector<shared_ptr<T>></tt>.
       */

  private:
      template <class S>
      using Pair_ = boost::fusion::pair<S, SafeSTLVector<std::shared_ptr<S>>>;

  public:
      typedef typename boost::fusion::result_of::as_map<
        typename boost::mpl::transform<T, Pair_<boost::mpl::_1>>::type>::type type;
  };

  /** Container for shared pointer of all the instantiated types. */
  using ObjectMapTypes_ = as_fusion_map<AllTypes_>::type;


private:

  /** @name Constructors*/
  ///@{
  /**
   * @brief Default constructor.
   *
   * Default constructor.
   */
  ObjectsContainer() = default;

  /**
   * Copy constructor.
   */
  ObjectsContainer(const self_t &container) = default;

  /**  Move constructor */
  ObjectsContainer(self_t &&container) = default;

public:
  /** Destructor */
  ~ObjectsContainer() = default;
  ///@}

public:
  /**
   * @name Creators
   */
  ///@{
  /**
   * Creates a objects container (non-const).
   */
  static std::shared_ptr<self_t> create();
  ///@}

private:
  /**
   * @name Assignment operators
   */
  ///@{
  /**
   * Copy assignment operator.
   */
  self_t &operator=(const self_t &container) = default;

  /**
   * Move assignment operator.
   */
  self_t &operator=(self_t &&container) = default;
  ///@}

public:

  /**
   * @brief Insert a pointer to an object into the container.
   *
   * Insert a @ref shared_ptr pointing to an object of type @ref T
   * into the container.
   *
   * In debug mode, before inserting it, it is checked that this object
   * is not already present.
   *
   * @tparam T Type of the object to be inserted.
   * @param[in] object Object to be inserted.
   */
  template <class T>
  void insert_object (const std::shared_ptr<T> object);

  /**
   * @brief Insert a pointer to a constant object into the container.
   *
   * Insert a @ref shared_ptr pointing to a constant object of type @ref T
   * into the container.
   *
   * In debug mode, before inserting it, it is checked that this object
   * is not already present.
   *
   * @tparam T Type of the object to be inserted.
   * @param[in] object Object to be inserted.
   */
  template <class T>
  void insert_const_object (const std::shared_ptr<const T> object);

  /**
   * @brief Retrieves a pointer to an object from the container.
   *
   * Retrieves a pointer to an object of type @ref T that has the
   * given unique @ref id number.
   *
   * In debug mode, it there no exists an object of type @ref T
   * with the given unique @ref id, an exception is thrown.
   *
   * @tparam T Type of the object to be retrieved.
   * @param[in] id Unique Id of the queried object.
   * @return Object retrieved.
   */
  template <class T> std::shared_ptr<T> get_object (const Index &id) const;

  /**
   * @brief Retrieves a pointer to a constant object from the container.
   *
   * Retrieves a pointer to a constant object of type @ref T that has the
   * given unique @ref id number.
   *
   * In debug mode, it there no exists a const object of type @ref T
   * with the given unique @ref id, an exception is thrown.
   *
   * @tparam T Type of the object to be retrieved.
   * @param[in] id Unique Id of the queried object.
   * @return Object retrieved.
   */
  template <class T> std::shared_ptr<const T> get_const_object (const Index &id) const;

  /**
   * @brief Checks if the a object of the type @ref T
   * with the given unique @ref id is contained inside.
   *
   * Checks if the a object of the type @ref T
   * with the given unique @ref id is contained inside.
   *
   * @tparam T Type of the objects to be checked.
   * @param[in] id Unique id to checked.
   * @return @p true if there is an object with the given unique @ref id,
   *         @p false elsewhere.
   */
  template <class T> bool is_object_present (const Index &id) const;

  /**
   * @brief Checks if the a constant object of the type @ref T
   * with the given unique @ref id is contained inside.
   *
   * Checks if the a constant object of the type @ref T
   * with the given unique @ref id is contained inside.
   *
   * @tparam T Type of the objects to be checked.
   * @param[in] id Unique id to checked.
   * @return @p true if there is an object with the given unique @ref id,
   *         @p false elsewhere.
   */
  template <class T> bool is_const_object_present (const Index &id) const;

  /**
   * @brief Prints the objects contained inside.
   *
   * For every object type prints the objects contained (if any) by
   * calling its @p print_info method.
   * Mostly used for testing and debugging purposes.
   *
   * @param[in] out Log stream for outputting the information.
   */
  void print_info (LogStream &out) const;

private:

  /** Container for the objects. */
  ObjectMapTypes_ objects_;

};

IGA_NAMESPACE_CLOSE

#endif /*__OBJECTS_CONTAINER_H_ */
