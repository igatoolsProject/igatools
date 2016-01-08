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

#ifdef XML_IO

#include <igatools/base/instantiated_types.inst>
#include <igatools/utils/safe_stl_vector.h>

#include <boost/fusion/container/vector.hpp>
#include <boost/fusion/container/map.hpp>

#include <boost/fusion/algorithm/transformation/transform.hpp>
#include <boost/fusion/include/transform.hpp>
#include <boost/mpl/copy.hpp>

IGA_NAMESPACE_OPEN

class LogStream;
template <class T> class SafeSTLSet;

/**
 * @brief Container class for high-level igatools objects.
 *
 * The stored objects are wrapped into shared pointers, that is,
 * the container is actually storing shared pointers.
 * And they can be stored as shared pointer to constant o non constant
 * objects.
 *
 * The objects that can be stored inside are:
 * - @ref Grid
 * - @ref SplineSpace
 * - @ref ReferenceSpaceBasis
 * - @ref GridFunction
 * - @ref Domain
 * - @ref PhysicalSpaceBasis
 * - @ref Function

 * All the instantiated dimensions of the above listed types can be
 * stored. In @p InstantiatedTypes a list of the valid types is defined.
 *
 * To every single type (constant or non constant) a @ref SafeSTLVector of
 * shared pointers to objects of the type is associated.
 *
 * The class provides methods for inserting and getting objects,
 * and also for querying if an object of a given type and unique
 * object @p Id number is present.
 *
 * This container can be created void and filled after by using the
 * @ref insert_object method, or it can be parsed directly from an XML
 * file by using the class @ref ObjectsContainerXMLReader.
 *
 * This class is not intended to provide multi-patch support, dofs
 * management, or similar features. It is just a container of pointers
 * for different types.
 *
 * @ingroup serializable
 * @see ObjectsContainerXMLReader
 * @see ObjectsContainerXMLWriter
 *
 * @author P. Antolin
 * @date 2015
 */
class ObjectsContainer
{
private:

  /**
   * Template-alias used to simplify the type definitions.
   */
  template <class A,class B>
  using boost_mpl_transform_t = typename boost::mpl::transform<A,B>::type;

  /**
   * Alias used to simplify the type definitions.
   */
  using boost_mpl_1 = typename boost::mpl::_1;


  /** Type for current class. */
  using self_t = ObjectsContainer;


  /**
   * This template-alias behaves like a functor and transforms a
   * <tt>boost::mpl::vector</tt> of types into a
   * <tt>boost::fusion::vector</tt> of <tt>shared_ptr</tt>s of the types <tt>T</tt>.
   */
  template< class T >
  using as_fusion_vector_shared_ptr_t =
    typename boost::fusion::result_of::as_vector<
    boost_mpl_transform_t<T,std::shared_ptr<boost_mpl_1>> >::type;


  /** Alias for all instantiated grids. */
  using Grids = typename InstantiatedTypes::Grids;

  /** Alias for all instantiated spline spaces. */
  using SpSpaces = typename InstantiatedTypes::SplineSpaces;

  /** Alias for all instantiated reference space bases. */
  using RefSpaces = typename InstantiatedTypes::RefSpaceBases;

  /** Alias for all instantiated grids functions. */
  using GridFunc = typename InstantiatedTypes::GridFunctions;

  /** Alias for all instantiated domains. */
  using Domains = typename InstantiatedTypes::Domains;

  /** Alias for all instantiated physical space bases. */
  using PhysSpaces = typename InstantiatedTypes::PhysSpaces;

  /** Alias for all instantiated  functions.*/
  using Functions = typename InstantiatedTypes::Functions;

public:
  /**
   * <tt>boost::fusion::vector</tt> of <tt>shared_ptr</tt> of all
   * instantiated grids.
   */
  using GridPtrs = as_fusion_vector_shared_ptr_t<Grids>;

  /**
   * <tt>boost::fusion::vector</tt> of <tt>shared_ptr</tt> of all
   * instantiated spline spaces.
   */
  using SpSpacePtrs = as_fusion_vector_shared_ptr_t<SpSpaces>;

  /**
   * <tt>boost::fusion::vector</tt> of <tt>shared_ptr</tt> of all
   * instantiated reference space bases.
   */
  using RefSpacePtrs = as_fusion_vector_shared_ptr_t<RefSpaces>;

  /**
   * <tt>boost::fusion::vector</tt> of <tt>shared_ptr</tt> of all
   * instantiated grid functions.
   */
  using GridFuncPtrs = as_fusion_vector_shared_ptr_t<GridFunc>;

  /**
   * <tt>boost::fusion::vector</tt> of <tt>shared_ptr</tt> of all
   * instantiated domains.
   */
  using DomainPtrs = as_fusion_vector_shared_ptr_t<Domains>;

  /**
   * <tt>boost::fusion::vector</tt> of <tt>shared_ptr</tt> of all
   * instantiated physical spaces.
   */
  using PhysSpacePtrs = as_fusion_vector_shared_ptr_t<PhysSpaces>;

  /**
   * <tt>boost::fusion::vector</tt> of <tt>shared_ptr</tt> of all
   * instantiated functions.
   */
  using FunctionPtrs = as_fusion_vector_shared_ptr_t<Functions>;


private:


  /**
   * Template-alias used to simplify the type definitions.
   */
  template <class A,class B>
  using boost_mpl_copy_t = typename boost::mpl::copy<A,boost::mpl::back_inserter<B>>::type;


  /**
   * <tt>boost::mpl::vector</tt> containing all the instantiated types
   * together.
   */
  typedef boost_mpl_copy_t<Grids,
          boost_mpl_copy_t<SpSpaces,
          boost_mpl_copy_t<RefSpaces,
          boost_mpl_copy_t<GridFunc,
          boost_mpl_copy_t<Domains,
          boost_mpl_copy_t<PhysSpaces,Functions>
          > > > > > JointTypes_;

  /** @ref JointTypes_ converted to constant types. */
  using JointConstTypes_ = boost_mpl_transform_t<JointTypes_,boost::add_const<boost_mpl_1> >;

  /**
   * @ref JointTypes_ and @ref JointConstTypes_ together in single
   * <tt>boost::mpl::vector</tt>.
   */
  using AllTypes_ = boost_mpl_copy_t<JointTypes_,JointConstTypes_>;


  template <class S>
  using Pair_ = boost::fusion::pair<S,SafeSTLVector<std::shared_ptr<S>>>;

  /**
   * This template-alias beahves like a functor and transforms
   * a sequence of types @ref T into a
   * <tt>boost::fusion::map</tt> composed of <tt>boost::fusion::pair</tt>s of the form
   * <tt>boost::fusion::pair<T, SafeSTLVector<shared_ptr<T>></tt>.
   */
  template <class T>
  using as_fusion_map_t =
    typename boost::fusion::result_of::as_map<
    boost_mpl_transform_t<T,Pair_<boost_mpl_1> > >::type;

  /** Container for shared pointers of all the instantiated types. */
  using ObjectMapTypes_ = as_fusion_map_t<AllTypes_>;


public:

  /** @name Constructors*/
  ///@{
  /** Default constructor. */
  ObjectsContainer() = default;

  /** Copy constructor. */
  ObjectsContainer(const self_t &container) = default;

  /**  Move constructor. */
  ObjectsContainer(self_t &&container) = default;

  /** Destructor. */
  ~ObjectsContainer() = default;
  ///@}

  /**
   * @name Creators
   */
  ///@{
  /** Creates an objects container wrapped into a shared pointer. */
  static std::shared_ptr<self_t> create();
  ///@}

private:
  /**
   * @name Assignment operators
   */
  ///@{
  /** Copy assignment operator. */
  self_t &operator=(const self_t &container) = default;

  /** Move assignment operator. */
  self_t &operator=(self_t &&container) = default;
  ///@}

public:

  /**
   * @brief Insert a @ref shared_ptr pointing to an object of type @p T
   * into the container.
   *
   * @tparam T Type of the object to be inserted.
   * @param[in] object Shared pointer of the object to be inserted.
   * @param[in] check_present If true, before inserting the object it is
   *            checked if it has been already inserted, throwing an
   *            assert it that case.
   */
  template <class T>
  void insert_object(const std::shared_ptr<T> object,
                     const bool check_present = false);

  /**
   * @brief Insert a @ref shared_ptr pointing to a constant object of
   * type @p T into the container.
   *
   * @tparam T Type of the object to be inserted.
   * @param[in] object Object to be inserted.
   * @param[in] check_present If true, before inserting the object it is
   *            checked if it has been already inserted, throwing an
   *            assert it that case.
   */
  template <class T>
  void insert_const_object(const std::shared_ptr<const T> object,
                           const bool check_present = false);

  /**
   * @brief Retrieves a pointer to an object of type @p T that has the
   * given unique @p id number.
   *
   * @warning In debug mode, it there no exists an object of type @p T
   * with the given unique @p id, an exception is thrown.
   *
   * @tparam T Type of the object to be retrieved.
   * @param[in] id Unique Id of the queried object.
   * @return Object retrieved.
   */
  template <class T> std::shared_ptr<T> get_object(const Index &id) const;

  /**
   * @brief Retrieves a pointer to a constant object of type @p T that
   * has the given unique @p id number.
   *
   * @warning In debug mode, it there no exists a const object of type
   * @p T with the given unique @p id, an exception is thrown.
   *
   * @tparam T Type of the object to be retrieved.
   * @param[in] id Unique Id of the queried object.
   * @return Object retrieved.
   */
  template <class T> std::shared_ptr<const T> get_const_object(const Index &id) const;

  /**
   * @brief Retrieves a @ref SafeSTLSet with the ids of all the non
   * constant objects contained of type @p T.
   *
   * If there is no object of the given type, it returns a void set.
   *
   * @tparam T Type of the object to be retrieved.
   * @return Set of object ids of the given type.
   */
  template <class T> SafeSTLSet<Index> get_object_ids() const;

  /**
   * @brief Retrieves a @ref SafeSTLSet with the ids of all the constant
   * objects contained of type @p T.
   *
   * If there is no object of the given type, it returns a void set.
   *
   * @tparam T Type of the object to be retrieved.
   * @return Set of object ids of the given type.
   */
  template <class T> SafeSTLSet<Index> get_const_object_ids() const;

  /**
   * @brief Checks if the a object of the type @p T
   * with the given unique @p id is contained inside.
   *
   * @tparam T Type of the objects to be checked.
   * @param[in] id Unique id to checked.
   * @return @p true if there is an object with the given unique @p id,
   *         @p false elsewhere.
   */
  template <class T> bool is_object_present(const Index &id) const;

  /**
   * @brief Checks if the a constant object of the type @p T
   * with the given unique @p id is contained inside.
   *
   * @tparam T Type of the objects to be checked.
   * @param[in] id Unique id to checked.
   * @return @p true if there is an object with the given unique @p id,
   *         @p false elsewhere.
   */
  template <class T> bool is_const_object_present(const Index &id) const;

  /**
   * @brief Check all the dependencies of the already inserted objects
   * and insert them too.
   */
  void fill_not_inserted_dependencies();

  /**
   * @brief Prints the objects contained inside.
   *
   * For every object type prints the objects contained (if any) by
   * calling its @p print_info method.
   * Mostly used for testing and debugging purposes.
   *
   * @param[in] out Log stream for outputting the information.
   */
  void print_info(LogStream &out) const;

  /**
   * @brief Checks if the container is empty.
   *
   * @return True if the container is empty, false elsewhere.
   */
  bool is_empty() const;

private:

  /** Container for the objects. */
  ObjectMapTypes_ objects_;


#ifdef SERIALIZATION
  /**
   * @name Functions needed for serialization
   * @see <a href="http://uscilab.github.io/cereal/index.html">Cereal serialization library</a>
   */
  ///@{
  friend class cereal::access;

  template<class Archive>
  void serialize(Archive &ar);

  ///@}
#endif // SERIALIZATION

};


IGA_NAMESPACE_CLOSE

#ifdef SERIALIZATION
#include <igatools/base/objects_container.serial>
#endif // SERIALIZATION

#endif // XML_IO

#endif /*__OBJECTS_CONTAINER_H_ */
