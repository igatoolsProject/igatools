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

#ifndef __OBJECTS_CONTAINER_WRITER_H_
#define __OBJECTS_CONTAINER_WRITER_H_

#include <igatools/base/config.h>

#include <set>

#ifdef XML_IO


IGA_NAMESPACE_OPEN

class XMLElement;
class ObjectsContainer;
//template <int dim> class Grid;
//template <class T> class SafeSTLVector;
//class IgCoefficients;
//class LogStream;
//template <class T1, class T2> class SafeSTLMap;

/**
 * @brief Helper class for creating an @ref ObjectsContainer parsed from
 * a file.
 *
 * It receives an input file that is parsed by means of
 * the @ref XMLFileParser for receiving a @ref XMLElement containing
 * all the information.
 *
 * Before starting to parse the @ref XMLElement, the XML document
 * is validated against an XML schema definition. This schema
 * definition controls the main structure of the file: which elements
 * should be present, attributes, types, enumerations, etc.
 * But the schema is unable to control other questions, e.g. if a certain
 * @ref Grid, used by a @ref SplineSpace, is defined in the file.
 * This schema is stored as a @ref std::string in the private member
 * @ref XML_SCHEMA_.
 *
 * Currently, this class is able to parse the following classes:
 * - @ref Grid
 * - @ref SplineSpace
 * - @ref BSpline
 * - @ref NURBS
 * - @ref grid_functions::IdentityGridFunction
 * - @ref grid_functions::ConstantGridFunction
 * - @ref IgGridFunction
 * - @ref Domain
 * - @ref PhysicalSpaceBasis
 * - @ref functions::ConstantFunction
 * - @ref IgFunction
 *
 * If any type different from the ones above is found, an exception
 * is thrown.
 *
 * @warning It could be unsafe to take a @ref XMLElement outside of
 * the class in order to be reused somewhere else. It could cause problems
 * with the @p Xerces-c memory deallocation. Not getter should be
 * provided for retrieving it.
 *
 * @see ObjectsContainer
 * @see XMLFileParser
 * @see XMLElement
 *
 * @author P. Antolin
 * @date 2015
 */
class ObjectsContainerWriter
{
private:

  /** @name Types and static values */
  ///@{

  /** Type for the current class. */
  typedef ObjectsContainerWriter Self_;

  /** Type for a shared pointer of the current class. */
  typedef std::shared_ptr<Self_> SelfPtr_;

//  /** Type for mapping between @p Id numbers local to the input file and
//   *  object unique @p Id. */
//  typedef SafeSTLMap<Index, Index> IdMap_;

//  /**
//   * Statically defined string defining the XML schema for validating the
//   * input files.
//   */
//  static const std::string XML_SCHEMA_;

  ///@}

  /** @name Constructors, destructor, assignment operators and creators */
  ///@{

  /**
   * @brief Default constructor.
   * @note Deleted, not allowed to be used.
   */
  ObjectsContainerWriter() = delete;

  /**
   * @brief Copy constructor.
   * @note Deleted, not allowed to be used.
   */
  ObjectsContainerWriter(const ObjectsContainerWriter &) = delete;

  /**
   * @brief Move constructor.
   * @note Deleted, not allowed to be used.
   */
  ObjectsContainerWriter(ObjectsContainerWriter &&) = delete;

  /**
   * @brief Copy assignment operator.
   * @note Deleted, not allowed to be used.
   */
  ObjectsContainerWriter &operator= (const ObjectsContainerWriter &) = delete;

  /**
   * @brief Move assignment operator.
   * @note Deleted, not allowed to be used.
   */
  ObjectsContainerWriter &operator= (ObjectsContainerWriter &&) = delete;

  ///@}

public:
  /**
   * @brief Creates an @ref ObjectsContainer by reading the objects
   * from a given @p file_path and inserting them.
   *
   * @param[in] file_path Path of the file to be parsed.
   * @return Objects container with all the objects stored inside.
   */
  static void write(const std::string &file_path,
                    const std::shared_ptr<ObjectsContainer> container);

private:
  static void write_grids (const std::shared_ptr<ObjectsContainer> container,
                           XMLElement &xml_elem);

  static void write_spline_spaces (const std::shared_ptr<ObjectsContainer> container,
                                   XMLElement &xml_elem);

  static void write_reference_space_bases (const std::shared_ptr<ObjectsContainer> container,
                                     XMLElement &xml_elem);

  static void write_grid_functions (const std::shared_ptr<ObjectsContainer> container,
                                    XMLElement &xml_elem);

  static void write_domains (const std::shared_ptr<ObjectsContainer> container,
                             XMLElement &xml_elem);

  static void write_physical_space_bases (const std::shared_ptr<ObjectsContainer> container,
                                          XMLElement &xml_elem);

  static void write_functions (const std::shared_ptr<ObjectsContainer> container,
                               XMLElement &xml_elem);

  template <class Grid>
  static void write_grid (const std::shared_ptr<Grid> grid,
                          XMLElement &xml_elem);

  template <class SpSpace>
  static void write_spline_space (const std::shared_ptr<SpSpace> spline_space,
                                  XMLElement &xml_elem);

  template <class BSpline>
  static void write_bspline (const std::shared_ptr<BSpline> bspline,
                             XMLElement &xml_elem);

  template <class NURBS>
  static void write_nurbs (const std::shared_ptr<NURBS> nurbs,
                           XMLElement &xml_elem);

  template <class IdGridFunc>
  static void write_identity_grid_function (const std::shared_ptr<IdGridFunc> id_func,
                                            XMLElement &xml_elem);

  template <class ConstGridFunc>
  static void write_constant_grid_function (const std::shared_ptr<ConstGridFunc> const_func,
                                            XMLElement &xml_elem);

  template <class IgGridFunc>
  static void write_ig_grid_function (const std::shared_ptr<IgGridFunc> ig_func,
                                      XMLElement &xml_elem);

  template <class LinearGridFunc>
  static void write_linear_grid_function (const std::shared_ptr<LinearGridFunc> linear_func,
                                          XMLElement &xml_elem);

  template <class Domain>
  static void write_domain (const std::shared_ptr<Domain> domain,
                            XMLElement &xml_elem);

  template <class PhysSpaceBasis>
  static void write_phys_space_basis (const std::shared_ptr<PhysSpaceBasis> phys_space,
                                      XMLElement &xml_elem);

  template <class IgFunction>
  static void write_ig_function (const std::shared_ptr<IgFunction> ig_function,
                                 XMLElement &xml_elem);

  template <class ConstantFunction>
  static void write_constant_function (const std::shared_ptr<ConstantFunction> const_function,
                                       XMLElement &xml_elem);

  template <class LinearFunc>
  static void write_linear_function (const std::shared_ptr<LinearFunc> linear_func,
                                     XMLElement &xml_elem);

};

IGA_NAMESPACE_CLOSE

#endif // XML_IO

#endif // __OBJECTS_CONTAINER_WRITER_H_
