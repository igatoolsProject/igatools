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

#ifndef __OBJECTS_CONTAINER_PARSER_H_
#define __OBJECTS_CONTAINER_PARSER_H_

#include <igatools/base/config.h>

#include <set>

#ifdef XML_IO


IGA_NAMESPACE_OPEN

class XMLFileParser;
class XMLElement;
class ObjectsContainer;
template <int dim> class Grid;
template <class T> class SafeSTLVector;
class IgCoefficients;
class LogStream;
template <class T1, class T2> class SafeSTLMap;

/**
 * @brief Helper class for creating an @ref ObjectsContainer parsed from
 * a file.
 *
 * This is a helper class for creating an @ref ObjectsContainer parsed
 * from a file. It receives an input file that is parsed by means of
 * the @ref XMLFileParser for receiving a @ref XMLElement containing
 * all the information.
 *
 * Before starting to parse the @ref XMLElement, the XML document
 * is validated against an XML schema definition. This schema
 * definition controls the main structure of the file: which elements
 * should be present, attributes, types, enumerations, etc.
 * But the schema is unable to control other questions, e.g. if a certain
 * @p Grid, used by a @ref SplineSpace, is defined in the file.
 *
 * Currently, this class is able to parse:
 * - @ref Grid
 * - @ref SplineSpace
 * - @ref BSpline
 * - @ref NURBS
 * - @ref IdentityGridFunction
 * - @ref ConstantGridFunction
 * - @ref IgGridFunction
 * - @ref Domain
 * - @ref PhysicalSpaceBasis
 * - @ref ConstantFunction
 * - @ref IgFunction
 *
 * If any type different from the ones above is found an exception
 * is thrown.
 *
 * @alert It could be unsafe to take a @ref XMLElement outside of
 * the class in order to be reused somewhere else. It could cause problem
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
class ObjectsContainerParser
{
private:

  /** @name Types and static values */
  ///@{

  /** Type for the current class. */
  typedef ObjectsContainerParser Self_;

  /** Type for a shared pointer of the current class. */
  typedef std::shared_ptr<Self_> SelfPtr_;

  /** Type for map between file local Id number and object unique Id. */
  typedef SafeSTLMap<Index, Index> IdMap_;

  /**
   * Statically defined string defining the XML schema for validating the
   * input files.
   */
  static const std::string XML_SCHEMA_;

  ///@}

  /** @name Constructors, destructor, assignment operators and creators */
  ///@{

  /**
   * @brief Deleted default constructor.
   *
   * Default constructor.
   * @note Deleted: not allowed.
   */
  ObjectsContainerParser() = delete;

  /**
   * @brief Deleted copy constructor.
   *
   * Copy constructor.
   * @note Deleted: not allowed.
   */
  ObjectsContainerParser(const ObjectsContainerParser &) = delete;

  /**
   * @brief Deleted move constructor.
   *
   * Move constructor.
   * @note Deleted: not allowed.
   */
  ObjectsContainerParser(ObjectsContainerParser &&) = delete;

  /**
   * @brief Deleted copy assignment operator.
   *
   * Copy assignment operator.
   * @note Deleted: not allowed.
   */
  ObjectsContainerParser &operator= (const ObjectsContainerParser &) = delete;

  /**
   * @brief Deleted move assignment operator.
   *
   * Move assignment operator.
   * @note Deleted: not allowed.
   */
  ObjectsContainerParser &operator= (ObjectsContainerParser &&) = delete;

  ///@}

public:
  /**
   * @brief Creates an @ref ObjectsContainer by reading the objects
   * from a given file.
   *
   * Creates an @ref ObjectsContainer by reading the objects
   * from a given @ref file_path and inserting them.
   *
   * @param[in] file_path Path of the file to be parsed.
   * @return Objects container with all the objects stored inside.
   */
  static std::shared_ptr<ObjectsContainer> parse(const std::string &file_path);

  /**
   * @brief Creates an @ref ObjectsContainer by reading the objects
   * from a given file and creating and storing them as constant.
   *
   * Creates an @ref ObjectsContainer by reading the objects
   * from a given @ref file_path and creating and inserting them as
   * constant.
   *
   * @param[in] file_path Path of the file to be parsed.
   * @return Objects container with all the objects stored inside.
   */
  static std::shared_ptr<ObjectsContainer> parse_const(const std::string &file_path);

private:
  /**
   * @brief Parses all the <tt>Grid</tt>s
   * contained into the XML document and stores them into the container.
   *
   * Parses all the <tt>Grid</tt>s
   * contained into the XML document @ref xml_elem and stores them into
   * the @ref container.
   *
   * @param[in] xml_elem XML element to be parsed.
   * @param[in,out] container Container for inserting the objects
   *                and also retrieving other ones needed.
   */
  static void parse_grids(const std::shared_ptr<XMLElement> xml_elem,
                          const bool parse_as_constant,
                          IdMap_ &id_map,
                          const std::shared_ptr<ObjectsContainer> container);

  /**
   * @brief Parses all the <tt>SplineSpace</tt>s
   * contained into the XML document and stores them into the container.
   *
   * Parses all the <tt>SplineSpace</tt>s
   * contained into the XML document @ref xml_elem and stores them into
   * the @ref container.
   *
   * @param[in] xml_elem XML element to be parsed.
   * @param[in,out] Map between the local Ids of the input file and
   *                the unique object Id.
   * @param[in,out] container Container for inserting the objects
   *                and also retrieving other ones needed.
   */
  static void parse_spline_spaces(const std::shared_ptr<XMLElement> xml_elem,
                                  const bool parse_as_constant,
                                  IdMap_ &id_map,
                                  const std::shared_ptr<ObjectsContainer> container);

  /**
   * @brief Parses all the <tt>BSpline</tt>s
   * contained into the XML document and stores them into the container.
   *
   * Parses all the <tt>BSpline</tt>s
   * contained into the XML document @ref xml_elem and stores them into
   * the @ref container.
   *
   * @param[in] xml_elem XML element to be parsed.
   * @param[in,out] Map between the local Ids of the input file and
   *                the unique object Id.
   * @param[in,out] container Container for inserting the objects
   *                and also retrieving other ones needed.
   */
  static void parse_bsplines(const std::shared_ptr<XMLElement> xml_elem,
                             const bool parse_as_constant,
                             IdMap_ &id_map,
                             const std::shared_ptr<ObjectsContainer> container);

  /**
   * @brief Parses all the <tt>NURBS</tt>s
   * contained into the XML document and stores them into the container.
   *
   * Parses all the <tt>NURBS</tt>s
   * contained into the XML document @ref xml_elem and stores them into
   * the @ref container.
   *
   * @param[in] xml_elem XML element to be parsed.
   * @param[in,out] Map between the local Ids of the input file and
   *                the unique object Id.
   * @param[in,out] container Container for inserting the objects
   *                and also retrieving other ones needed.
   */
  static void parse_nurbs(const std::shared_ptr<XMLElement> xml_elem,
                          const bool parse_as_constant,
                          IdMap_ &id_map,
                          const std::shared_ptr<ObjectsContainer> container);

  /**
   * @brief Parses all the <tt>GridFunction</tt>s and <tt>NURBS</tt>s
   * contained into the XML document and stores them into the container.
   *
   * Parses all the <tt>GridFunction</tt>s and <tt>NURBS</tt>s
   * contained into the XML document @ref xml_elem and stores them into
   * the @ref container.
   *
   * They are parsed in proper order to avoid conflicts in the declaration
   * of <tt>NURBS</tt>s that take <tt>IgGridFunction</tt>s as arguments.
   *
   * @param[in] xml_elem XML element to be parsed.
   * @param[in,out] Map between the local Ids of the input file and
   *                the unique object Id.
   * @param[in,out] container Container for inserting the objects
   *                and also retrieving other ones needed.
   */
  static void parse_grid_functions_and_nurbs(const std::shared_ptr<XMLElement> xml_elem,
                                             const bool parse_as_constant,
                                             IdMap_ &id_map,
                                             const std::shared_ptr<ObjectsContainer> container);

  /**
   * @brief Parses all the <tt>IdentityGridFunction</tt>s
   * contained into the XML document and stores them into the container.
   *
   * @param[in] xml_elem XML element to be parsed.
   * @param[in,out] Map between the local Ids of the input file and
   *                the unique object Id.
   * @param[in,out] container Container for inserting the objects
   *                and also retrieving other ones needed.
   */
  static void parse_identity_grid_functions(const std::shared_ptr<XMLElement> xml_elem,
                                            const bool parse_as_constant,
                                            IdMap_ &id_map,
                                            const std::shared_ptr<ObjectsContainer> container);

  /**
   * @brief Parses all the <tt>ConstantGridFunction</tt>s
   * contained into the XML document and stores them into the container.
   *
   * @param[in] xml_elem XML element to be parsed.
   * @param[in,out] Map between the local Ids of the input file and
   *                the unique object Id.
   * @param[in,out] container Container for inserting the objects
   *                and also retrieving other ones needed.
   */
  static void parse_constant_grid_functions(const std::shared_ptr<XMLElement> xml_elem,
                                            const bool parse_as_constant,
                                            IdMap_ &id_map,
                                            const std::shared_ptr<ObjectsContainer> container);

  /**
   * @brief Parses all the <tt>IgGridFunction</tt>s
   * contained into the XML document and stores them into the container.
   *
   * Parses all the <tt>IgGridFunction</tt>s
   * contained into the XML document @ref xml_elem and stores them into
   * the @ref container.
   *
   * @param[in] xml_elem XML element to be parsed.
   * @param[in,out] Map between the local Ids of the input file and
   *                the unique object Id.
   * @param[in,out] container Container for inserting the objects
   *                and also retrieving other ones needed.
   */
  static void parse_ig_grid_functions(const std::shared_ptr<XMLElement> xml_elem,
                                      const bool parse_as_constant,
                                      const bool &first_parsing,
                                      IdMap_ &id_map,
                                      const std::shared_ptr<ObjectsContainer> container);

  /**
   * @brief Parses all the <tt>Domain</tt>s
   * contained into the XML document and stores them into the container.
   *
   * Parses all the <tt>Domain</tt>s
   * contained into the XML document @ref xml_elem and stores them into
   * the @ref container.
   *
   * @param[in] xml_elem XML element to be parsed.
   * @param[in,out] Map between the local Ids of the input file and
   *                the unique object Id.
   * @param[in,out] container Container for inserting the objects
   *                and also retrieving other ones needed.
   */
  static void parse_domains(const std::shared_ptr<XMLElement> xml_elem,
                            const bool parse_as_constant,
                            IdMap_ &id_map,
                            const std::shared_ptr<ObjectsContainer> container);

  /**
   * @brief Parses all the <tt>PhysicalSpaceBasis</tt>s
   * contained into the XML document and stores them into the container.
   *
   * Parses all the <tt>PhysicalSpaceBasis</tt>s
   * contained into the XML document @ref xml_elem and stores them into
   * the @ref container.
   *
   * @param[in] xml_elem XML element to be parsed.
   * @param[in,out] Map between the local Ids of the input file and
   *                the unique object Id.
   * @param[in,out] container Container for inserting the objects
   *                and also retrieving other ones needed.
   */
  static void parse_phys_spaces(const std::shared_ptr<XMLElement> xml_elem,
                                const bool parse_as_constant,
                                IdMap_ &id_map,
                                const std::shared_ptr<ObjectsContainer> container);

  /**
   * @brief Parses all the <tt>Function</tt>s
   * contained into the XML document and stores them into the container.
   *
   * Parses all the <tt>Function</tt>s
   * contained into the XML document @ref xml_elem and stores them into
   * the @ref container.
   *
   * @param[in] xml_elem XML element to be parsed.
   * @param[in,out] Map between the local Ids of the input file and
   *                the unique object Id.
   * @param[in,out] container Container for inserting the objects
   *                and also retrieving other ones needed.
   */
  static void parse_functions(const std::shared_ptr<XMLElement> xml_elem,
                              const bool parse_as_constant,
                              IdMap_ &id_map,
                              const std::shared_ptr<ObjectsContainer> container);

  /**
   * @brief Parses all the <tt>IgFunction</tt>s
   * contained into the XML document and stores them into the container.
   *
   * Parses all the <tt>IgFunction</tt>s
   * contained into the XML document @ref xml_elem and stores them into
   * the @ref container.
   *
   * @param[in] xml_elem XML element to be parsed.
   * @param[in,out] Map between the local Ids of the input file and
   *                the unique object Id.
   * @param[in,out] container Container for inserting the objects
   *                and also retrieving other ones needed.
   */
  static void parse_ig_functions(const std::shared_ptr<XMLElement> xml_elem,
                                 const bool parse_as_constant,
                                 IdMap_ &id_map,
                                 const std::shared_ptr<ObjectsContainer> container);

  /**
   * @brief Parses all the <tt>ConstantFunction</tt>s
   * contained into the XML document and stores them into the container.
   *
   * Parses all the <tt>ConstantFunction</tt>s
   * contained into the XML document @ref xml_elem and stores them into
   * the @ref container.
   *
   * @param[in] xml_elem XML element to be parsed.
   * @param[in,out] Map between the local Ids of the input file and
   *                the unique object Id.
   * @param[in,out] container Container for inserting the objects
   *                and also retrieving other ones needed.
   */
  static void parse_constant_functions(const std::shared_ptr<XMLElement> xml_elem,
                                       const bool parse_as_constant,
                                       IdMap_ &id_map,
                                       const std::shared_ptr<ObjectsContainer> container);

  /**
   * @brief Parses a <tt>Grid</tt> XML element and inserts
   * it into the objects container.
   *
   * Parses a <tt>Grid</tt> XML element contained in
   * @ref xml_elem and inserts it into the objects @ref container.
   *
   * @tparam dim Dimension of the Grid.
   * @param[in] xml_elem XML element to be parsed.
   * @param[in,out] Map between the local Ids of the input file and
   *                the unique object Id.
   * @param[in,out] container Container for inserting the object
   *                and also retrieving other ones needed.
   */
  template <int dim>
  static void parse_grid(const std::shared_ptr<XMLElement> xml_elem,
                         const bool parse_as_constant,
                         IdMap_ &id_map,
                         const std::shared_ptr<ObjectsContainer> container);

  /**
   * @brief Parses a <tt>SplineSpace</tt> XML element and inserts
   * it into the objects container.
   *
   * Parses a <tt>SplineSpace</tt> XML element contained in
   * @ref xml_elem and inserts it into the objects @ref container.
   *
   * @tparam dim Dimension of the SplineSpace.
   * @tparam range Range of the SplineSpace.
   * @tparam rank Rank of the SplineSpace.
   * @param[in] xml_elem XML element to be parsed.
   * @param[in,out] Map between the local Ids of the input file and
   *                the unique object Id.
   * @param[in,out] container Container for inserting the object
   *                and also retrieving other ones needed.
   */
  template <int dim, int range, int rank>
  static void parse_spline_space(const std::shared_ptr<XMLElement> xml_elem,
                                 const bool parse_as_constant,
                                 IdMap_ &id_map,
                                 const std::shared_ptr<ObjectsContainer> container);

  /**
   * @brief Parses a <tt>BSpline</tt> XML element and inserts
   * it into the objects container.
   *
   * Parses a <tt>BSpline</tt> XML element contained in
   * @ref xml_elem and inserts it into the objects @ref container.
   *
   * @tparam dim Dimension of the BSpline space basis.
   * @tparam range Range of the BSpline space basis.
   * @tparam rank Rank of the BSpline space basis.
   * @param[in] xml_elem XML element to be parsed.
   * @param[in,out] Map between the local Ids of the input file and
   *                the unique object Id.
   * @param[in,out] container Container for inserting the object
   *                and also retrieving other ones needed.
   */
  template <int dim, int range, int rank>
  static void parse_bspline(const std::shared_ptr<XMLElement> xml_elem,
                            const bool parse_as_constant,
                            IdMap_ &id_map,
                            const std::shared_ptr<ObjectsContainer> container);

  /**
   * @brief Parses a <tt>NURBS</tt> XML element and inserts
   * it into the objects container.
   *
   * Parses a <tt>NURBS</tt> XML element contained in
   * @ref xml_elem and inserts it into the objects @ref container.
   *
   * @tparam dim Dimension of the NURBS space basis.
   * @tparam range Range of the NURBS space basis.
   * @tparam rank Rank of the NURBS space basis.
   * @param[in] xml_elem XML element to be parsed.
   * @param[in,out] Map between the local Ids of the input file and
   *                the unique object Id.
   * @param[in,out] container Container for inserting the object
   *                and also retrieving other ones needed.
   */
  template <int dim, int range, int rank>
  static void parse_nurbs(const std::shared_ptr<XMLElement> xml_elem,
                          const bool parse_as_constant,
                          IdMap_ &id_map,
                          const std::shared_ptr<ObjectsContainer> container);

  /**
   * @brief Parses a <tt>IdentityGridFunction</tt> XML element and inserts
   * it into the objects container.
   *
   * Parses a <tt>IdentityGridFunction</tt> XML element contained in
   * @ref xml_elem and inserts it into the objects @ref container.
   *
   * @tparam dim Dimension of the grid function.
   * @param[in] xml_elem XML element to be parsed.
   * @param[in,out] Map between the local Ids of the input file and
   *                the unique object Id.
   * @param[in,out] container Container for inserting the object
   *                and also retrieving other ones needed.
   */
  template <int dim>
  static void parse_identity_grid_function(const std::shared_ptr<XMLElement> xml_elem,
                                           const bool parse_as_constant,
                                           IdMap_ &id_map,
                                           const std::shared_ptr<ObjectsContainer> container);

  /**
   * @brief Parses a <tt>ConstantGridFunction</tt> XML element and inserts
   * it into the objects container.
   *
   * Parses a <tt>ConstantGridFunction</tt> XML element contained in
   * @ref xml_elem and inserts it into the objects @ref container.
   *
   * @tparam dim Dimension of the grid function.
   * @tparam space_dim Space dimension of the grid function.
   * @param[in] xml_elem XML element to be parsed.
   * @param[in,out] Map between the local Ids of the input file and
   *                the unique object Id.
   * @param[in,out] container Container for inserting the object
   *                and also retrieving other ones needed.
   */
  template <int dim, int space_dim>
  static void parse_constant_grid_function(const std::shared_ptr<XMLElement> xml_elem,
                                           const bool parse_as_constant,
                                           IdMap_ &id_map,
                                           const std::shared_ptr<ObjectsContainer> container);

  /**
   * @brief Parses a <tt>IgGridFunction</tt> XML element and inserts
   * it into the objects container.
   *
   * Parses a <tt>IgGridFunction</tt> XML element contained in
   * @ref xml_elem and inserts it into the objects @ref container.
   *
   * This function is called twice:
   * - During the first time (<tt>first_parsing == true</tt>) the
   *   grid functions based on @p BSpline are parsed.
   * - In the second time (<tt>first_parsing == false</tt>) the
   *   grid functions based on @p NURBS are parsed.
   *
   * @tparam dim Dimension of the grid function.
   * @tparam space_dim Space dimension of the grid function.
   * @param[in] xml_elem XML element to be parsed.
   * @param[in,out] Map between the local Ids of the input file and
   *                the unique object Id.
   * @param[in] first_parsing Indicates if the function is called for
   *                          first time, or not.
   * @param[in,out] container Container for inserting the object
   *                and also retrieving other ones needed.
   */
  template <int dim, int space_dim>
  static void parse_ig_grid_function(const std::shared_ptr<XMLElement> xml_elem,
                                     const bool parse_as_constant,
                                     const bool &first_parsing,
                                     IdMap_ &id_map,
                                     const std::shared_ptr<ObjectsContainer> container);

  /**
   * @brief Parses a <tt>Domain</tt> XML element and inserts
   * it into the objects container.
   *
   * Parses a <tt>Domain</tt> XML element contained in
   * @ref xml_elem and inserts it into the objects @ref container.
   *
   * @tparam dim Dimension of the domain.
   * @tparam codim Codimension of the domain.
   * @param[in] xml_elem XML element to be parsed.
   * @param[in,out] Map between the local Ids of the input file and
   *                the unique object Id.
   * @param[in,out] container Container for inserting the object
   *                and also retrieving other ones needed.
   */
  template <int dim, int codim>
  static void parse_domain(const std::shared_ptr<XMLElement> xml_elem,
                           const bool parse_as_constant,
                           IdMap_ &id_map,
                           const std::shared_ptr<ObjectsContainer> container);

  /**
   * @brief Parses a <tt>PhysicalSpaceBasis</tt> XML element and inserts
   * it into the objects container.
   *
   * Parses a <tt>PhysicalSpaceBasis</tt> XML element contained in
   * @ref xml_elem and inserts it into the objects @ref container.
   *
   * @tparam dim Dimension of the physical space basis.
   * @tparam range Range of the physical space basis.
   * @tparam rank Rank of the physical space basis.
   * @tparam codim Codimension of the physical space basis.
   * @param[in] xml_elem XML element to be parsed.
   * @param[in,out] Map between the local Ids of the input file and
   *                the unique object Id.
   * @param[in,out] container Container for inserting the object
   *                and also retrieving other ones needed.
   */
  template <int dim, int codim, int range, int rank>
  static void parse_phys_space(const std::shared_ptr<XMLElement> xml_elem,
                               const bool parse_as_constant,
                               IdMap_ &id_map,
                               const std::shared_ptr<ObjectsContainer> container);

  /**
   * @brief Parses an <tt>IgFunction</tt> XML element and inserts it into
   * the objects container.
   *
   * Parses an <tt>IgFunction</tt> XML element contained in @ref xml_elem
   * and inserts it into the objects @ref container.
   *
   * @tparam dim Dimension of the function.
   * @tparam codim Codimension of the function.
   * @tparam range Range of the function.
   * @tparam rank Rank of the function.
   * @param[in] xml_elem XML element to be parsed.
   * @param[in,out] Map between the local Ids of the input file and
   *                the unique object Id.
   * @param[in,out] container Container for inserting the object
   *                and also retrieving other ones needed.
   */
  template <int dim, int codim, int range, int rank>
  static void parse_ig_function(const std::shared_ptr<XMLElement> xml_elem,
                                const bool parse_as_constant,
                                IdMap_ &id_map,
                                const std::shared_ptr<ObjectsContainer> container);

  /**
   * @brief Parses an <tt>ConstantFunction</tt> XML element and inserts it into
   * the objects container.
   *
   * Parses an <tt>ConstantFunction</tt> XML element contained in @ref xml_elem
   * and inserts it into the objects @ref container.
   *
   * @tparam dim Dimension of the function.
   * @tparam codim Codimension of the function.
   * @tparam range Range of the function.
   * @tparam rank Rank of the function.
   * @param[in] xml_elem XML element to be parsed.
   * @param[in,out] Map between the local Ids of the input file and
   *                the unique object Id.
   * @param[in,out] container Container for inserting the object
   *                and also retrieving other ones needed.
   */
  template <int dim, int codim, int range, int rank>
  static void parse_constant_function(const std::shared_ptr<XMLElement> xml_elem,
                                      const bool parse_as_constant,
                                      IdMap_ &id_map,
                                      const std::shared_ptr<ObjectsContainer> container);

  /**
   * @brief Parses a <tt>Name</tt> XML element.
   * Parses a <tt>Name</tt> XML element contained in @ref xml_elem.
   *
   * @param[in] xml_elem XML element to be parsed.
   * return Name string parsed.
   */
  static std::string parse_name(const std::shared_ptr<XMLElement> xml_elem);

  /**
   * @brief Parses a <tt>DofsProperty</tt> XML element.
   * Parses a <tt>DofsProperty</tt> XML element contained in @ref xml_elem.
   *
   * @param[in] xml_elem XML element to be parsed.
   * return Dofs property string parsed.
   */
  static std::string parse_dofs_property(const std::shared_ptr<XMLElement> xml_elem);

  /**
   * @brief Produces an string of the given type and dimensions and the
   * object id associated to it.
   *
   * Produces an string of the given type and dimensions and the object
   * id associated to it.
   *
   * It will produce a string like <tt>Type<X, Y, Z> (IgaObjectId (W))</tt>.
   * This is used for error messaging purposes.
   *
   * @param[in] object_type String with the object type.
   * @param[in] object_id Id associated to the object.
   * @param[in] dims Vector with the dimensions of the object.
   * @return String containing the object type and its id
   */
  static std::string get_type_id_string(const std::string &object_type,
                                        const Index &object_id,
                                        const SafeSTLVector<int> &dims);

  /**
   * @todo Parses an IgCoefficients vector.
   *
   * Parses an IgCoefficients vector from the given @p xml_elem
   * The indices of the parsed coefficients are checked with
   * @p space_global_dofs, that are the global indices of the dofs
   * of the space associated to the coefficients.
   *
   * @param[in] xml_elem XML element to be parsed.
   * @param[in] parsing_msg String containing information about the
   *             XML element document to produce richer error messages.
   * @param[in] space_global_dofs Global indices of the dofs of
   *            space to which the coefficients are associated to.
   * @return IgCoefficients shared pointer vector.
   */
  // TODO: maybe this should be done with a shared pointer, in
  // order to avoid the copy.
  static std::shared_ptr<IgCoefficients>
  parse_ig_coefficients(const std::shared_ptr<XMLElement> xml_elem,
                        const std::string &parsing_msg,
                        const std::set<Index> &space_global_dofs);

};

IGA_NAMESPACE_CLOSE

#endif // XML_IO

#endif // __OBJECTS_CONTAINER_PARSER_H_
