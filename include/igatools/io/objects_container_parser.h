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

/**
 * @brief Class for parsing input files.
 *
 * This is a class for parsing XML input files and validate them against
 * a XML Schema grammar.
 *
 * This class provides the capability of creating a @p XercesDOMParser,
 * checking the validity of the input file (if it exists, it is corrupted,
 * etc), and retrieving a XML @p DOMDocument containing an input file.
 *
 * This class uses a @ref XMLParserErrorHandler for managing the possible
 * errors than can appear during the parsing process.
 *
 * The destructor of the class is in charge of of deleting @ref parser_
 * and to shutdown all the @p xerces active process.
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
   * @todo To be documented.
   */
  static std::shared_ptr<ObjectsContainer> parse(const std::string &file_path,
                                                 const std::string &schema_file);

private:
  /**
   * @todo To be documented.
   * @param[in] xml_elem XML element to be parsed.
   * @param[in,out] container Objects container to be filled with the parsed object.
   */
  static void parse_grids(const std::shared_ptr<XMLElement> xml_elem,
                   const std::shared_ptr<ObjectsContainer> container);

  /**
   * @todo To be documented.
   * @param[in] xml_elem XML element to be parsed.
   * @param[in,out] container Objects container to be filled with the parsed object.
   */
  static void parse_spline_spaces(const std::shared_ptr<XMLElement> xml_elem,
                           const std::shared_ptr<ObjectsContainer> container);

  /**
   * @todo To be documented.
   * @param[in] xml_elem XML element to be parsed.
   * @param[in,out] container Objects container to be filled with the parsed object.
   */
  static void parse_bsplines(const std::shared_ptr<XMLElement> xml_elem,
                      const std::shared_ptr<ObjectsContainer> container);

  /**
   * @todo To be documented.
   * @param[in] xml_elem XML element to be parsed.
   * @param[in,out] container Objects container to be filled with the parsed object.
   */
  static void parse_nurbs(const std::shared_ptr<XMLElement> xml_elem,
                   const std::shared_ptr<ObjectsContainer> container);

  /**
   * @todo To be documented.
   * @param[in] xml_elem XML element to be parsed.
   * @param[in,out] container Objects container to be filled with the parsed object.
   */
  static void parse_grid_functions_and_nurbs(const std::shared_ptr<XMLElement> xml_elem,
                           const std::shared_ptr<ObjectsContainer> container);

  /**
   * @todo To be documented.
   * @param[in] xml_elem XML element to be parsed.
   * @param[in,out] container Objects container to be filled with the parsed object.
   */
  static void parse_ig_grid_functions(const std::shared_ptr<XMLElement> xml_elem,
                               const std::shared_ptr<ObjectsContainer> container,
                               const bool &based_on_nurbs);

  /**
   * @todo To be documented.
   * @param[in] xml_elem XML element to be parsed.
   * @param[in,out] container Objects container to be filled with the parsed object.
   */
  static void parse_domains(const std::shared_ptr<XMLElement> xml_elem,
                     const std::shared_ptr<ObjectsContainer> container);

  /**
   * @todo To be documented.
   * @param[in] xml_elem XML element to be parsed.
   * @param[in,out] container Objects container to be filled with the parsed object.
   */
  static void parse_phys_spaces(const std::shared_ptr<XMLElement> xml_elem,
                         const std::shared_ptr<ObjectsContainer> container);

  /**
   * @todo To be documented.
   * @param[in] xml_elem XML element to be parsed.
   * @param[in,out] container Objects container to be filled with the parsed object.
   */
  static void parse_functions(const std::shared_ptr<XMLElement> xml_elem,
                       const std::shared_ptr<ObjectsContainer> container);

  /**
   * @todo To be documented.
   * @param[in] xml_elem XML element to be parsed.
   * @param[in,out] container Objects container to be filled with the parsed object.
   */
  template <int dim>
  static void parse_grid(const std::shared_ptr<XMLElement> xml_elem,
                  const std::shared_ptr<ObjectsContainer> container);

  /**
   * @todo To be documented.
   * @param[in] xml_elem XML element to be parsed.
   * @param[in,out] container Objects container to be filled with the parsed object.
   */
  template <int dim, int range, int rank>
  static void parse_spline_space(const std::shared_ptr<XMLElement> xml_elem,
                          const std::shared_ptr<ObjectsContainer> container);

  /**
   * @todo To be documented.
   * @param[in] xml_elem XML element to be parsed.
   * @param[in,out] container Objects container to be filled with the parsed object.
   */
  template <int dim, int range, int rank>
  static void parse_bspline(const std::shared_ptr<XMLElement> xml_elem,
                     const std::shared_ptr<ObjectsContainer> container);

  /**
   * @todo To be documented.
   * @param[in] xml_elem XML element to be parsed.
   * @param[in,out] container Objects container to be filled with the parsed object.
   */
  template <int dim, int range, int rank>
  static void parse_nurbs(const std::shared_ptr<XMLElement> xml_elem,
                   const std::shared_ptr<ObjectsContainer> container);

  /**
   * @todo To be documented.
   * @param[in] xml_elem XML element to be parsed.
   * @param[in,out] container Objects container to be filled with the parsed object.
   */
  template <int dim, int space_dim>
  static void parse_ig_grid_function(const std::shared_ptr<XMLElement> xml_elem,
                              const std::shared_ptr<ObjectsContainer> container,
                              const bool &first_parsing);

  /**
   * @todo To be documented.
   * @param[in] xml_elem XML element to be parsed.
   * @param[in,out] container Objects container to be filled with the parsed object.
   */
  template <int dim, int codim>
  static void parse_domain(const std::shared_ptr<XMLElement> xml_elem,
                    const std::shared_ptr<ObjectsContainer> container);

  /**
   * @todo To be documented.
   * @param[in] xml_elem XML element to be parsed.
   * @param[in,out] container Objects container to be filled with the parsed object.
   */
  template <int dim, int codim, int range, int rank>
  static void parse_phys_space(const std::shared_ptr<XMLElement> xml_elem,
                        const std::shared_ptr<ObjectsContainer> container);

  /**
   * @todo To be documented.
   * @param[in] xml_elem XML element to be parsed.
   * @param[in,out] container Objects container to be filled with the parsed object.
   */
  template <int dim, int codim, int range, int rank>
  static void parse_ig_function(const std::shared_ptr<XMLElement> xml_elem,
                         const std::shared_ptr<ObjectsContainer> container);

  /**
   * @todo To be documented.
   * @param[in] xml_elem XML element to be parsed.
   * @param[in,out] container Objects container to be filled with the parsed object.
   */
  static std::string parse_name(const std::shared_ptr<XMLElement> xml_elem);

  /**
   * @todo To be documented.
   * @param[in] xml_elem XML element to be parsed.
   * @param[in,out] container Objects container to be filled with the parsed object.
   */
  static std::string parse_dofs_property(const std::shared_ptr<XMLElement> xml_elem);

  /**
   * @todo To be documented.
   * @param[in] xml_elem XML element to be parsed.
   * @param[in,out] container Objects container to be filled with the parsed object.
   */
  static std::string get_type_id_string(const std::string &object_type,
                                 const Index &object_id,
                                 const SafeSTLVector<int> &dims);

  /**
   * @todo To be documented.
   * @param[in] xml_elem XML element to be parsed.
   * @param[in,out] container Objects container to be filled with the parsed object.
   */
  static std::string get_type_dimensions_string(const std::string &object_type,
                                         const SafeSTLVector<int> &dims);

  /**
   * @todo To be documented.
   * @param[in] xml_elem XML element to be parsed.
   * @param[in,out] container Objects container to be filled with the parsed object.
   */
  // TODO: maybe this should be done with a shared pointer, in
  // order to astatic void the copy
  static IgCoefficients parse_ig_coefficients(const std::shared_ptr<XMLElement> xml_elem,
                                       const std::string &parsing_msg,
                                       const std::set<Index> &space_global_dofs);

};

IGA_NAMESPACE_CLOSE

#endif // XML_IO

#endif // __OBJECTS_CONTAINER_PARSER_H_
