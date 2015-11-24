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

#ifdef XML_IO


IGA_NAMESPACE_OPEN

class XMLFileParser;
class XMLElement;
class ObjectsContainer;
template <int dim> class Grid;

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

public:
  /**
   * @brief Constructor taking the error handler and the file path.
   *
   * Constructor taking the error handler and the path of the file to be
   * parsed.
   *
   * Inside the constructor, the validity of the input file is checked by means
   * of @ref check_file, the XML platform is initialized and @ref parser_ is
   * created.
   *
   * @param[in] file_path Path of the file to be parsed.
   */
  ObjectsContainerParser(const std::string &file_path);

private:
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


public:

  /**
   * @brief Returns a new instance wrapped into a shared pointer.
   *
   * Builds and returns a new instance of the class wrapped into a
   * shared pointer. It uses the above defined default constructor.
   *
   * @param[in] file_path Path of the file to be parsed.
   * @return A shared pointer with a new instance of the class.
   */
  static SelfPtr_ create(const std::string &file_path);

  ///@}


private:

  /// Object to parse the file.
  /** Object to parse the file. */
  const std::shared_ptr<XMLFileParser> file_parser_;

public:
  /**
   * @todo To be documented.
   */
  std::shared_ptr<ObjectsContainer> parse(const std::string &schema_file) const;

private:
  void parse_grids(const std::shared_ptr<XMLElement> xml_elem,
                   const std::shared_ptr<ObjectsContainer> container) const;

  template <int dim>
  void parse_grid(const std::shared_ptr<XMLElement> xml_elem,
                  const std::shared_ptr<ObjectsContainer> container) const;

  void parse_spline_spaces(const std::shared_ptr<XMLElement> xml_elem,
                           const std::shared_ptr<ObjectsContainer> container) const;

  template <int dim, int range, int rank>
  void parse_spline_space(const std::shared_ptr<XMLElement> xml_elem,
                          const std::shared_ptr<ObjectsContainer> container) const;

  void parse_bsplines(const std::shared_ptr<XMLElement> xml_elem,
                      const std::shared_ptr<ObjectsContainer> container) const;

  template <int dim, int range, int rank>
  void parse_bspline(const std::shared_ptr<XMLElement> xml_elem,
                     const std::shared_ptr<ObjectsContainer> container) const;


};

IGA_NAMESPACE_CLOSE

#endif // XML_IO

#endif // __OBJECTS_CONTAINER_PARSER_H_
