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

#ifndef __XML_FILE_PARSER_H_
#define __XML_FILE_PARSER_H_

#include <igatools/base/config.h>

#ifdef XML_IO

#include <xercesc/util/XercesDefs.hpp>
#include <xercesc/sax/ErrorHandler.hpp>

XERCES_CPP_NAMESPACE_BEGIN
class XercesDOMParser;
class DOMDocument;
XERCES_CPP_NAMESPACE_END


IGA_NAMESPACE_OPEN

class XMLElement;
class XMLParserErrorHandler;
class LogStream;

/**
 * @brief Class for parsing XML input files.
 *
 * This class provides the capability of creating a @p XercesDOMParser,
 * checking the validity of the input file (if it exists, it is corrupted,
 * etc), and retrieving a  @ref XMLElement wrapping the XML document
 * contained in the input file.
 *
 * Optionally, the parser validates the read XML document against a XML
 * Schema grammar.
 *
 * Once the class is instantiated by means of the @ref create static
 * method, the file is parsed by calling the @ref parse method.
 * This @ref parse method can be called as many times as needed for
 * parsing more than one file.
 *
 * This class uses a @ref XMLParserErrorHandler for managing the possible
 * errors than can appear during the parsing process.
 *
 * The destructor of the class is in charge of of deleting @ref parser_
 * and to shutdown all the @p Xerces-c active process.
 *
 * @note This class uses the @p Xerces-c library.
 *
 * @see XMLParserErrorHandler
 *
 * @author P. Antolin
 * @date 2015
 */
class XMLFileParser
{
private:

  /** @name Types and static values */
  ///@{

  /** Type for the current class. */
  typedef XMLFileParser Self_;

  /** Type for a shared pointer of the current class. */
  typedef std::shared_ptr<Self_> SelfPtr_;

  ///@}

  /** @name Constructors, destructor, assignment operators and creators */
  ///@{

  /**
   * @brief Constructor.
   *
   * This constructor initializes the @p Xerces-c utils and creates
   * the XML parser.
   */
  XMLFileParser();

  /**
   * @brief Copy constructor.
   * @note Deleted, not allowed to be used.
   */
  XMLFileParser(const XMLFileParser &) = delete;

  /**
   * @brief Move constructor.
   * @note Deleted, not allowed to be used.
   */
  XMLFileParser(XMLFileParser &&) = delete;

  /**
   * @brief Copy assignment operator.
   * @note Deleted, not allowed to be used.
   */
  XMLFileParser &operator= (const XMLFileParser &) = delete;

  /**
   * @brief Move assignment operator.
   * @note Deleted, not allowed to be used.
   */
  XMLFileParser &operator= (XMLFileParser &&) = delete;


public:

  /**
   * @brief Destructor.
   *
   * It finalizes the parsing process: deletes the @ref parser_
   * and terminates the platform utils.
   */
  ~XMLFileParser();

  /**
   * @brief Build and returns a new instance wrapped into a shared pointer.
   *
   * It uses the above defined constructor.
   *
   * @return A shared pointer with a new instance of the class.
   */
  static SelfPtr_ create();

  ///@}

public:
  /**
   * @brief Parses the input file and returns a XML document object.
   *
   * Before parsing the file of the input file is checked by means
   * of the static method @ref check_file.
   *
   * This method validates the content of @p file_path against the
   * XML schema defined in the @p grammar_definition string. This is a
   * string containing the full definition of the grammar.
   *
   * @warning If there is any problem parsing the input file, error
   * messages and exceptions will be thrown.
   *
   * @param[in] file_path Path of the file to be parsed.
   * @param[in] grammar_definition Definition of the XML schema
   *                               contained into a string.
   * @return XML element containing the parsed document.
   */
  std::shared_ptr<XMLElement> parse(const std::string &file_path,
                                    const std::string &grammar_definition) const;

  /**
   * @brief Parses the input file and returns a XML document object.
   *
   * Parses the input file and returns a XML document object.
   * Before parsing the file the validity of the file is checked by
   * calling the static method @ref check_file.
   *
   * @warning If there is any problem parsing the input file, error
   * messages and exceptions will be thrown.
   *
   * @param[in] file_path Path of the file to be parsed.
   * @return XML element containing the parsed document.
   */
  std::shared_ptr<XMLElement> parse(const std::string &file_path) const;

private:

  /** Parser error handler. */
  const std::shared_ptr<XMLParserErrorHandler> error_handler_;

  /** DOM parser. */
  xercesc::XercesDOMParser *parser_;

  /**
   * @brief Checks if the file can be read.
   *
   * It throws an error message if the file can not be read properly.
   *
   * @param[in] file_path Path of the file to be checked.
   */
  static void check_file(const std::string &file_path);

};

IGA_NAMESPACE_CLOSE

#endif // XML_IO

#endif // __XML_FILE_PARSER_H_
