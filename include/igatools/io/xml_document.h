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

#ifndef __XML_DOCUMENT_H_
#define __XML_DOCUMENT_H_

#include <igatools/base/config.h>

#ifdef XML_IO

#include <xercesc/util/XercesDefs.hpp>
#include <xercesc/dom/DOMDocument.hpp>

XERCES_CPP_NAMESPACE_BEGIN
class DOMDocument;
class DOMImplementation;
class DOMElement;
XERCES_CPP_NAMESPACE_END


IGA_NAMESPACE_OPEN

template <class T> class SafeSTLVector;
class LogStream;
class XMLElement;

/**
 * @brief Class for managing XML DOM elements of @p Xerces-c.
 *
 * This is wrapper class of the XML <tt>Xerces-c DOMElement</tt> class.
 *
 * It takes a pointer to a @p DOMElement as argument for the constructor
 * and provides some functionalities for accessing data and querying
 * information of the XML element.
 *
 * The class provides methods for extracting children elements
 * of the element, that are wrapped inside new instances of the class
 * and returned.
 *
 * @note This class uses @p Xerces-c library.
 * @note The class is not in charge of freeing the @p DOMElement pointer.
 *
 * @see XMLFileParser
 *
 * @author P. Antolin
 * @date 2015
 */
class XMLDocument
{
private:

  /** @name Types and static values */
  ///@{

  /** Type for the current class. */
  typedef XMLDocument Self_;

  /** Type for a shared pointer of the current class. */
  typedef std::shared_ptr<Self_> SelfPtr_;

  /** Type for a pointer of the @p Xerces-c DOM document class. */
  typedef xercesc::DOMDocument *DOMDocPtr_;

  ///@}

  /** @name Constructors, destructor, assignment operators and creators */
  ///@{

  /**
   * Constructor.
   *
   * Constructs a new instance of the class by taking the @p name
   * of the first XML element in the document.
   *
   * @param[in] name Name of the XML element.
   */
  XMLDocument(const char *name);

  XMLDocument(const std::string &file_path);

  XMLDocument(const std::string &file_path,
              const std::string &grammar_definition);

  /**
   * @brief Default constructor.
   * @note Deleted, not allowed to be used.
   */
  XMLDocument() = delete;

  /**
   * @brief Copy constructor.
   * @note Deleted, not allowed to be used.
   */
  XMLDocument(const XMLDocument &) = delete;

  /**
   * @brief Move constructor.
   * @note Deleted, not allowed to be used.
   */
  XMLDocument(XMLDocument &&) = delete;

  /**
   * @brief Copy assignment operator.
   * @note Deleted, not allowed to be used.
   */
  XMLDocument &operator= (const XMLDocument &) = delete;

  /**
   * @brief Move assignment operator.
   * @note Deleted, not allowed to be used.
   */
  XMLDocument &operator= (XMLDocument &&) = delete;

public:

  /**
   * @brief Destructor.
   *
   * It finalizes the parsing process: deletes @ref xml_doc_ and
   * @ref dom_impl_ and terminates the platform utils.
   */
  ~XMLDocument();

  /**
   * @brief Creates a new document with a single XML element with
   * the given @p name.
   *
   * @param[in] name Name of the XML element.
   * @return A shared pointer with a new instance of the class.
   */
  static SelfPtr_ create_void_document(const std::string &name);

  /**
   * @brief Parses the input file and returns a instance of the class
   * wrapped into a shared pointer.
   *
   * Before parsing the file the validity of the file is checked by
   * calling the static method @ref check_file.
   *
   * @warning If there is any problem parsing the input file, error
   * messages and exceptions will be thrown.
   *
   * @param[in] file_path Path of the file to be parsed.
   * @return XML document containing the parsed document.
   */
  static SelfPtr_ parse_from_file(const std::string &file_path);

  /**
   * @brief Parses the input file and returns a instance of the class
   * wrapped into a shared pointer.
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
   * @return XML document containing the parsed document.
   */
  static SelfPtr_ parse_from_file(const std::string &file_path,
                                  const std::string &grammar_definition);

  ///@}

  /**
   * @brief Retrieves the document element wrapped into a @ref XMLElement.
   *
   * @return DOM document element.
   */
  std::shared_ptr<XMLElement> get_document_element() const;

  /**
   * @brief Creates a new @ref XMLElement with the given name.
   *
   * @param[in] name Name of the new element created.
   * @return XML element.
   */
  std::shared_ptr<XMLElement> create_new_element(const std::string &name) const;

  /**
   * @brief Creates a new @ref XMLElement with the given @p name and
   * a the @p text.
   *
   * @param[in] name Name of the new element created.
   * @param[in] text Text contained in the element.
   * @return XML element.
   */
  std::shared_ptr<XMLElement> create_new_text_element(const std::string &name,
                                                      const std::string &text) const;

  /**
   * @brief Creates a new @ref XMLElement with the given name.
   *
   * @param[in] name Name of the new element created.
   * @return XML element.
   */
  template <class T>
  std::shared_ptr<XMLElement> create_size_dir_vector_element
      (const std::string &name,
       const SafeSTLVector<T> &vec,
       const Index &dir) const;

  /**
   * @brief Creates a new @ref XMLElement with the given name.
   *
   * @param[in] name Name of the new element created.
   * @return XML element.
   */
  template <class T>
  std::shared_ptr<XMLElement> create_vector_element
      (const std::string &name, const SafeSTLVector<T> &vec) const;

  /**
   * @brief Prints the XML document content.
   *
   * Mostly used for testing and debugging purposes.
   * @param[in] out Log stream for printing information.
   */
  void print_info(LogStream &out) const;

  /**
   * @brief Writes the XML document content to a file.
   *
   * @param[in] file_path File Path of the file.
   */
  void write_to_file(const std::string &file_path) const;

private:

  /** @p Xerces-c DOM document pointer. */
  DOMDocPtr_ xml_doc_;

  /** @p Xerces-c DOM implementation pointer. */
  xercesc::DOMImplementation *dom_impl_;

  /**
   * @brief Checks if the file can be read.
   *
   * It throws an error message if the file can not be read properly.
   *
   * @param[in] file_path Path of the file to be checked.
   */
  static void check_file(const std::string &file_path);

  void initialize_xml();

public:
  template <class T>
  static std::string create_string_from_vector (const SafeSTLVector<T> &vec);

};


IGA_NAMESPACE_CLOSE

#endif // XML_IO

#endif // __XML_DOCUMENT_H_
