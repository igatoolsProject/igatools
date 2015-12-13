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
#include <xercesc/sax/ErrorHandler.hpp>

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
 * @brief Manage runtime errors occurred during the parsing process of
 * input files.
 *
 * This class is in charge of managing the runtime errors than can
 * occurred during the parsing process of input files.
 *
 * The class derives from the @p Xerces-c class @ref ParserErrorHandler,
 * in such a way that it takes care of error messages thrown by
 * @p Xerces-c parser.
 *
 * @warning When an error or warning is thrown, an exception raises,
 * and, if not caught, the execution is automatically finished.
 *
 * @see XMLFileParser
 *
 * @author P. Antolin
 * @date 2015
 */
class XMLParserErrorHandler : public xercesc::ErrorHandler
{
private:

  /** @name Types and static values */
  ///@{

  /** Type for the current class. */
  typedef XMLParserErrorHandler Self_;

  /** Type for a shared pointer of the current class. */
  typedef std::shared_ptr<Self_> SelfPtr_;

  ///@}

  /** @name Constructors, destructor, assignment operators and creators */
  ///@{

  /**
   * @brief Default constructor.
   * @return New instance of the class.
   */
  XMLParserErrorHandler() = default;

  /**
   * @brief Copy constructor.
   * @note Deleted, not allowe to be used.
   */
  XMLParserErrorHandler(const XMLParserErrorHandler &) = delete;

  /**
   * @brief Move constructor.
   * @note Deleted, not allowed to be used.
  */
  XMLParserErrorHandler(XMLParserErrorHandler &&) = delete;

  /**
   * @brief Copy assignment operator.
   * @note Deleted, not allowed to be used.
   */
  XMLParserErrorHandler &operator= (const XMLParserErrorHandler &) = delete;

  /**
   * @brief Move assignment operator.
   * @note Deleted, not allowed to be used.
   */
  XMLParserErrorHandler &operator= (XMLParserErrorHandler &&) = delete;

public:
  /**
   * @brief Returns a new instance wrapped into a shared pointer.
   *
   * It uses the above defined default constructor.
   *
   * @return A shared pointer with a new instance of the class.
   */
  static SelfPtr_ create();

  ///@}

  /** @name xercesc ParserErrorHandler pure virtual methods */
  ///@{

  /**
   * @brief Receive notification of a recoverable error.
   * This corresponds to the definition of "error" in section 1.2 of the
   * W3C XML 1.0 Recommendation. For example, a validating parser would
   * use this callback to report the violation of a validity constraint.
   * The default behaviour is to take no action.
   *
   * The SAX parser must continue to provide normal parsing events after
   * invoking this method: it should still be possible for the application
   * to process the document through to the end. If the application cannot
   * do so, then the parser should report a fatal error even if the XML
   * 1.0 recommendation does not require it to do so.
   *
   * An error message is thrown containing additional information about
   * the line and the column where it occurred.
   *
   * @param[in] ex The warning information encapsulated in a SAX parse
   *               exception.
   *
   * @exception SAXException Any SAX exception, possibly wrapping another
   *                         exception.
   *
   * @note The documentation has been partially extracted from the parent
   * <tt>Xerces-c ParserErrorHandler</tt> class.
   */
  virtual void error(const xercesc::SAXParseException &ex) override final;

  /**
   * @brief Receive notification of a warning.
   *
   * SAX parsers will use this method to report conditions that are not
   * errors or fatal errors as defined by the XML 1.0 recommendation.
   * The default behaviour is to take no action.
   * The SAX parser must continue to provide normal parsing events after
   * invoking this method: it should still be possible for the application
   * to process the document through to the end.
   *
   * \attention Currently the warning message is treated as an error
   * message, i.e. the execution will be finished after throwing the
   * message.
   *
   * An error message is thrown containing additional information about
   * the line and the column where it occurred.
   *
   * @param[in] ex The warning information encapsulated in a SAX parse
   *               exception.
   *
   * @exception SAXException Any SAX exception, possibly wrapping another
   *                         exception.
   *
   * @note The documentation has been partially extracted from the parent
   * <tt>Xerces-c ParserErrorHandler</tt> class.
   */
  virtual void warning(const xercesc::SAXParseException &ex) override final;

  /**
   * @brief Receive notification of a non-recoverable error.
   *
   * This corresponds to the definition of "fatal error" in section 1.2 of
   * the W3C XML 1.0 Recommendation. For example, a parser would use this
   * callback to report the violation of a well-formedness constraint.
   *
   * The application must assume that the document is unusable after the
   * parser has invoked this method, and should continue (if at all) only
   * for the sake of collecting addition error messages: in fact, SAX
   * parsers are free to stop reporting any other events once this method
   * has been invoked.
   *
   * An error message is thrown containing additional information about
   * the line and the column where it occurred.
   *
   * @param[in] ex The error information encapsulated in a SAX parse
   *               exception.
   *
   * @exception SAXException  Any SAX exception, possibly wrapping another
   *                          exception.
   *
   * @note The documentation has been partially extracted from the parent
   * <tt>Xerces-c ParserErrorHandler</tt> class.
   * @note Currently, this function is doing nothing.
   */
  virtual void fatalError(const xercesc::SAXParseException &ex) override final;

  /**
   * @brief Reset the Error handler object on its reuse.
   *
   * This method helps in reseting the Error handler object implementation
   * defaults each time the Error handler is begun.
   *
   * @note The documentation has been extracted from the parent
   * <tt>xercesc ParserErrorHandler</tt> class.
   *
   * @note Currently, this function is doing nothing.
   */
  virtual void resetErrors() override final;

  ///@}

};



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
