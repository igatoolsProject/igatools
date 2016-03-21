//-+--------------------------------------------------------------------
// Igatools a general purpose Isogeometric analysis library.
// Copyright (C) 2012-2016  by the igatools authors (see authors.txt).
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

#ifdef IGATOOLS_WITH_XML_IO

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
 * @brief Class for managing XML DOM documents of @p Xerces-c.
 *
 * This is the main class for dealing with XML documents.
 * It is a wrapper class for the XML <tt>Xerces-c DOMDocument</tt> class.
 *
 * It allows to create new XML DOM void documents (by means of the method
 * @ref create_void_document) or to parse documents from input files
 * by calling @ref parse_from_file).
 *
 * Parsed files can be validated against a XSD schema.
 *
 * Once the document is created, the main XML element can be extracted
 * wrapped into a @ref XMLElement by calling @ref get_document_element.
 *
 * XML document can be extended new elements (using methods
 * @ref create_new_element, @ref create_new_text_element and
 * @ref create_vector_element) that will be append to the existing elements.
 *
 * XML documents can be written to a file by calling the method
 * @ref write_fo_file.
 *
 *
 * @note This class uses @p Xerces-c library.
 * @note The class is not in charge of freeing the @p DOMDocument pointer.
 *
 * @see XMLElement
 *
 * @author P. Antolin
 * @date 2015
 */
class XMLDocument
{
private:

  class ParserErrorHandler;

  /** @name Types and static values */
  ///@{

  /** Type for the current class. */
  typedef XMLDocument Self_;

  /** Type for a shared pointer of the current class. */
  typedef std::shared_ptr<Self_> SelfPtr_;

  /** Type for a pointer of the @p Xerces-c DOM document class. */
  typedef xercesc::DOMDocument *DOMDocPtr_;

  /** Type for a shared pointer of @ref XMLElement. */
  typedef std::shared_ptr<XMLElement> XMLElemPtr_;

  /** Default precision for writing vectors of real numbers. */
  static const int default_precision_;

  ///@}

  /** @name Constructors, destructor, assignment operators and creators */
  ///@{

  /**
   * @brief Constructs a new void XML document.
   *
   * It creates a new document with a single XML element with
   * the given @p name.
   * Once the document has been created, it can be filled by creating
   * new XML elements.
   *
   * @param[in] name Name of the XML element.
   */
  XMLDocument(const char *name);

  /**
   * @brief Constructs a new XML document by parsing it from the
   * @p file_path.
   *
   * Before parsing the file the validity of the file is checked by
   * calling the static method @ref check_file.
   *
   * Once the document has been created, the information can be
   * extracted and the document can be populated.
   *
   * @warning If there is any problem parsing the input file, error
   * messages and exceptions will be thrown.
   *
   * @param[in] file_path File to be parsed.
   */
  XMLDocument(const std::string &file_path);

  /**
   * @brief Constructs a new XML document by parsing it from the
   * @p file_path and validating it against @p grammar_definition.
   *
   * Before parsing the file the validity of the file is checked by
   * calling the static method @ref check_file.
   *
   * Once the document has been created, the information can be
   * extracted and the document can be populated.
   *
   * @warning If there is any problem parsing the input file, error
   * messages and exceptions will be thrown.
   *
   * @param[in] file_path File to be parsed.
   * @param[in] grammar_definition String containing the XSD grammar
   * definition to validate the XML document.
   */
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
   * It finalizes all the @p Xerces-c running process and frees
   * all the structures.
   */
  ~XMLDocument();

  /**
   * @brief Creates a new void document with a single XML element with
   * the given @p name.
   *
   * @param[in] name Name of the XML element.
   * @return A shared pointer with a new instance of the class.
   */
  static SelfPtr_ create_void_document(const std::string &name);

  /**
   * @brief Parses the @p file_path and returns a instance of the class
   * wrapped into a shared pointer.
   *
   * @param[in] file_path Path of the file to be parsed.
   * @return XML document containing the parsed document.
   */
  static SelfPtr_ parse_from_file(const std::string &file_path);

  /**
   * @brief Parses the @p file_path, validates it against the
   * @p grammar_definition and returns a instance of the class
   * wrapped into a shared pointer.
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
   * @brief Retrieves the main document element wrapped into a
   * @ref XMLElement.
   *
   * @return DOM document element.
   */
  XMLElemPtr_ get_document_element() const;

  /**
   * @brief Creates a new @ref XMLElement with the given @p name.
   *
   * @param[in] name Name of the new element created.
   * @return XML element.
   */
  XMLElemPtr_ create_new_element(const std::string &name) const;

  /**
   * @brief Creates a new @ref XMLElement with the given @p name and
   * containing the given @p text.
   *
   * @param[in] name Name of the new element created.
   * @param[in] text Text contained in the element.
   * @return XML element.
   */
  XMLElemPtr_ create_new_text_element(const std::string &name,
                                      const std::string &text) const;

  /**
   * @brief Creates a new @ref XMLElement with the given @p name and
   * containing a @p vector of values of the type @p T
   *
   * @tparam T Type of the entries of the vector to be written.
   * @param[in] name Name of the new element created.
   * @param[in] vector Vector of values to be added.
   * @return XML element.
   */
  template <class T>
  XMLElemPtr_ create_vector_element(const std::string &name,
                                    const SafeSTLVector<T> &vector,
                                    const int &precision = default_precision_,
                                    const bool scientific_format = true) const;

  /**
   * @brief Prints the XML document content.
   *
   * Mostly used for testing and debugging purposes.
   * @param[in] out Log stream for printing information.
   */
  void print_info(LogStream &out) const;

  /**
   * @brief Writes the XML document content to the given @p file_path.
   *
   * @param[in] file_path File Path of the file.
   * @param[in] pretty_print Flag indicating if the file must
   * written with a pretty print format.
   */
  void write_to_file(const std::string &file_path,
                     const bool pretty_print = false) const;

  /**
   * @brief Checks if the file can be read.
   *
   * It throws an error message if the file can not be read properly.
   *
   * @param[in] file_path Path of the file to be checked.
   */
  static void check_file(const std::string &file_path);

private:

  /// @p Xerces-c DOM document pointer.
  DOMDocPtr_ xml_doc_;

  /// @p Xerces-c DOM implementation pointer.
  xercesc::DOMImplementation *dom_impl_;

  /**
   * @brief Method for initializing all the @p Xerces-c processes.
   *
   * This is called just once, by the constructor.
   */
  void initialize_xml();

//public:

  /**
   * @brief Writes the passed @p vector as a string with spaces between
   * the entries.
   *
   * @tparam T Type of the entries of the vector.
   * @param[in] vector Vector to be written as string.
   */
  template <class T>
  static std::string create_string_from_vector(const SafeSTLVector<T> &vector,
                                               const int &precision = default_precision_,
                                               const bool scientific_format = true);

private:

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
   * @see XMLDocument
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

};


IGA_NAMESPACE_CLOSE

#endif // IGATOOLS_WITH_XML_IO

#endif // __XML_DOCUMENT_H_
