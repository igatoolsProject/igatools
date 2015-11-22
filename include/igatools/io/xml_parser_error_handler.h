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

#ifndef XML_PARSER_ERROR_HANDLER_H_
#define XML_PARSER_ERROR_HANDLER_H_

#include <igatools/base/config.h>

#ifdef XML_IO

#include <xercesc/sax/ErrorHandler.hpp>

IGA_NAMESPACE_OPEN


/**
 * @brief Manage runtime errors occurred during the parsing process of
 * input files.
 *
 * This class is in charge of managing the runtime errors than can
 * occurred during the parsing process of input files.
 *
 * The class derives from the @p xercesc class @ref ParserErrorHandler, in such
 * a way, it takes care of error messages thrown by @p xercesc parser.
 *
 * @attention When an error or warning is thrown, an exception raises,
 * and the execution is automatically finished.
 *
 * @author P. Antolin
 * @date 2015
 */
class XMLParserErrorHandler : public xercesc::ErrorHandler
{
private:

  /** @name Types and static values */
  ///@{

  /// Type for the current class.
  /** Type for the current class. */
  typedef XMLParserErrorHandler Self_;

  /// Type for a shared pointer of the current class.
  /** Type for a shared pointer of the current class. */
  typedef std::shared_ptr<Self_> SelfPtr_;

  ///@}

private:
  /** @name Constructors, destructor, assignment operators and creators */
  ///@{

  /**
   * @brief Default constructor
   *
   * Default constructor.
   * @return New instance of the class.
   */
  XMLParserErrorHandler() = default;

  /**
   * @brief Deleted copy constructor.
   *
   * Copy constructor.
   * @note Deleted: not allowed.
   */
  XMLParserErrorHandler(const XMLParserErrorHandler &) = delete;

  /**
   * @brief Deleted move constructor.
   *
  * Move constructor.
  * @note Deleted: not allowed.
  */
  XMLParserErrorHandler(XMLParserErrorHandler &&) = delete;

  /**
   * @brief Deleted copy assignment operator.
   *
   * Copy assignment operator.
   * @note Deleted: not allowed.
   */
  XMLParserErrorHandler &operator= (const XMLParserErrorHandler &) = delete;

  /**
   * @brief Deleted move assignment operator.
   *
   * Move assignment operator.
   * @note Deleted: not allowed.
   */
  XMLParserErrorHandler &operator= (XMLParserErrorHandler &&) = delete;

public:
  /**
   * @brief Returns a new instance wrapped into a shared pointer.
   *
   * Builds and returns a new instance of the class wrapped into a
   * shared pointer. It uses the above defined default constructor.
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
   * <tt>xercesc ParserErrorHandler</tt> class.
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
   * <tt>xercesc ParserErrorHandler</tt> class.
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
   * <tt>xercesc ParserErrorHandler</tt> class.
   * @note This function is doing nothing currently.
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
   * @note The documentation has been extracted from the parent
   * <tt>xercesc ParserErrorHandler</tt> class.
   *
   * @note Currently, this function is doing nothing.
   */
  virtual void resetErrors() override final;

  ///@}

};


IGA_NAMESPACE_CLOSE

#endif // XML_IO

#endif // XML_PARSER_ERROR_HANDLER_H_
