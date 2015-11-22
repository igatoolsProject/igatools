//-+--------------------------------------------------------------------
// Isolde (Isogeometric Solver Environment) is a software for
// analyzing continuum mechanics problems by means of Isogeometric Methods.
// Copyright (C) 2012-2014 by the isolde authors (see authors.txt).
//
// This file is part of the isolde software.
//
// Isolde is property of the University of Pavia and IMATI / CNR,
// Italy. It can not be neither redistributed nor modify without
// an express authorization of the authors.
//
// This software is based on igatools library, that is distributed
// under GNU General Public License version 3.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
//-+--------------------------------------------------------------------

#ifndef XML_PARSER_ERROR_HANDLER_H_
#define XML_PARSER_ERROR_HANDLER_H_

#if XML_IO

#include <igatools/base/config.h>
#include <igatools/utils/safe_stl_vector.h>

#include <xercesc/sax/ErrorHandler.hpp>


IGA_NAMESPACE_OPEN


/**
 * Forward declarations.
 */
//class InputTreeLog;
class LogManager;


/**
 * @brief Manage runtime errors occurred during the parser of the input files.
 *
 * This class is in charge of managing the runtime error than can occurred
 * during the parser of the input files.
 *
 * The class derives from the @p xercesc class @ref ErrorHandler, in such a way,
 * it takes care of error messages thrown by @p xercesc parser.
 *
 * It also allows to throw error occurred during the parsing process by means
 * of the functions @ref throw_error. This function makes use of the
 * @ref InputTreeLog object contained inside to produce richer error messages.
 *
 * This class can track different nested input files, with their respective
 * @ref InputTreeLog
 *
 * @attention When an error or warning is thrown, an exception raises,
 * and the execution is automatically finished.
 *
 * \author antolin, 2015
 */
class XMLFileParserHandler : public xercesc::ErrorHandler
{
private:

  /** @name Types and static values */
  ///@{

  /// Type for the current class.
  /** Type for the current class. */
  typedef XMLFileParserHandler Self_;

  /// Type for a shared pointer of the current class.
  /** Type for a shared pointer of the current class. */
  typedef std::shared_ptr<Self_> SelfPtr_;

  ///@}


  /** @name Constructors, destructor, assignment operators and creators */
  ///@{

  /**
   * @brief Constructor taking a input tree log and the file path.
   *
   * Constructor taking a input tree log and the file path.
   *
   * @param[in] log Input tree log to be used for logging purposes.
   * @param[in] file_path Path of the file being parsed.
   */
//  XMLFileParserHandler(const std::shared_ptr<InputTreeLog> log,
//                     const std::shared_ptr<LogManager> log_manager,
//                     const std::string &file_path);
  XMLFileParserHandler(const std::shared_ptr<LogManager> log_manager,
                     const std::string &file_path);

public:
  /**
   * @brief Returns a new instance wrapped into a @p shared_ptr
   *
   * Builds and returns a new instance of the class wrapped into a shared pointer.
   * It uses the above defined constructor for @p log and @p file_path.
   *
   * @param[in] log Input tree log to be used for logging purposes.
   * @param[in] file_path Path of the file being parsed.
   */
//  static SelfPtr_ create(const std::shared_ptr<InputTreeLog> log,
//                         const std::shared_ptr<LogManager> log_manager,
//                         const std::string &file_path);
  static SelfPtr_ create(const std::shared_ptr<LogManager> log_manager,
                         const std::string &file_path);


private:
  /**
   * @brief Deleted default constructor.
   *
   * Default constructor.
   * @note Deleted: not allowed.
   */
  XMLFileParserHandler() = delete;

  /**
   * @brief Deleted copy constructor.
   *
   * Copy constructor.
   * @note Deleted: not allowed.
   */
  XMLFileParserHandler(const XMLFileParserHandler &) = delete;

  /**
   * @brief Deleted move constructor.
   *
  * Move constructor.
  * @note Deleted: not allowed.
  */
  XMLFileParserHandler(XMLFileParserHandler &&) = delete;

  /**
   * @brief Deleted copy assignment operator.
   *
   * Copy assignment operator.
   * @note Deleted: not allowed.
   */
  XMLFileParserHandler &operator= (const XMLFileParserHandler &) = delete;

  /**
   * @brief Deleted move assignment operator.
   *
   * Move assignment operator.
   * @note Deleted: not allowed.
   */
  XMLFileParserHandler &operator= (XMLFileParserHandler &&) = delete;

  ///@}

private:

  /// Input tree logs used for throwing error messages.
  /**
   * Vector of input tree logs used for throwing error messages.
   * Every log corresponds to an input file.
   */
//  SafeSTLVector<std::shared_ptr<InputTreeLog>> logs_;

  /// Log manager used for outputting messages.
  /** Log manager used for outputting messages. */
  const std::shared_ptr<LogManager> log_manager_;

  /// Paths of the input files being parsed.
  /**
   *  Paths of the input files being parsed.
   *
   *  They are store in order, first vector entry corresponds to the main
   *  input file. Subsequent elements, to the nested input files.
   */
  SafeSTLVector<std::string> file_paths_;


  /** @name xercesc ErrorHandler pure virtual methods */
  ///@{

public:

  /**
   * @brief Receive notification of a recoverable error.
   * This corresponds to the definition of "error" in section 1.2 of the W3C
   * XML 1.0 Recommendation. For example, a validating parser would use this
   * callback to report the violation of a validity constraint. The default
   * behaviour is to take no action.
   *
   * The SAX parser must continue to provide normal parsing events after
   * invoking this method: it should still be possible for the application
   * to process the document through to the end. If the application cannot
   * do so, then the parser should report a fatal error even if the XML 1.0
   * recommendation does not require it to do so.
   *
   * The error message is thrown by means of the method @ref error.
   * It contains information about the input file, line and column, and
   * the XML tree location.
   *
   * @param[in] ex The warning information encapsulated in a SAX parse
   *               exception.
   *
   * @exception SAXException Any SAX exception, possibly wrapping another
   *                         exception.
   *
   * @note The documentation has been partially extracted from the parent
   * <tt>xercesc ErrorHandler</tt> class.
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
   * \attention Currently the warning message is treated as an error message,
   * i.e. the execution will be finished after throwing the message.
   *
   * The error message is thrown by means of the method @ref error.
   * It contains information about the input file, line and column, and
   * the XML tree location.
   *
   * @param[in] ex The warning information encapsulated in a SAX parse
   *               exception.
   *
   * @exception SAXException Any SAX exception, possibly wrapping another
   *                         exception.
   *
   * @note The documentation has been partially extracted from the parent
   * <tt>xercesc ErrorHandler</tt> class.
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
   * for the sake of collecting addition error messages: in fact, SAX parsers are free to stop reporting any other events once this method has been invoked.
   *
   * The error message is thrown by means of the method @ref error.
   * It contains information about the input file, line and column, and
   * the XML tree location.
   *
   * @param[in] ex The error information encapsulated in a SAX parse exception.
   *
   * @exception SAXException  Any SAX exception, possibly wrapping another
   *                          exception.
   *
   * @note The documentation has been partially extracted from the parent
   * <tt>xercesc ErrorHandler</tt> class.
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
   * <tt>xercesc ErrorHandler</tt> class.
   *
   * @note The documentation has been extracted from the parent
   * <tt>xercesc ErrorHandler</tt> class.
   * @note This function is doing nothing currently.
   */
  virtual void resetErrors() override final;

  ///@}


  /**
   * @brief Throws an exception with an error message if the @p condition
   * is not fulfilled.
   *
   * Throws an exception error if the @p condition is not fulfilled.
   * The error message is generated by included the @p file_path_,
   * the @p line_number and @p column_number, the XML tree location provided
   * by the last of the @p logs_ and the provided @p message.
   *
   * @attention After the error message, the execution is finished.
   *
   * @param[in] condition If it is not fulfilled, the error is thrown.
   * @param[in] message Information for the error message.
   * @param[in] line_number Line number in the input file where the error
   *            occurs. The -1 default value is provided, if this value
   *            is not overridden, the line number is set as @p unknown
   *            in the error message.
   * @param[in] column_number Column number in the input file where the error
   *            occurs. The -1 default value is provided, if this value
   *            is not overridden, the column number is set as @p unknown
   *            in the error message.
   */
  void error(const bool condition,
             const std::string &message,
             const int &line_number = -1,
             const int &column_number = -1) const;

  /**
   * @brief Throws an exception with an error message.
   *
   * Throws an exception error.
   * The error message is generated by included the @p file_path_,
   * the XML tree location provided by the last of the @p logs_ and the
   * provided @p message.
   *
   * @attention After the error message, the execution is finished.
   *
   * @param[in] message Information for the error message.
   */
  void error(const std::string &msg) const;

  /**
   * @brief Throws an exception with an error message.
   *
   * Throws an exception error.
   * The error message is generated by included the @p file_path_,
   * the @p line_number and @p column_number, the XML tree location provided
   * by the last of the @p logs_ and the provided @p message.
   *
   * @attention After the error message, the execution is finished.
   *
   * @param[in] message Information for the error message.
   * @param[in] line_number Line number in the input file where the error
   *            occurs. The -1 default value is provided, if this value
   *            is not overridden, the line number is set as @p unknown
   *            in the error message.
   * @param[in] column_number Column number in the input file where the error
   *            occurs. The -1 default value is provided, if this value
   *            is not overridden, the column number is set as @p unknown
   *            in the error message.
   */
  void error(const std::string &msg,
             const int &line_number,
             const int &column_number) const;


  /** @name Functions for adding/removing tags from the Input Tree Logs */
  ///@{
  /**
   * @brief Adds an @p id to the last tag name defined.
   *
   * Adds an @p id number to the last tag name defined in the last of the
   * @ref logs_.
   *
   * @param[in] id Number to be added to the last tag name.
   *
   * @note It calls the function @p add_id of the @ref @InputTreeLog.
   */
  void add_id(const Index &id);

  /**
   * @brief Appends an @p tag_name at the end of the tree.
   *
   * Appends a new @p tag_name at the end of the tag names tree in the last
   * of the @ref logs_.
   *
   * @param[in] tag_name Name of the tag to be included at the beginning of
   *                     the tree.
   *
   * @note It calls the function @p add_tag_name of the @ref @InputTreeLog.
   */
  void add_tag_name(const std::string &tag_name);

  /**
   * @brief Removes the last tag name of the tree.
   *
   * Removes the last tag name appended to the tree in the last
   * of the @ref logs_.
   *
   * @note It calls the function @p pop_tag_name of the @ref @InputTreeLog.
   */
  void pop_tag_name();
  ///@}

  /**
   * @brief Appends a new nested input file with its corresponding XML tag.
   *
   * Appends a new nested input file with its corresponding XML tag for the
   * main XML element of the file.
   *
   * @param[in] file_path Path of the file being parsed.
   * @param[in] tag_name Name of the tag to be included at the beginning of
   *                     the tree.
   *
   * @note It creates a new @p log and appends it to @ref logs_.
   */
  void add_file(const std::string &file_path, const std::string &tag_name);

  /**
   * @brief Removes the last input file and tree log.
   *
   * Removes the last input file of @ref file_paths_ and the last log of
   * @ref logs_.
   */
  void pop_file();

  /**
   * @brief Prints information of the class to the log stream.
   *
   * Prints information of the class to the log stream.
   *
   * @param[in] out Log in which the information is streamed.
   */
  void print_info(iga::LogStream &out) const;

};


IGA_NAMESPACE_CLOSE

#endif // XML_IO

#endif // XML_IGATOOLS_OBJECTS_PARSER_H_
