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

#ifndef XML_FILE_PARSER_H_
#define XML_FILE_PARSER_H_

#if XML_IO

#include <igatools/base/config.h>

#include <xercesc/util/XercesDefs.hpp>
#include <xercesc/sax/ErrorHandler.hpp>

XERCES_CPP_NAMESPACE_BEGIN
class XercesDOMParser;
class DOMDocument;
XERCES_CPP_NAMESPACE_END


IGA_NAMESPACE_OPEN

/**
 * Forward declarations.
 */

class XMLParserErrorHandler;

/**
 * @brief Class for parsing input files.
 *
 * This is a class for parsing XML input files and validate them againts
 * a XML Schema grammar.
 *
 * This class provides the capability of creating a @p XercesDOMParser,
 * checking the validity of the input file (if it exists, it is corrupted, etc),
 * and retrieve and XML @p DOMDocument where the input file the information is
 * contained.
 *
 * This class uses a @ref XMLParserErrorHandler for managing the possible errors
 * than can appear during the parsing process.
 *
 * The destructor of the class is in charge of of deleting @ref parser_ and
 * shutdown all the @p xerces active process.
 *
 * \author antolin, 2015
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

private:

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
   * @param[in] eh Error handler for managing parsing errors.
   * @param[in] file_path Path of the file to be parsed.
   */
  XMLFileParser(const std::shared_ptr<XMLParserErrorHandler> eh,
             const std::string &file_path);

  /**
   * @brief Deleted default constructor.
   *
   * Default constructor.
   * @note Deleted: not allowed.
   */
  XMLFileParser() = delete;

  /**
   * @brief Deleted copy constructor.
   *
   * Copy constructor.
   * @note Deleted: not allowed.
   */
  XMLFileParser(const XMLFileParser &) = delete;

  /**
   * @brief Deleted move constructor.
   *
   * Move constructor.
   * @note Deleted: not allowed.
   */
  XMLFileParser(XMLFileParser &&) = delete;

  /**
   * @brief Deleted copy assignment operator.
   *
   * Copy assignment operator.
   * @note Deleted: not allowed.
   */
  XMLFileParser &operator= (const XMLFileParser &) = delete;

  /**
   * @brief Deleted move assignment operator.
   *
   * Move assignment operator.
   * @note Deleted: not allowed.
   */
  XMLFileParser &operator= (XMLFileParser &&) = delete;


public:

  /**
   * @brief Destructor.
   *
   * Destructor: finalizes the parsing process: deletes the @ref parser_ and terminates
   * the platform utils.
   */
  ~XMLFileParser();

  /**
   * @brief Returns a new instance wrapped into a shared pointer.
   *
   * Builds and returns a new instance of the class wrapped into a shared pointer.
   * It uses the above defined default constructor.
   *
   * @param[in] eh Error handler for managing parsing errors.
   * @param[in] file_path Path of the file to be parsed.
   * @return A shared pointer with a new instance of the class.
   */
  static SelfPtr_ create(const std::shared_ptr<XMLParserErrorHandler> eh,
                         const std::string &file_path);

  ///@}


private:

  /// Path of the file to be parsed.
  /** Path of the file to be parsed. */
  const std::string file_path_;

  /// DOM parser.
  /** DOM parser. */
  xercesc::XercesDOMParser *parser_;

  /// Error handler for managing parsing errors.
  /** Error handler for managing parsing errors. */
  const std::shared_ptr<XMLParserErrorHandler> eh_;


public:
  /**
   * @brief Parses the file and returns a XML document object.
   *
   * Parses the file and returns a XML document object. Before parsing the file,
   * all the required configuration flags are set and the needed Schema
   * grammars are loaded.
   *
   * @attention If there is any problem parsing the input file, error messages
   * and exceptions will be thrown.
   *
   * @param[in] load_grammar Pointer to a function in charge of loading the
   *                         XML Schema grammar into the @ref parser_.
   *
   * @return XML document object.
   */
  xercesc::DOMDocument *
  parse(void (*load_grammar)(xercesc::XercesDOMParser *const));

private:
  /**
   * @brief Checks if the @ref file_path_ can be read
   *
   * Checks if the @ref file_path_ can be read. It throws an error message
   * if the file can not be read properly.
   */
  void check_file() const;

};

IGA_NAMESPACE_CLOSE

#endif // XML_IO

#endif // XML_FILE_PARSER_H_

