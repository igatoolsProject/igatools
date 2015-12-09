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

#ifndef __XML_ELEMENT_H_
#define __XML_ELEMENT_H_

#include <igatools/base/config.h>

#ifdef XML_IO

#include <xercesc/util/XercesDefs.hpp>
#include <xercesc/dom/DOMElement.hpp>
XERCES_CPP_NAMESPACE_BEGIN
class DOMElement;
class DOMText;
XERCES_CPP_NAMESPACE_END


IGA_NAMESPACE_OPEN

template <class T> class SafeSTLVector;
class LogStream;

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
class XMLElement
{
private:

  /** @name Types and static values */
  ///@{

  /** Type for the current class. */
  typedef XMLElement Self_;

  /** Type for a shared pointer of the current class. */
  typedef std::shared_ptr<Self_> SelfPtr_;

  /** Type for a pointer of the @p Xerces-c DOM element class. */
  typedef xercesc::DOMElement *DOMElemPtr_;

  ///@}

  /** @name Constructors, destructor, assignment operators and creators */
  ///@{

  /**
   * Constructor.
   *
   * Constructs a new instance of the class by taking a pointer
   * to a <tt>Xerces-c DOMElement</tt>.
   *
   * @param[in] dom_elem @p Xerces-c XML element object.
   */
  XMLElement(const DOMElemPtr_ dom_elem);

  /**
   * @brief Default constructor.
   * @note Deleted, not allowed to be used.
   */
  XMLElement() = delete;

  /**
   * @brief Copy constructor.
   * @note Deleted, not allowed to be used.
   */
  XMLElement(const XMLElement &) = delete;

  /**
   * @brief Move constructor.
   * @note Deleted, not allowed to be used.
   */
  XMLElement(XMLElement &&) = delete;

  /**
   * @brief Copy assignment operator.
   * @note Deleted, not allowed to be used.
   */
  XMLElement &operator= (const XMLElement &) = delete;

  /**
   * @brief Move assignment operator.
   * @note Deleted, not allowed to be used.
   */
  XMLElement &operator= (XMLElement &&) = delete;

public:
  /**
   * @brief Builds and returns a new instance wrapped into a shared pointer.
   *
   * It uses the above defined constructor.
   *
   * @param[in] dom_elem @p Xerces-c XML element object.
   * @return A shared pointer with a new instance of the class.
   */
  static SelfPtr_ create(const DOMElemPtr_ dom_elem);

  ///@}

  /**
   * @brief Returns a vector containing all the first level children nodes of
   * type element present in current element.
   *
   * @return Vector containing the children elements.
   */
  SafeSTLVector<SelfPtr_> get_children_elements() const;

  /**
   * @brief Returns a vector containing all the first level children nodes of
   * type element present in current element with a given tag @p name.
   *
   * @param[in] name Name of the children elements to be extracted.
   * @return Vector containing the children elements.
   */
  SafeSTLVector<SelfPtr_> get_children_elements(const std::string &name) const;

  /**
   * @brief Checks if there exists at least one attribute in the element
   * with the given @p name.
   *
   * @param[in] name Name of the attribute to be checked.
   * @return @p true if the attribute is present, @p false elsewhere.
   */
  bool has_attribute(const std::string &name) const;

  /**
   * @brief Checks if there exists at least one child element in the
   * element with the given @p name.
   *
   * @param[in] name Name of the child element to be checked.
   * @return @p true if the child element is present, @p false elsewhere.
   */
  bool has_element(const std::string &name) const;

  /**
   * @brief Returns the tag @p name of the element.
   * @return Name of the node.
   */
  std::string get_name() const;

  /**
   * @brief Return the value of type @p T of the element.
   *
   * @tparam T Type of value returned.
   * @return Parsed value returned.
   *
   * @warning It throws an exception if the type of the value
   * is not the specified one.
   */
  template <class T> T get_value() const;

  /**
   * @brief Returns the value of type @p T of an attribute with the given
   * @p name in the element.
   *
   * @tparam T Type of value returned.
   * @param[in] name Name of the attribute.
   * @return Value of the attribute.
   *
   * @warning In debug mode, if the attribute is not present,
   * an error is thrown.
   */
  template <class T> T get_attribute(const std::string &name) const;

  /**
   * @brief Returns a vector of values with type @p T that are contained
   * in the element.
   *
   * @tparam T Type of value returned in the vector.
   * @return Vector containing the extracted numerical values.
   *
   * @warning It will throw an exception if is not able to parse the vector.
   */
  template <class T>
  SafeSTLVector<T> get_values_vector() const;

  /**
   * @brief Returns the only one child element contained in the element.
   *
   * @return Single extracted element.
   *
   * @warning In debug mode, if more than one children element, an error
   * is thrown.
   */
  SelfPtr_ get_single_element();

  /**
   * @brief Returns the only one child element contained in the element
   * with the given tag @p name.
   *
   * @param[in] name Name of the element to be extracted.
   * @return Single extracted element.
   *
   * @warning In debug mode, if more than one children element that has
   * the given @p name, an error is thrown.
   */
  SelfPtr_ get_single_element(const std::string &name);

  /**
   * @brief Prints the XML element content.
   *
   * Mostly used for testing and debugging purposes.
   * @param[in] out Log stream for printing information.
   */
  void print_info(LogStream &out) const;

private:

  /** @p Xerces-c DOM element pointer. */
  const DOMElemPtr_ root_elem_;

  /**
   * @brief Returns the only one text child element contained in the
   * element.
   *
   * @return Single extracted element pointer.
   *
   * @warning In debug mode, if more than one children of type @p DOMText
   * are present the element, or there is none, an error is thrown.
   */
  xercesc::DOMText *get_single_text_element() const;

};


IGA_NAMESPACE_CLOSE

#endif // XML_IO

#endif // __XML_ELEMENT_H_
