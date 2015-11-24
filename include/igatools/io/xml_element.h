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

/**
 * @todo To be documented.
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

  /** Type for a shared pointer of the current class. */
  typedef xercesc::DOMElement * DOMElemPtr_;

  ///@}

  /** @name Constructors, destructor, assignment operators and creators */
  ///@{

  /**
   * @todo To be documented.
   * @param[in] dom_elem Xerces-C XML element object.
   */
  XMLElement(const DOMElemPtr_ dom_elem);

  /**
   * @brief Deleted default constructor.
   *
   * Default constructor.
   * @note Deleted: not allowed.
   */
  XMLElement() = delete;

  /**
   * @brief Deleted copy constructor.
   *
   * Copy constructor.
   * @note Deleted: not allowed.
   */
  XMLElement(const XMLElement &) = delete;

  /**
   * @brief Deleted move constructor.
   *
   * Move constructor.
   * @note Deleted: not allowed.
   */
  XMLElement(XMLElement &&) = delete;

  /**
   * @brief Deleted copy assignment operator.
   *
   * Copy assignment operator.
   * @note Deleted: not allowed.
   */
  XMLElement &operator= (const XMLElement &) = delete;

  /**
   * @brief Deleted move assignment operator.
   *
   * Move assignment operator.
   * @note Deleted: not allowed.
   */
  XMLElement &operator= (XMLElement &&) = delete;

public:
  /**
   * @brief Returns a new instance wrapped into a shared pointer.
   *
   * Builds and returns a new instance of the class wrapped into a
   * shared pointer. It uses the above defined constructor.
   *
   * @param[in] dom_elem Xerces-C XML element object.
   * @return A shared pointer with a new instance of the class.
   */
  static SelfPtr_ create(const DOMElemPtr_ dom_elem);

  ///@}

  /**
   * @brief Returns all the first level children elements contained in
   * the current element.
   *
   * Returns a vector containing all the first level children nodes of
   * type element present in current element.
   *
   * @return Vector containing the children elements.
   */
  SafeSTLVector<SelfPtr_> get_children_elements() const;



  /**
   * @brief Returns all the first level children elements contained in
   * the current element with a given tag @p name.
   *
   * Returns a vector containing all the first level children nodes of
   * type element present in current element with a given tag @p name.
   *
   * @param[in] name Name of the children elements to be extracted.
   * @return Vector containing the children elements.
   */
  SafeSTLVector<SelfPtr_> get_children_elements(const std::string &name) const;


  /**
   * @brief Checks if a certain attribute is the element.
   *
   * Checks if there exists at least one attribute in the element with
   * the given @p name.
   *
   * @param[in] name Name of the attribute to be checked.
   * @return @p true if the attribute is present, @p false elsewhere.
   */
  bool
  has_attribute(const std::string &name) const;


  /**
   * @brief Checks if a certain child element is the element.
   *
   * Checks if there exists at least one child element in the element with
   * the given @p name.
   *
   * @param[in] name Name of the child element to be checked.
   * @return @p true if the child element is present, @p false elsewhere.
   */
  bool
  has_element(const std::string &name) const;

  /**
   * @brief Returns the tag @p name of the element.
   *
   * Returns the tag @p name of the element.
   *
   * @return Name of the node.
   */
  std::string get_name() const;


  /**
   * @brief Returns the value of the element.
   *
   * Return the value of the element given the template type.
   *
   * @return Parsed value returned.
   * @tparam T Type of value returned.
   */
  template <class T> T get_value() const;

  /**
   * @brief Returns the only one text child element contained in the
   * element.
   *
   * Returns the only one text child element of type @p DOMText contained
   * in the element.
   *
   * @return Single extracted element wrapped into a shared pointer.
   *
   * @warning In debug mode, if more than one children of type @p DOMText
   * are present the element, or there is anyone, an error is thrown.
   */
  xercesc::DOMText *get_single_text_element() const;

  /**
   * @brief Returns the value of the attribute with the given @p name.
   *
   * Returns the value of type @p T of an attribute with the given @p name
   * in the element.
   *
   * @param[in] name Name of the attribute.
   * @return Value of the attribute.
   *
   * @tparam T Type of value returned.
   * @note In debug mode, if the attribute is not present,
   * an error is thrown.
   */
  template <class T> T get_attribute(const std::string &name) const;

  /**
   * @brief Returns a vector of numerical values contained in the element.
   *
   * Returns a vector of numerical values with type @p T that
   * are contained in the element.
   *
   * @tparam T Type of value returned in the vector.
   * @return Vector containing the extracted numerical values.
   */
  template <class T>
  SafeSTLVector<T> get_values_vector();

  /**
   * @brief Returns the only one child element contained in the element.
   *
   * Returns the only one child element contained in the element.
   *
   * @return Single extracted element.
   *
   * @note In debug mode, if more than one children element, an error
   * is thrown.
   */
  SelfPtr_ get_single_element();

  /**
   * @brief Returns the only one child element contained in the element
   * with the given tag @p name.
   *
   * Returns the only one child element contained in the element that has
   * the given tag @p name.
   *
   * @param[in] name Name of the element to be extracted.
   * @return Single extracted element.
   *
   * @note In debug mode, if more than one children element that has
   * the given @p name, an error is thrown.
   */
  SelfPtr_ get_single_element(const std::string &name);

private:
  /**
   * @todo To be documented.
   */
  const DOMElemPtr_ root_elem_;

};


IGA_NAMESPACE_CLOSE

#endif // XML_IO

#endif // __XML_ELEMENT_H_

