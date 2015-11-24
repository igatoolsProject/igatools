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
XERCES_CPP_NAMESPACE_BEGIN
class DOMDocument;
XERCES_CPP_NAMESPACE_END


IGA_NAMESPACE_OPEN

/**
 * @todo To be documented.
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

  ///@}

  /** @name Constructors, destructor, assignment operators and creators */
  ///@{

  /**
   * @todo To be documented.
   */
  XMLDocument(const std::shared_ptr<xercesc::DOMDocument> dom_doc);

  /**
   * @brief Deleted default constructor.
   *
   * Default constructor.
   * @note Deleted: not allowed.
   */
  XMLDocument() = delete;

  /**
   * @brief Deleted copy constructor.
   *
   * Copy constructor.
   * @note Deleted: not allowed.
   */
  XMLDocument(const XMLDocument &) = delete;

  /**
   * @brief Deleted move constructor.
   *
   * Move constructor.
   * @note Deleted: not allowed.
   */
  XMLDocument(XMLDocument &&) = delete;

  /**
   * @brief Deleted copy assignment operator.
   *
   * Copy assignment operator.
   * @note Deleted: not allowed.
   */
  XMLDocument &operator= (const XMLDocument &) = delete;

  /**
   * @brief Deleted move assignment operator.
   *
   * Move assignment operator.
   * @note Deleted: not allowed.
   */
  XMLDocument &operator= (XMLDocument &&) = delete;

public:
  /**
   * @brief Returns a new instance wrapped into a shared pointer.
   *
   * Builds and returns a new instance of the class wrapped into a
   * shared pointer. It uses the above defined constructor.
   *
   * @param[in] dom_doc Xerces-C XML document object.
   * @return A shared pointer with a new instance of the class.
   */
  static SelfPtr_ create(const std::shared_ptr<xercesc::DOMDocument> dom_doc);

  ///@}

};

IGA_NAMESPACE_CLOSE

#endif // XML_IO

#endif // __XML_DOCUMENT_H_
