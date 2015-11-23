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

#ifndef __OBJECTS_CONTAINER_H_
#define __OBJECTS_CONTAINER_H_

#include <igatools/base/config.h>

IGA_NAMESPACE_OPEN

/**
 * @brief To be documented.
 *
 * @author P. Antolin 2015
 *
 * @todo document more
 */
class ObjectsContainer
{
  /** Type for current class. */
  using self_t = ObjectsContainer;

  /** @name Constructors*/
  ///@{
public:
  /**
   * To document
   */
  ObjectsContainer();


public:
  /**
   * Copy constructor.
   */
  ObjectsContainer(const self_t &container) = default;

  /**  Move constructor */
  ObjectsContainer(self_t &&container) = default;

  /** Destructor */
  ~ObjectsContainer() = default;
  ///@}

  /**
   * @name Creators
   * @note The functions here return:
   * - a <b>non-const</b> Container object wrapped by a std::shared_ptr
   * in the case of <b>create()</b> functions;
   * - a <b>const</b> Container object wrapped by a std::shared_ptr
   * in the case of <b>const_create()</b> functions.
   */
  ///@{
  /**
   * Creates a objects container (non-const).
   */
  static std::shared_ptr<self_t> create();

  /**
   * Creates a objects container (const).
   */
  static std::shared_ptr<const self_t> const_create();
  ///@}

  /**
   * @name Assignment operators
   */
  ///@{

  /**
   * Copy assignment operator.
   */
  self_t &operator=(const self_t &container) = default;

  /**
   * Move assignment operator.
   */
  self_t &operator=(self_t &&container) = default;
  ///@}


};

IGA_NAMESPACE_CLOSE

#endif /*OBJECTS_CONTAINER_H_ */
