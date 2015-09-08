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

#ifndef __SAFE_STL_VECTOR_H_
#define __SAFE_STL_VECTOR_H_

#include <igatools/base/config.h>
#include <igatools/utils/safe_stl_container.h>
#include <vector>

IGA_NAMESPACE_OPEN

/**
 * @brief iga version of std::vector.
 * It can be used as a std::vector but in Debug mode
 * it provides bounds checking.
 */
template<class T>
class SafeSTLVector :
  public SafeSTLContainer<std::vector<T>>
{
  using base_t = SafeSTLContainer<std::vector<T>>;
public :
  /** Inherit the constructors of the base class. */
  using SafeSTLContainer<std::vector<T>>::SafeSTLContainer;


private:
#ifdef SERIALIZATION
  /**
   * @name Functions needed for serialization
   * @see <a href="http://uscilab.github.io/cereal/serialization_functions.html">Cereal serialization</a>
   */
  ///@{
  friend class cereal::access;

  template<class Archive>
  void serialize(Archive &ar)
  {
    ar &make_nvp("SafeSTLContainer_Vector",
                 base_class<std::vector<T>>(this));
  }
  ///@}
#endif // SERIALIZATION

};

IGA_NAMESPACE_CLOSE

#endif // SAFE_STL_VECTOR_H_
