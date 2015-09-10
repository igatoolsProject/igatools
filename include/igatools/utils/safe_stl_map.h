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

#ifndef __SAFE_STL_MAP_H_
#define __SAFE_STL_MAP_H_

#include <igatools/base/config.h>
#include <igatools/base/print_info_utils.h>

#include <map>

IGA_NAMESPACE_OPEN

/**
 *
 * @ingroup serializable
 */
template <class Key, class T>
class SafeSTLMap : public std::map<Key, T>
{
public:
  /** Inherit the constructors of the base class. */
  using std::map<Key, T>::map;

  /**
   * Returns the number of entries in the container.
   */
  Size size() const noexcept
  {
    return std::map<Key, T>::size();
  }

  /**
   * @name Printing info
   */
  ///@{
private:
  template <class A>
  EnableIf<has_print_info<A>(0), void>
  t_print_info(LogStream &out) const
  {
    for (auto &entry : *this)
    {
      out << "{";
      entry.first.print_info(out);
      out << ", ";
      entry.second.print_info(out);
      out << "}" << std::endl;
    }
  }

  template <class A>
  EnableIf<(!has_print_info<A>(0)), void>
  t_print_info(LogStream &out) const
  {
    for (auto &entry : *this)
    {
      out << "{ " << entry.first;
      out << ", " << entry.second  << "}";
      out << std::endl;
    }
  }

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
    ar &make_nvp("SafeSTLContainer_Map",
                 base_class<std::map<Key,T>>(this));

  }
  ///@}
#endif // SERIALIZATION

public:

  /**
   * Prints the content of the vector on the LogStream @p out.
   * Its use is intended mainly for testing and debugging purpose.
   */
  void print_info(LogStream &out) const
  {
    t_print_info<T>(out);
  }
  ///@}
};


IGA_NAMESPACE_CLOSE


#endif
