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

#ifndef SAFE_STL_ARRAY_H_
#define SAFE_STL_ARRAY_H_

#include <igatools/base/config.h>
#include <igatools/utils/safe_stl_container.h>

#include <array>

IGA_NAMESPACE_OPEN



/**
 * @brief iga version of std::array.
 * It can be used as a std::array but in Debug mode
 * it provides bounds checking.
 */
template<class T,int N>
class SafeSTLArray :
  public SafeSTLContainer<std::array<T,N>>
{
  using base_t = SafeSTLContainer<std::array<T,N>>;
public :

  /** Inherit the constructors of the base class. */
  using SafeSTLContainer<std::array<T,N>>::SafeSTLContainer;

  /**
   * Default Constructor. It fills the <tt>N</tt> entries of the array with <tt>N</tt> copies
   * of the default-construced object of type <tt>T</tt>.
   */
  SafeSTLArray() {};


  /**
   * Constructor. It fills the <tt>N</tt> entries of the array with <tt>N</tt> copies
   * of the input argument <tt>val</tt>.
   */
  SafeSTLArray(const T &val)
  {
    std::array<T,N>::fill(val);
  };


  SafeSTLArray(std::initializer_list<T> list)
  {
    Assert(list.size() == N, ExcDimensionMismatch(list.size(),N));

    for (int i = 0 ; i < N ; ++i)
      (*this)[i] = list.begin()[i] ;
  };

private:

#ifdef SERIALIZATION
  /**
   * @name Functions needed for boost::serialization
   * @see <a href="http://www.boost.org/doc/libs/release/libs/serialization/">boost::serialization</a>
   */
  ///@{
  friend class boost::serialization::access;

  template<class Archive, int Size = N>
  void serialize(Archive &ar, const unsigned int version, EnableIf<(Size > 0)> * = 0)
  {
//    ar &boost::serialization::make_nvp("SafeSTLContainer_Array",
//                                       boost::serialization::base_object<base_t>(*this));
    ar &boost::serialization::make_nvp("SafeSTLContainer_Array",
                                       boost::serialization::base_object<std::array<T,N>>(*this));
  }

  template<class Archive, int Size = N>
  void serialize(Archive &ar, const unsigned int version, EnableIf<!(Size > 0)> * = 0)
  {}
  ///@}
#endif // SERIALIZATION

};

IGA_NAMESPACE_CLOSE


#endif // SAFE_STL_ARRAY_H_
