//-+--------------------------------------------------------------------
// Igatools a general purpose Isogeometric analysis library.
// Copyright (C) 2012-2016  by the igatools authors (see authors.txt).
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

template <class> class SafeSTLVector;
template <class,int> class DynamicMultiArray;
template <int> class TensorSize;
template <int> class TensorIndex;
template <class,int> class CartesianProductArray;

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
  SafeSTLArray()
  {};

  /**
   * Constructor. It fills the <tt>N</tt> entries of the array with <tt>N</tt> copies
   * of the input argument <tt>val</tt>.
   */
  explicit SafeSTLArray(const T &val)
  {
    std::array<T,N>::fill(val);
  };


  SafeSTLArray(std::initializer_list<T> list)
  {
    Assert(list.size() == N, ExcDimensionMismatch(list.size(),N));

    for (int i = 0 ; i < N ; ++i)
      (*this)[i] = list.begin()[i] ;
  };

  /**
   * Returns a reference to the <tt>n</tt>-th entry of the container.
   * @note In Debug mode the value of <tt>n</tt> is checked if within the valid bounds of the container.
   */
  T &operator[](Size n)
  {
    Assert(n < this->size(), ExcIndexRange(n, 0, this->size()));
    return std::array<T,N>::operator[](n);
  }

  /**
   * Returns a const-reference to the <tt>n</tt>-th entry of the container.
   * @note In Debug mode the value of <tt>n</tt> is checked if within the valid bounds of the container.
   */
  const T &operator[](Size n) const
  {
    Assert(n < this->size(), ExcIndexRange(n, 0, this->size()));
    return std::array<T,N>::operator[](n);
  }

private:

#ifdef IGATOOLS_WITH_SERIALIZATION
  /**
   * @name Functions needed for serialization
   * @see <a href="http://uscilab.github.io/cereal/serialization_functions.html">Cereal serialization</a>
   */
  ///@{
  friend class cereal::access;

  template<class Archive>
  void serialize(Archive &ar)
  {
    ar &make_nvp("SafeSTLContainer_Array",
                 base_class<std::array<T,N>>(this));
  }
  ///@}
#endif // IGATOOLS_WITH_SERIALIZATION

};


IGA_NAMESPACE_CLOSE



#ifdef IGATOOLS_WITH_SERIALIZATION
#include <igatools/utils/safe_stl_array.serial>
#endif // IGATOOLS_WITH_SERIALIZATION


#endif // SAFE_STL_ARRAY_H_
