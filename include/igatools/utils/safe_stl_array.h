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

template <class> class SafeSTLVector;
template <class,int> class DynamicMultiArray;
template <int> class TensorSize;


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
#endif // SERIALIZATION

};

IGA_NAMESPACE_CLOSE



#ifdef SERIALIZATION
using SafeSTLArrayAliasInt0 = iga::SafeSTLArray<int,0>;
CEREAL_SPECIALIZE_FOR_ALL_ARCHIVES(SafeSTLArrayAliasInt0,cereal::specialization::member_serialize);
using SafeSTLArrayAliasInt1 = iga::SafeSTLArray<int,1>;
CEREAL_SPECIALIZE_FOR_ALL_ARCHIVES(SafeSTLArrayAliasInt1,cereal::specialization::member_serialize);
using SafeSTLArrayAliasInt2 = iga::SafeSTLArray<int,2>;
CEREAL_SPECIALIZE_FOR_ALL_ARCHIVES(SafeSTLArrayAliasInt2,cereal::specialization::member_serialize);
using SafeSTLArrayAliasInt3 = iga::SafeSTLArray<int,3>;
CEREAL_SPECIALIZE_FOR_ALL_ARCHIVES(SafeSTLArrayAliasInt3,cereal::specialization::member_serialize);
using SafeSTLArrayAliasInt4 = iga::SafeSTLArray<int,4>;
CEREAL_SPECIALIZE_FOR_ALL_ARCHIVES(SafeSTLArrayAliasInt4,cereal::specialization::member_serialize);
using SafeSTLArrayAliasInt6 = iga::SafeSTLArray<int,6>;
CEREAL_SPECIALIZE_FOR_ALL_ARCHIVES(SafeSTLArrayAliasInt6,cereal::specialization::member_serialize);


using Vec = iga::SafeSTLVector<iga::Real>;
using SafeSTLArrayAliasVec0 = iga::SafeSTLArray<Vec,0>;
CEREAL_SPECIALIZE_FOR_ALL_ARCHIVES(SafeSTLArrayAliasVec0,cereal::specialization::member_serialize);
using SafeSTLArrayAliasVec1 = iga::SafeSTLArray<Vec,1>;
CEREAL_SPECIALIZE_FOR_ALL_ARCHIVES(SafeSTLArrayAliasVec1,cereal::specialization::member_serialize);
using SafeSTLArrayAliasVec2 = iga::SafeSTLArray<Vec,2>;
CEREAL_SPECIALIZE_FOR_ALL_ARCHIVES(SafeSTLArrayAliasVec2,cereal::specialization::member_serialize);
using SafeSTLArrayAliasVec3 = iga::SafeSTLArray<Vec,3>;
CEREAL_SPECIALIZE_FOR_ALL_ARCHIVES(SafeSTLArrayAliasVec3,cereal::specialization::member_serialize);


// The next ones are used by DofDistribution
using SafeSTLArrayDMArrayAliasVec0 = iga::SafeSTLArray<iga::DynamicMultiArray<int,0>,1>;
CEREAL_SPECIALIZE_FOR_ALL_ARCHIVES(SafeSTLArrayDMArrayAliasVec0,cereal::specialization::member_serialize);
using SafeSTLArrayDMArrayAliasVec1 = iga::SafeSTLArray<iga::DynamicMultiArray<int,1>,1>;
CEREAL_SPECIALIZE_FOR_ALL_ARCHIVES(SafeSTLArrayDMArrayAliasVec1,cereal::specialization::member_serialize);
using SafeSTLArrayDMArrayAliasVec2 = iga::SafeSTLArray<iga::DynamicMultiArray<int,2>,1>;
CEREAL_SPECIALIZE_FOR_ALL_ARCHIVES(SafeSTLArrayDMArrayAliasVec2,cereal::specialization::member_serialize);
using SafeSTLArrayDMArrayAliasVec3 = iga::SafeSTLArray<iga::DynamicMultiArray<int,3>,1>;
CEREAL_SPECIALIZE_FOR_ALL_ARCHIVES(SafeSTLArrayDMArrayAliasVec3,cereal::specialization::member_serialize);

// The next ones are used by SplineSpace::ComponentContainer
using SafeSTLArrayTSizeAlias0 = iga::SafeSTLArray<iga::TensorSize<0>,1>;
CEREAL_SPECIALIZE_FOR_ALL_ARCHIVES(SafeSTLArrayTSizeAlias0,cereal::specialization::member_serialize);
using SafeSTLArrayTSizeAlias1 = iga::SafeSTLArray<iga::TensorSize<1>,1>;
CEREAL_SPECIALIZE_FOR_ALL_ARCHIVES(SafeSTLArrayTSizeAlias1,cereal::specialization::member_serialize);
using SafeSTLArrayTSizeAlias2 = iga::SafeSTLArray<iga::TensorSize<2>,1>;
CEREAL_SPECIALIZE_FOR_ALL_ARCHIVES(SafeSTLArrayTSizeAlias2,cereal::specialization::member_serialize);
using SafeSTLArrayTSizeAlias3 = iga::SafeSTLArray<iga::TensorSize<3>,1>;
CEREAL_SPECIALIZE_FOR_ALL_ARCHIVES(SafeSTLArrayTSizeAlias3,cereal::specialization::member_serialize);


//#include <igatools/utils/safe_stl_array.serialization>
#endif // SERIALIZATION


#endif // SAFE_STL_ARRAY_H_
