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


using VecReal = iga::SafeSTLVector<iga::Real>;
using SafeSTLArrayVecRealAlias0 = iga::SafeSTLArray<VecReal,0>;
CEREAL_SPECIALIZE_FOR_ALL_ARCHIVES(SafeSTLArrayVecRealAlias0,cereal::specialization::member_serialize);
using SafeSTLArrayVecRealAlias1 = iga::SafeSTLArray<VecReal,1>;
CEREAL_SPECIALIZE_FOR_ALL_ARCHIVES(SafeSTLArrayVecRealAlias1,cereal::specialization::member_serialize);
using SafeSTLArrayVecRealAlias2 = iga::SafeSTLArray<VecReal,2>;
CEREAL_SPECIALIZE_FOR_ALL_ARCHIVES(SafeSTLArrayVecRealAlias2,cereal::specialization::member_serialize);
using SafeSTLArrayVecRealAlias3 = iga::SafeSTLArray<VecReal,3>;
CEREAL_SPECIALIZE_FOR_ALL_ARCHIVES(SafeSTLArrayVecRealAlias3,cereal::specialization::member_serialize);

// The next ones are used by SplineSpace
using VecInt = iga::SafeSTLVector<int>;
//using SafeSTLArrayVecIntAlias0 = iga::SafeSTLArray<VecInt,0>;
//CEREAL_SPECIALIZE_FOR_ALL_ARCHIVES(SafeSTLArrayVecIntAlias0,cereal::specialization::member_serialize);
using SafeSTLArrayVecIntAlias1 = iga::SafeSTLArray<VecInt,1>;
CEREAL_SPECIALIZE_FOR_ALL_ARCHIVES(SafeSTLArrayVecIntAlias1,cereal::specialization::member_serialize);
using SafeSTLArrayVecIntAlias2 = iga::SafeSTLArray<VecInt,2>;
CEREAL_SPECIALIZE_FOR_ALL_ARCHIVES(SafeSTLArrayVecIntAlias2,cereal::specialization::member_serialize);
using SafeSTLArrayVecIntAlias3 = iga::SafeSTLArray<VecInt,3>;
CEREAL_SPECIALIZE_FOR_ALL_ARCHIVES(SafeSTLArrayVecIntAlias3,cereal::specialization::member_serialize);


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
using SafeSTLArrayTSizeAlias0_1 = iga::SafeSTLArray<iga::TensorSize<0>,1>;
CEREAL_SPECIALIZE_FOR_ALL_ARCHIVES(SafeSTLArrayTSizeAlias0_1,cereal::specialization::member_serialize);
using SafeSTLArrayTSizeAlias1_1 = iga::SafeSTLArray<iga::TensorSize<1>,1>;
CEREAL_SPECIALIZE_FOR_ALL_ARCHIVES(SafeSTLArrayTSizeAlias1_1,cereal::specialization::member_serialize);
using SafeSTLArrayTSizeAlias2_1 = iga::SafeSTLArray<iga::TensorSize<2>,1>;
CEREAL_SPECIALIZE_FOR_ALL_ARCHIVES(SafeSTLArrayTSizeAlias2_1,cereal::specialization::member_serialize);
using SafeSTLArrayTSizeAlias3_1 = iga::SafeSTLArray<iga::TensorSize<3>,1>;
CEREAL_SPECIALIZE_FOR_ALL_ARCHIVES(SafeSTLArrayTSizeAlias3_1,cereal::specialization::member_serialize);
using SafeSTLArrayTSizeAlias2_2 = iga::SafeSTLArray<iga::TensorSize<2>,2>;
CEREAL_SPECIALIZE_FOR_ALL_ARCHIVES(SafeSTLArrayTSizeAlias2_2,cereal::specialization::member_serialize);



// The next ones are used by SplineSpace
using SafeSTLArrayBoolAlias1 = iga::SafeSTLArray<bool,1>;
CEREAL_SPECIALIZE_FOR_ALL_ARCHIVES(SafeSTLArrayBoolAlias1,cereal::specialization::member_serialize);
using SafeSTLArrayBoolAlias2 = iga::SafeSTLArray<bool,2>;
CEREAL_SPECIALIZE_FOR_ALL_ARCHIVES(SafeSTLArrayBoolAlias2,cereal::specialization::member_serialize);
using SafeSTLArrayBoolAlias3 = iga::SafeSTLArray<bool,3>;
CEREAL_SPECIALIZE_FOR_ALL_ARCHIVES(SafeSTLArrayBoolAlias3,cereal::specialization::member_serialize);


// The next ones are used by SplineSpace
using SafeSTLArrayBoolArrayAlias1_1 = iga::SafeSTLArray<iga::SafeSTLArray<bool,1>,1>;
CEREAL_SPECIALIZE_FOR_ALL_ARCHIVES(SafeSTLArrayBoolArrayAlias1_1,cereal::specialization::member_serialize);
using SafeSTLArrayBoolArrayAlias2_1 = iga::SafeSTLArray<iga::SafeSTLArray<bool,2>,1>;
CEREAL_SPECIALIZE_FOR_ALL_ARCHIVES(SafeSTLArrayBoolArrayAlias2_1,cereal::specialization::member_serialize);
using SafeSTLArrayBoolArrayAlias3_1 = iga::SafeSTLArray<iga::SafeSTLArray<bool,3>,1>;
CEREAL_SPECIALIZE_FOR_ALL_ARCHIVES(SafeSTLArrayBoolArrayAlias3_1,cereal::specialization::member_serialize);
using SafeSTLArrayBoolArrayAlias2_2 = iga::SafeSTLArray<iga::SafeSTLArray<bool,2>,2>;
CEREAL_SPECIALIZE_FOR_ALL_ARCHIVES(SafeSTLArrayBoolArrayAlias2_2,cereal::specialization::member_serialize);



// The next ones are used by SplineSpace
using SafeSTLVecTIAlias1 = iga::SafeSTLVector<iga::TensorIndex<1>>;
CEREAL_SPECIALIZE_FOR_ALL_ARCHIVES(SafeSTLVecTIAlias1,cereal::specialization::member_serialize);
using SafeSTLVecTIAlias2 = iga::SafeSTLVector<iga::TensorIndex<2>>;
CEREAL_SPECIALIZE_FOR_ALL_ARCHIVES(SafeSTLVecTIAlias2,cereal::specialization::member_serialize);
using SafeSTLVecTIAlias3 = iga::SafeSTLVector<iga::TensorIndex<3>>;
CEREAL_SPECIALIZE_FOR_ALL_ARCHIVES(SafeSTLVecTIAlias3,cereal::specialization::member_serialize);


// The next ones are used by SplineSpace
template <int N>
using VecTI = iga::SafeSTLVector<iga::TensorIndex<N>>;

using SafeSTLArrayVecTIAlias1_1 = iga::SafeSTLArray<VecTI<1>,1>;
CEREAL_SPECIALIZE_FOR_ALL_ARCHIVES(SafeSTLArrayVecTIAlias1_1,cereal::specialization::member_serialize);
using SafeSTLArrayVecTIAlias2_1 = iga::SafeSTLArray<VecTI<2>,1>;
CEREAL_SPECIALIZE_FOR_ALL_ARCHIVES(SafeSTLArrayVecTIAlias2_1,cereal::specialization::member_serialize);
using SafeSTLArrayVecTIAlias3_1 = iga::SafeSTLArray<VecTI<3>,1>;
CEREAL_SPECIALIZE_FOR_ALL_ARCHIVES(SafeSTLArrayVecTIAlias3_1,cereal::specialization::member_serialize);
using SafeSTLArrayVecTIAlias2_2 = iga::SafeSTLArray<VecTI<2>,2>;
CEREAL_SPECIALIZE_FOR_ALL_ARCHIVES(SafeSTLArrayVecTIAlias2_2,cereal::specialization::member_serialize);


// The next ones are used by SplineSpace
using SafeSTLArrayTIndexAlias1_1 = iga::SafeSTLArray<iga::TensorIndex<1>,1>;
CEREAL_SPECIALIZE_FOR_ALL_ARCHIVES(SafeSTLArrayTIndexAlias1_1,cereal::specialization::member_serialize);
using SafeSTLArrayTIndexAlias2_1 = iga::SafeSTLArray<iga::TensorIndex<2>,1>;
CEREAL_SPECIALIZE_FOR_ALL_ARCHIVES(SafeSTLArrayTIndexAlias2_1,cereal::specialization::member_serialize);
using SafeSTLArrayTIndexAlias3_1 = iga::SafeSTLArray<iga::TensorIndex<3>,1>;
CEREAL_SPECIALIZE_FOR_ALL_ARCHIVES(SafeSTLArrayTIndexAlias3_1,cereal::specialization::member_serialize);
using SafeSTLArrayTIndexAlias2_2 = iga::SafeSTLArray<iga::TensorIndex<2>,2>;
CEREAL_SPECIALIZE_FOR_ALL_ARCHIVES(SafeSTLArrayTIndexAlias2_2,cereal::specialization::member_serialize);



// The next ones are used by SplineSpace
//using SafeSTLArrayCPArrayIntAlias0 = iga::SafeSTLArray<iga::CartesianProductArray<int,0>,1>;
//CEREAL_SPECIALIZE_FOR_ALL_ARCHIVES(SafeSTLArrayCPArrayIntAlias0,cereal::specialization::member_serialize);
using SafeSTLArrayCPArrayIntAlias1_1 = iga::SafeSTLArray<iga::CartesianProductArray<int,1>,1>;
CEREAL_SPECIALIZE_FOR_ALL_ARCHIVES(SafeSTLArrayCPArrayIntAlias1_1,cereal::specialization::member_serialize);
using SafeSTLArrayCPArrayIntAlias2_1 = iga::SafeSTLArray<iga::CartesianProductArray<int,2>,1>;
CEREAL_SPECIALIZE_FOR_ALL_ARCHIVES(SafeSTLArrayCPArrayIntAlias2_1,cereal::specialization::member_serialize);
using SafeSTLArrayCPArrayIntAlias3_1 = iga::SafeSTLArray<iga::CartesianProductArray<int,3>,1>;
CEREAL_SPECIALIZE_FOR_ALL_ARCHIVES(SafeSTLArrayCPArrayIntAlias3_1,cereal::specialization::member_serialize);
using SafeSTLArrayCPArrayIntAlias2_2 = iga::SafeSTLArray<iga::CartesianProductArray<int,2>,2>;
CEREAL_SPECIALIZE_FOR_ALL_ARCHIVES(SafeSTLArrayCPArrayIntAlias2_2,cereal::specialization::member_serialize);

//#include <igatools/utils/safe_stl_array.serialization>
#endif // SERIALIZATION


#endif // SAFE_STL_ARRAY_H_
