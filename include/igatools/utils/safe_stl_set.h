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

#ifndef __SAFE_STL_SET_H_
#define __SAFE_STL_SET_H_

#include <igatools/base/config.h>
#include <igatools/utils/safe_stl_container.h>
#include <igatools/utils/tensor_index.h>

#include <set>

IGA_NAMESPACE_OPEN




/**
 * @brief iga version of std::set.
 * It can be used as a std::set with print_info support
 */

template<class T>
class SafeSTLSet :
  public SafeSTLContainer<std::set<T>>
{
  using base_t = SafeSTLContainer<std::set<T>>;
public :
  /** Inherit the constructors of the base class. */
  using SafeSTLContainer<std::set<T>>::SafeSTLContainer;

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
    ar &make_nvp("SafeSTLContainer_Set",
                 base_class<std::set<T>>(this));
  }
  ///@}
#endif // SERIALIZATION


};

IGA_NAMESPACE_CLOSE


#ifdef SERIALIZATION
using SafeSTLSetIntAlias = iga::SafeSTLSet<int>;
CEREAL_SPECIALIZE_FOR_ALL_ARCHIVES(SafeSTLSetIntAlias,cereal::specialization::member_serialize);

using SafeSTLSetTensorIndexAlias0 = iga::SafeSTLSet<iga::TensorIndex<0>>;
CEREAL_SPECIALIZE_FOR_ALL_ARCHIVES(SafeSTLSetTensorIndexAlias0,cereal::specialization::member_serialize);
using SafeSTLSetTensorIndexAlias1 = iga::SafeSTLSet<iga::TensorIndex<1>>;
CEREAL_SPECIALIZE_FOR_ALL_ARCHIVES(SafeSTLSetTensorIndexAlias1,cereal::specialization::member_serialize);
using SafeSTLSetTensorIndexAlias2 = iga::SafeSTLSet<iga::TensorIndex<2>>;
CEREAL_SPECIALIZE_FOR_ALL_ARCHIVES(SafeSTLSetTensorIndexAlias2,cereal::specialization::member_serialize);
using SafeSTLSetTensorIndexAlias3 = iga::SafeSTLSet<iga::TensorIndex<3>>;
CEREAL_SPECIALIZE_FOR_ALL_ARCHIVES(SafeSTLSetTensorIndexAlias3,cereal::specialization::member_serialize);

//#include <igatools/utils/safe_stl_set.serialization>
#endif // SERIALIZATION


#endif // SAFE_STL_SET_H_
