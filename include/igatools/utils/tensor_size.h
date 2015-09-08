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

#ifndef __TENSOR_SIZE_H_
#define __TENSOR_SIZE_H_

#include <igatools/base/config.h>
#include <igatools/utils/tensor_index.h>

IGA_NAMESPACE_OPEN

/**
 * @brief Type for the size of a tensor-like container.
 *
 * It is a list of rank number of sizes.
 * The main difference with the TensorIndex is due to the fact that TensorSize
 * has the method flat_size() that returns the multiplication of the sizes
 * along each direction.
 *
 * @ingroup serializable
 *
 * @author M. Martinelli
 * @date 21 Jan 2014
 */
template <int rank>
class TensorSize : public TensorIndex<rank>
{
public:
  /** @name Constructors */
  ///@{

  /** Default constructor. Initializes all the direction indices to zero. */
  explicit TensorSize(const Size value = 0) noexcept ;

  /** Constructor using an SafeSTLArray. */
  explicit TensorSize(const SafeSTLArray<Size,rank> &arr) noexcept;

  /** Copy constructor converting a TensorIndex<rank> @p arr. */
  explicit TensorSize(const TensorIndex<rank> &arr) noexcept;

  /** Constructor using an initializer-list. */
  TensorSize(std::initializer_list<Size> list) noexcept;

  /** Copy constructor. */
  TensorSize(const TensorSize<rank> &arr) = default;

  /** Move constructor. */
  TensorSize(TensorSize<rank> &&arr) = default;

  /** Destructor. */
  ~TensorSize() = default;
  ///@}

  /** @name Assignment operators */
  ///@{
  /**
   * Copy assignment operator
   */
  TensorSize<rank> &operator=(const TensorSize<rank> &arr) = default;

  /**
   * Move assignment operator
   */
  TensorSize<rank> &operator=(TensorSize<rank> &&arr) = default;
  ///@}

  /**
   * Return the flat size, i.e. the multiplications of the sizes along
   * each direction.
   */
  Size flat_size() const noexcept ;

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
    ar &make_nvp("TensorIndex",
                 base_class<TensorIndex<rank>>(this));
  }
  ///@}
#endif // SERIALIZATION
};



/**
 * Returns the tensor size with components given by the difference of
 *  @p index with @p j in all directions.
 *
 * @relates TensorSize
 */
template <int rank>
TensorSize<rank>
operator-(const TensorSize<rank> &index,const Index j);



/**
 * Output operator for TensorIndex.
 *
 * @relates TensorIndex
*/
template <int rank>
LogStream &
operator<<(LogStream &out, const TensorSize<rank> &tensor_size);

IGA_NAMESPACE_CLOSE


#ifdef SERIALIZATION
using TensorSizeAlias0 = iga::TensorSize<0>;
CEREAL_SPECIALIZE_FOR_ALL_ARCHIVES(TensorSizeAlias0,cereal::specialization::member_serialize);
using TensorSizeAlias1 = iga::TensorSize<1>;
CEREAL_SPECIALIZE_FOR_ALL_ARCHIVES(TensorSizeAlias1,cereal::specialization::member_serialize);
using TensorSizeAlias2 = iga::TensorSize<2>;
CEREAL_SPECIALIZE_FOR_ALL_ARCHIVES(TensorSizeAlias2,cereal::specialization::member_serialize);
using TensorSizeAlias3 = iga::TensorSize<3>;
CEREAL_SPECIALIZE_FOR_ALL_ARCHIVES(TensorSizeAlias3,cereal::specialization::member_serialize);

//#include <igatools/utils/tensor_size.serialization>
#endif // SERIALIZATION


#endif // #ifndef TENSOR_SIZE_H_
