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

#ifndef __TENSOR_INDEX_H_
#define __TENSOR_INDEX_H_

#include <igatools/base/config.h>
#include <igatools/base/logstream.h>
#include <igatools/utils/safe_stl_vector.h>
#include <igatools/base/array_utils.h>

IGA_NAMESPACE_OPEN

/**
 * @brief Type for a tensor index.
 *
 * It is a list of rank number of non-negative indices.
 *
 * This class makes possible the
 * rank independent treatment of tensor type containers.
 *
 * @ingroup serializable
 *
 * @author martinelli 2014
 * @author pauletti 2015
 */
template <int rank>
class TensorIndex : public SafeSTLArray<Index, rank>
{
public:
  /** @name Constructors */
  ///@{

  /**
   * Default constructor. Initializes all the direction indices to
   * @p value. If no @p value is provided, the initialization is made with zeros.
   */
  explicit TensorIndex(const Size value = 0) noexcept ;

  /** Constructor using an SafeSTLArray. */
  explicit TensorIndex(const SafeSTLArray<Index,rank> &arr) noexcept;

  /** Constructor using an initializer-list. */
  TensorIndex(std::initializer_list<Index> list) noexcept;

  /** Copy constructor. */
  TensorIndex(const TensorIndex<rank> &arr) = default;

  /** Move constructor. */
  TensorIndex(TensorIndex<rank> &&arr) = default;

  /** Destructor. */
  ~TensorIndex() = default;
  ///@}

  /** @name Assignment operators */
  ///@{

  /**
   * Copy assignment operator
   */
  TensorIndex<rank> &operator=(const TensorIndex<rank> &arr) = default;

  /**
   * Move assignment operator
   */
  TensorIndex<rank> &operator=(TensorIndex<rank> &&arr) = default;

  ///@}

  /**
   * Returns a rank k tensor index formed by the component
   * in the k indices of sub_indices
   */
  template<int k>
  TensorIndex<k> get_sub_tensor(const TensorIndex<k> &sub_indices) const
  {
    TensorIndex<k> res;
    int j = 0;
    for (auto i : sub_indices)
      res[j++] = (*this)[i];

    return res;
  }

  /**
   * @name Increment/decrement operators.
   */
  ///@{
  /** Increment the tensor index by the values @p idx.*/
  TensorIndex<rank> &operator +=(const TensorIndex<rank> &idx) noexcept ;

  /** Decrement the tensor index by the values @p idx.*/
  TensorIndex<rank> &operator -=(const TensorIndex<rank> &idx) noexcept ;

  /** Increment the indices in all directions by the value @p j. */
  TensorIndex<rank> &operator +=(const Index j) noexcept ;

  /** Decrement the indices in all directions by the value @p j. */
  TensorIndex<rank> &operator -=(const Index j) noexcept ;
  ///@}

  /**
   * Returns the memory used to define the object's data members.
   */
  std::size_t memory_consumption() const;

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
    ar &make_nvp("SafeSTLArray",
                 base_class<SafeSTLArray<Index,rank>>(this));
  }
  ///@}
#endif // SERIALIZATION
};


/**
 * Returns the tensor index with components given by the sum of @p index_a with @p index_b.
 *
 * @relates TensorIndex
 */
template <int rank>
TensorIndex<rank>
operator+(const TensorIndex<rank> &index_a,const TensorIndex<rank> &index_b) ;

/**
 * Returns the tensor index with components given by the sum of @p index with @p j in all directions.
 *
 * @relates TensorIndex
 */
template <int rank>
TensorIndex<rank>
operator+(const TensorIndex<rank> &index,const Index j) ;

/**
 * Returns the tensor index with components given by the difference of
 *  @p index with @p j in all directions.
 *
 * @relates TensorIndex
 */
template <int rank>
TensorIndex<rank>
operator-(const TensorIndex<rank> &index,const Index j) ;

#if 0 //already implemente in std::array
template <int rank>
bool
operator<(const TensorIndex<rank> &index_a,const TensorIndex<rank> &index_b)
{
  for (int j=0; j<rank; ++j)
    if (index_a[j] >= index_b[j])
      return false;
  return true;
}


template <int rank>
bool
operator<=(const TensorIndex<rank> &index_a,const TensorIndex<rank> &index_b)
{
  for (int j=0; j<rank; ++j)
    if (index_a[j] > index_b[j])
      return false;
  return true;
}


template <int rank>
bool
operator>(const TensorIndex<rank> &index_a,const TensorIndex<rank> &index_b)
{
  return (!(index_a <= index_b));
}


template <int rank>
bool
operator>=(const TensorIndex<rank> &index_a,const TensorIndex<rank> &index_b)
{
  return (!(index_a < index_b));
}
#endif
/**
 * Output operator for TensorIndex.
 *
 * @relates TensorIndex
*/
template <int rank>
LogStream &
operator<<(LogStream &out, const TensorIndex<rank> &tensor_index) ;

#if 0
/**
 * Generates a vector with the tensor indices of the given
 * rectangular range.
 *
 *  @relates TensorIndex
 *  @author pauletti 2015
 */
template<int k>
SafeSTLVector<TensorIndex<k>> tensor_range(TensorIndex<k> first, TensorIndex<k> last)
{
  Assert(first <= last, ExcMessage("first bigger than last"));
  SafeSTLVector<TensorIndex<k>> result;
  TensorIndex<k-1> ind(sequence<k-1>());
  auto vec = tensor_range<k-1>(first.get_sub_tensor(ind), last.get_sub_tensor(ind));

  for (int i=first[k-1]; i<last[k-1]; ++i)
  {
    for (auto &t_k_1 : vec)
    {
      TensorIndex<k> t_k;
      for (int j=0; j<k-1; ++j)
        t_k[j] = t_k_1[j];
      t_k[k-1] = i;
      result.insert(t_k);
    }
  }
  return result;
}

template<>
SafeSTLVector<TensorIndex<1>> tensor_range(TensorIndex<1> first, TensorIndex<1> last);

template<>
SafeSTLVector<TensorIndex<0>> tensor_range(TensorIndex<0> first, TensorIndex<0> last);
#endif

IGA_NAMESPACE_CLOSE


#ifdef NDEBUG
#include <igatools/utils/tensor_index-inline.h>
#endif



#ifdef SERIALIZATION
using TensorIndexAlias0 = iga::TensorIndex<0>;
CEREAL_SPECIALIZE_FOR_ALL_ARCHIVES(TensorIndexAlias0,cereal::specialization::member_serialize)
using TensorIndexAlias1 = iga::TensorIndex<1>;
CEREAL_SPECIALIZE_FOR_ALL_ARCHIVES(TensorIndexAlias1,cereal::specialization::member_serialize)
using TensorIndexAlias2 = iga::TensorIndex<2>;
CEREAL_SPECIALIZE_FOR_ALL_ARCHIVES(TensorIndexAlias2,cereal::specialization::member_serialize)
using TensorIndexAlias3 = iga::TensorIndex<3>;
CEREAL_SPECIALIZE_FOR_ALL_ARCHIVES(TensorIndexAlias3,cereal::specialization::member_serialize)

//#include <igatools/utils/tensor_index.serialization>
#endif // SERIALIZATION



#endif // TENSOR_INDEX_H_



