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

// QualityAssurance: martinelli, 28 Jan 2014

#ifndef TENSOR_SIZED_CONTAINER_H_
#define TENSOR_SIZED_CONTAINER_H_

#include <igatools/base/config.h>
#include <igatools/utils/tensor_size.h>
#include <igatools/utils/tensor_index.h>

IGA_NAMESPACE_OPEN


/**
 * @brief This class represent a rank-dimensional tensor-sized container.
 *
 * It is intended to be used as base class for tensor-like containers,
 * e.g. MultiArray, StaticMultiArray, CartesianProductArray, TensorProductArray, etc.
 * in order to have an unified treatment for the size information and for
 * the index transformations through the member functions
 * flat_to_tensor() and tensor_to_flat().
 *
 * @author M. Martinelli
 * @date 27 Jan 2014
 *
 * @ingroup serializable
 *
 */
template <int rank>
class TensorSizedContainer
{
public:
  /** @name Constructors and destructor */
  ///@{

  /** Default constructor. Sets the size of the container to be 0 in each dimension. */
  TensorSizedContainer();

  /**
   * Constuctor. Sets the size of the container with size (in each dimension) specified by
   * the input argument @p size.
   */
  explicit TensorSizedContainer(const TensorSize<rank> &size);

  /**
   * Constuctor. Sets the size of the container with the same size in each dimension,
   *  specified by the input argument @p size.
   */
  explicit TensorSizedContainer(const Size size);

  /** Copy constructor. */
  TensorSizedContainer(const TensorSizedContainer<rank> &in) = default;

  /** Move constructor. */
  TensorSizedContainer(TensorSizedContainer<rank> &&in) = default;

  /** Destructor.*/
  ~TensorSizedContainer() = default;
  ///@}



  /** @name Assignment operators */
  ///@{
  /** Copy assignment operator .*/
  TensorSizedContainer<rank> &operator=(const TensorSizedContainer<rank> &in) = default;

  /** Move assignment operator .*/
  TensorSizedContainer<rank> &operator=(TensorSizedContainer<rank> &&in) = default;
  ///@}

  /** @name Functions for getting the informations about the size of the container */
  ///@{
  /** Return the tensor-size of the container. */
  TensorSize<rank> tensor_size() const;


  /**
   * Return the flat-size (i.e. the total size) of the container.
   *
   * If the direction size are for example (2,5,3) this
   * function would return 2x5x3 = 30.
   */
  Size flat_size() const;
  ///@}

  void print_info(LogStream &out) const;

  /** @name Functions for the index transformations */
  ///@{
  bool valid_index(const TensorIndex<rank> &tensor_index) const;

  /**
   * Transformation from a tensor-index to a flat-index.
   */
  Index tensor_to_flat(const TensorIndex<rank> &tensor_index) const;

  /**
   * Transformation from a flat-index to a tensor-index.
   */
  TensorIndex<rank> flat_to_tensor(const Index flat_index) const;
  ///@}


  /**
   * Returns the memory used to define the object's data members.
   */
  auto memory_consumption() const
  {
    return size_.memory_consumption() + weight_.memory_consumption();
  }

protected:
  /**
   * Reset the size_ member in order to be equal to the input argument @p size.
   *
   * @note This function is protected because it should be used in the resize() functions
   * of a derived class.
   */
  void reset_size(const TensorSize<rank> &size);


private:
  /** Size of the container along the different directions. */
  TensorSize<rank> size_;

  /**
   * Weights for the index conversion.
   */
  TensorIndex<rank> weight_;


#ifdef SERIALIZATION
  /**
   * @name Functions needed for boost::serialization
   * @see <a href="http://www.boost.org/doc/libs/release/libs/serialization/">boost::serialization</a>
   */
  ///@{
  friend class boost::serialization::access;

  template<class Archive>
  void serialize(Archive &ar, const unsigned int version)
  {
    ar &make_nvp("size_",size_);

    ar &make_nvp("weight_",weight_);
  }

  ///@}
#endif // SERIALIZATION
};


IGA_NAMESPACE_CLOSE


#ifdef NDEBUG
#include <igatools/utils/tensor_sized_container-inline.h>
#endif

#endif // #ifndef TENSOR_SIZED_CONTAINER_H_
