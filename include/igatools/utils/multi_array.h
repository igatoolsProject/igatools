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


// QualityAssurance: martinelli, 31 Jan 2014

#ifndef MULTI_ARRAY_H_
#define MULTI_ARRAY_H_

#include <igatools/base/config.h>
#include <igatools/base/print_info_utils.h>
#include <igatools/utils/tensor_sized_container.h>
#include <igatools/utils/multi_array_iterator.h>
#include <igatools/utils/container_view.h>


IGA_NAMESPACE_OPEN

/**
 * @brief Base class representing a multi-array, i.e. a tensor-like array
 * container of fixed or variable dimension and fixed @p rank.
 *
 * This container represent collections of objects that can be referred with
 * \f$d\f$-dimensional multi-indices, e.g.
 * \f[
 *   \mathbf{W} = \{ \mathbf{w}_{i_1,\dots,i_d} \in \mathcal{W},
 *   \quad i_k=1,\dots,n_k , \quad k=1,\dots,d \} \;.
 * \f]
 * The type of container for storing the multi-array entries is specified by the
 * template parameter @p STLContainer (and therefore it specifies implicilty the type of the entries).
 * This permits to use the MultiArray class as base
 * class for static and dynamic multi-arrays (StaticMultiArray and DynamicMultiArray respectively).
 *
 * Moreover, this permits to unify the implementation the entry access operators and the
 * functions returning the multi-array iterator.
 *
 * Due to the generality of the underlying @p STLContainer,
 * we provide the entries access operators () (const and non-const version)
 * accepting both multi-indices (through its representation via the class TensorIndex)
 * and flat indices.
 *
 * Its main use is in the range independent treatment through
 * the use of the flat index.
 *
 * - @p rank == 0 is a single element
 * - @p rank == 1 is a vector
 * - @p rank == 2 a tensor
 *
 *
 * ## Iterator
 *
 * For the MultiArray class it is also defined
 * <a href="http://www.cplusplus.com/reference/iterator/RandomAccessIterator/">random access iterator</a> called MultiArrayIterator,
 * that can be used to access the MultiArray entries with iterator techniques
 * (e.g. the range-for loop).
 *
 * In order to get the iterator pointing to the first entry of the container you can use the
 * MultiArray::begin() functions, returning the iterator in the const or non-const version,
 * or the function MultiArray::cbegin(), returning a const iterator. To get the iterator pointing to
 * one-pass the end of the container, you can use MultiArray::end() (const and non-const version) or
 * MultiArray::cend() (const version).
 *
 *
 * ## View
 * Another way to access the elements of a MultiArray is using a ContainerView
 * (or its <tt>const</tt> counterpart ConstContainerView).
 * A ContainerView is basically an object that:
 * - store two iterators (one pointing the begin of the container
 * and the other pointing to one-pass-end the container;
 * - have some methods for accessing the entries of the container.
 * It is important to note that a ContainerView is a lightweight object
 * (compared to the container from which the View refers from) because the memory used to store
 * the data is allocated on the container and not on the ContainerView,
 * and that the entry access are managed by the iterators.
 *
 * @ingroup multi_array_containers
 * @ingroup serializable
 * @author M. Martinelli
 * @date 2014, 2015
 */
//TODO(pauletti, May 31, 2014): we should provide a tensorindex type for loop
template<class STLContainer, int rank>
class MultiArray : public TensorSizedContainer<rank>
{
public:
  /** Type of the entries stored in the STL container. */
  using value_type = typename STLContainer::value_type;

  /** Type for the reference in the STL container. */
  using reference = typename STLContainer::reference;

  /** Type for the const_reference in the STL container. */
  using const_reference = typename STLContainer::const_reference;

  /** @name Constructors and destructor */
  ///@{
  /**
   * Construct an empty multiarray.
   */
  MultiArray();

  /**
   * Construct a square multiarray of zeros with @p dim entries in each array dimension.
   */
  MultiArray(const Size dim);

  /**
   * Construct a rectangular multiarray of zeros with @p dim[i] entries in the i-th array dimension.
   */
  MultiArray(const TensorSize<rank> &dim);

  /** Copy constructor. */
  MultiArray(const MultiArray<STLContainer,rank> &data) = default;

  /** Move constructor. */
  MultiArray(MultiArray<STLContainer,rank> &&data) = default;

  /** Destructor. */
  ~MultiArray() = default;
  ///@}


  /** @name Assignment operators */
  ///@{
  /** Copy assignment operator. */
  MultiArray<STLContainer,rank> &operator=(const MultiArray<STLContainer,rank> &data) = default;

  /** Move assignment operator. */
  MultiArray<STLContainer,rank> &operator=(MultiArray<STLContainer,rank> &&data) = default;
  ///@}

  /** @name Access operators */
  ///@{

  /**
   * Flat index access operator (non-const version).
   * @note In Debug mode the index @p i is checked in order to be
   * in the bounds of the container.
   */
  reference operator[](const Index i);

  /**
   * Flat index access operator (const version).
   * @note In Debug mode the index @p i is checked in order to be
   * in the bounds of the container.
   */
  const_reference operator[](const Index i) const;


  /**
   *  Tensor index access operator (non-const version).
   */
  reference operator()(const TensorIndex<rank> &i);

  /**
   *  Tensor index access operator (const version).
   */
  const_reference operator()(const TensorIndex<rank> &i) const;

  /** Return the entries of the multiarray as unidimensional STLContainer. */
  const STLContainer &get_data() const;
  ///@}


  /** @name Dealing with the iterators */
  ///@{

  /** Type of the const iterator. */
  using const_iterator = MultiArrayConstIterator<MultiArray<STLContainer,rank>>;

  /** Type of the iterator. */
  using iterator = MultiArrayIterator<MultiArray<STLContainer,rank>>;


  /** Returns a const_iterator pointing to the first element in the container. */
  const_iterator cbegin() const;

  /** Returns a const_iterator pointing to the to one-pass the end in the container. */
  const_iterator cend() const;

  /** Returns a const_iterator pointing to the first element in the container. */
  const_iterator begin() const;

  /** Returns a const_iterator pointing to the to one-pass the end in the container. */
  const_iterator end() const;

  /** Returns an iterator pointing to the first element in the container. */
  iterator begin();

  /** Returns an iterator pointing to the to one-pass the end in the container. */
  iterator end();

  ///@}


  /** @name Getting a view of the data*/
  ///@{
  /** Returns a ContainerView of the MultiArray. */
  ContainerView<MultiArray<STLContainer,rank>> get_view();

  /** Returns a ConstContainerView of the MultiArray. */
  ConstContainerView<MultiArray<STLContainer,rank>> get_const_view() const;

  /**
   * Returns a ContainerView of the underlying STLContainer
   * used to store the entries of the MultiArray.
   */
  ContainerView<STLContainer> get_flat_view();

  /**
   * Returns a ConstContainerView of the underlying STLContainer
   * used to store the entries of the MultiArray.
   */
  ConstContainerView<STLContainer> get_flat_const_view() const;
  ///@}

  /** @name Functions to easily fill the multiarray */
  ///@{
  /**
   * Fills the multiarray with an arithmetic progression starting with the value @p init_value.
   */
  void fill_progression(const value_type &init_value = {});


  /** Fills the multiarray copying in each entry the content of @p value. */
  void fill(const value_type &value);
  ///@}


  /**
   * @name Printing info
   */
  ///@{
private:
  template <class A>
  EnableIf<has_print_info<A>(0), void>
  t_print_info(LogStream &out) const
  {
    data_.print_info(out);
  }

  template <class A>
  EnableIf<(!has_print_info<A>(0)), void>
  t_print_info(LogStream &out) const
  {
    for (const auto &v : data_)
      out << v << " ";
  }
public:
  /**
   * Prints the content of the MultiArray on the LogStream @p out.
   * Its use is intended mainly for testing and debugging purpose.
   */
  void print_info(LogStream &out) const;
  ///@}


  /**
   * Returns an estimate of the memory used to define the object.
   */
  std::size_t memory_consumption() const;

protected:
  /**
   * @name Functions for changing the size or the shape of the MultiArray
   * (it makes sense only if STLContainer is resizable).
   */
  ///@{

  /**
   * Resize the MultiArray as square container with @p dim entries in each
   * array dimension.
   */
  void resize(const Size dim);


  /**
   * Resize the MultiArray as rectangular container with <p>dim[i]</p> entries
   * in the i-th array dimension.
   */
  void resize(const TensorSize<rank> &dim);


  /**
   * Resize the MultiArray as rectangular container with <p>dim[i]</p> entries
   * in the i-th array dimension.
   *
   * All the entries are initialized to the value @p val.
   */
  void resize(const TensorSize<rank> &dim, const typename STLContainer::value_type &val);
  ///@}


private:

  /**
   * Data of type Entry stored in a STL container.
   */
  STLContainer data_;



#ifdef IGATOOLS_WITH_SERIALIZATION
  /**
   * @name Functions needed for boost::serialization
   * @see <a href="http://www.boost.org/doc/libs/release/libs/serialization/">boost::serialization</a>
   */
  ///@{
  friend class cereal::access;

  template<class Archive>
  void
  serialize(Archive &ar)
  {
    ar &make_nvp("MultiArray_base_t",
                 base_class<TensorSizedContainer<rank>>(this));

    ar &make_nvp("data_",data_);
  }
  ///@}
#endif // IGATOOLS_WITH_SERIALIZATION
};


IGA_NAMESPACE_CLOSE


#include <igatools/utils/multi_array-inline.h>

#endif
