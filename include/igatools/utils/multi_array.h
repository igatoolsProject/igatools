//-+--------------------------------------------------------------------
// Igatools a general purpose Isogeometric analysis library.
// Copyright (C) 2012-2014  by the igatools authors (see authors.txt).
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
#include <igatools/utils/tensor_sized_container.h>
#include <igatools/utils/multi_array_iterator.h>



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
 * accepting both multi-indices (trhough its representation via the class TensorIndex)
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
 * @ingroup multi_array_containers
 *
 * @author M. Martinelli
 * @date 31 Jan 2014
 */
template<class STLContainer, int rank>
class MultiArray : public TensorSizedContainer<rank>
{
public:
    /** Type of the entries stored in the STL container. */
    using Entry = typename STLContainer::value_type;


    /** @name Constructors and destructor */
    ///@{
    /**
     * Construct an empty multiarray.
     */
    MultiArray() = default;

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
     *  Flat index access operator (non-const version).
     */
    Entry &operator()(const Index i);

    /**
     *  Flat index access operator (const version).
     */
    const Entry &operator()(const Index i) const;


    /**
     *  Tensor index access operator (non-const version).
     */
    Entry &operator()(const TensorIndex<rank> &i);

    /**
     *  Tensor index access operator (const version).
     */
    const Entry &operator()(const TensorIndex<rank> &i) const;

    /** Return the entries of the multiarray as unidimensional std::vector. */
    const STLContainer &get_data() const;
    ///@}


    /** @name Dealing with the iterators */
    ///@{

    /** Type of the const iterator. */
    using const_iterator = MultiArrayIterator<const MultiArray<STLContainer,rank>>;

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


    /** @name Functions to easily fill the multiarray */
    ///@{
    /**
     * Fills the multiarray with an arithmetic progression starting with the value @p init_value.
     */
    void fill_progression(const Entry &init_value = {});


    /** Fills the multiarray copying in each entry the content of @p value. */
    void fill(const Entry &value);
    ///@}

protected:

    /**
     * Data of type Entry stored in a STL container.
     * @note This member variable is protected because should be resized by DynamicMultiArray.
     */
    STLContainer data_;
};



IGA_NAMESPACE_CLOSE


#include <igatools/utils/multi_array-inline.h>

#endif
