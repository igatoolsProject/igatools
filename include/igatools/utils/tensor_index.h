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

// QualityAssurance: martinelli, 21 Jan 2014

#ifndef TENSOR_INDEX_H_
#define TENSOR_INDEX_H_

#include <igatools/base/config.h>
#include <igatools/base/logstream.h>




IGA_NAMESPACE_OPEN

/**
 * @brief Type for a tensor index.
 *
 * It is a list of rank number of indices.
 *
 * This class makes possible the
 * rank independent treatment of tensor type containers.
 *
 * @author M.Martinelli
 * @date 20 Jan 2014
 */
template <int rank>
class TensorIndex : public std::array<Index, rank>
{
public:

    /** @name Constructors */
    ///@{

    /** Default constructor. Initializes all the direction indices to zero. */
    TensorIndex(const Index value = 0) noexcept ;

    /** Constructor using an std::array. */
    TensorIndex(const std::array<Index,rank> &arr) noexcept;

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
     * @name Increment/decrement operators.
     */
    ///@{

    /** Increment the tensor index by the values @p idx.*/
    TensorIndex<rank> &operator +=(const TensorIndex<rank> &idx) noexcept ;

    /** Increment the indices in all directions by the value @p j. */
    TensorIndex<rank> &operator +=(const Index j) noexcept ;

    /** Decrement the indices in all directions by the value @p j. */
    TensorIndex<rank> &operator -=(const Index j) noexcept ;
    ///@}
};



/**
 * Output operator for TensorIndex.
 *
 * @relates TensorIndex
*/
template <int rank>
LogStream &
operator<<(LogStream &out, const TensorIndex<rank> &tensor_index) ;

IGA_NAMESPACE_CLOSE

#ifdef NDEBUG
#include <igatools/utils/tensor_index-inline.h>
#endif

#endif // TENSOR_INDEX_H_



