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

// QualityAssurance: martinelli, 21 Jan 2014

#ifndef TENSOR_INDEX_H_
#define TENSOR_INDEX_H_

#include <igatools/base/config.h>
#include <igatools/base/logstream.h>
#include <igatools/utils/vector.h>
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
 * This class inherits publicly from std::array, but reimplements the access operator[]
 * for bounds checking.
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

    /**
     * Default constructor. Initializes all the direction indices to
     * @p value. If no @p value is provided, the initialization is made with zeros.
     */
    explicit TensorIndex(const Size value = 0) noexcept ;

    /** Constructor using an std::array. */
    explicit TensorIndex(const std::array<Index,rank> &arr) noexcept;

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

    /** @name Access (read/write) operators */
    ///@{
    /**
     * Read/write operator. Returns the reference to the i-th size.
     * @note In Debug mode the index @p i is checked in order to be
     * in the bounds of the array containing the different direction sizes.
     */
    Index &operator[](const Index i);

    /**
     * Read operator. Returns the const-reference to the i-th size.
     * @note In Debug mode the index @p i is checked in order to be
     * in the bounds of the array containing the different direction sizes.
     */
    const Index &operator[](const Index i) const;
    ///@}


    template<int k>
    TensorIndex<k> get_sub_tensor(const TensorIndex<k> &index) const
    {
        TensorIndex<k> res;
        int j = 0;
        for (auto i : index)
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

template <int rank>
bool
operator<(const TensorIndex<rank> &index_a,const TensorIndex<rank> &index_b)
{
    for (int j=0; j<rank; ++j)
        if (index_a[j] >= index_b[j])
            return false;
    return true;
}
/**
 * Output operator for TensorIndex.
 *
 * @relates TensorIndex
*/
template <int rank>
LogStream &
operator<<(LogStream &out, const TensorIndex<rank> &tensor_index) ;


/**
 * Generates a vector with the tensor indices of the given
 * rectangular range.
 *
 */
template<int k>
vector<TensorIndex<k>> tensor_range(TensorIndex<k> first, TensorIndex<k> last)
{
    Assert(first < last, ExcMessage("first not smaller than last"));
    vector<TensorIndex<k>> result;
    TensorIndex<k-1> ind(arr::sequence<k-1>());
    auto vec = tensor_range<k-1>(first.get_sub_tensor(ind), last.get_sub_tensor(ind));

    for (int i=first[k-1]; i<last[k-1]; ++i)
    {
        for (auto &t_k_1 : vec)
        {
            TensorIndex<k> t_k;
            for (int j=0; j<k-1; ++j)
                t_k[j] = t_k_1[j];
            t_k[k-1] = i;
            result.emplace_back(t_k);
        }
    }
    return result;
}

template<>
vector<TensorIndex<1>> tensor_range(TensorIndex<1> first, TensorIndex<1> last);

template<>
vector<TensorIndex<0>> tensor_range(TensorIndex<0> first, TensorIndex<0> last);

IGA_NAMESPACE_CLOSE

#ifdef NDEBUG
#include <igatools/utils/tensor_index-inline.h>
#endif

#endif // TENSOR_INDEX_H_



