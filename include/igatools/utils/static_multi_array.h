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

// QualityAssurance: martinelli, 31 Jan 2014

#ifndef STATIC_MULTI_ARRAY_H_
#define STATIC_MULTI_ARRAY_H_

#include <igatools/base/config.h>
#include <igatools/utils/multi_array.h>
#include <igatools/utils/safe_stl_array.h>

IGA_NAMESPACE_OPEN


/**
 * @brief Static multi-dimensional array container, i.e. a tensor-like array
 * container of fixed dimension and fixed @p rank.
 *
 * Basically is a MultiArray in which the STL container is a SafeSTLArray.
 *
 * @see MultiArray
 *
 * @ingroup multi_array_containers
 *
 * @author M. Martinelli
 * @date 31 Jan 2014
 *
 */
template< class T, int dim, int rank >
class StaticMultiArray : public MultiArray<SafeSTLArray<T,constexpr_pow(dim,rank)>,rank>
{

public:

    /** Number of entries of type T in the StaticMultiArray */
    static constexpr Size n_entries = constexpr_pow(dim,rank);



    /** @name Constructors and destructor */
    ///@{
    /**
     * Construct a StaticMultiArray in which its elements are in undefined state.
     */
    StaticMultiArray();

    /**
     * Construct a StaticMultiArray in which its elements are set to be equal
     * to the argument value @p val.
     */
    explicit StaticMultiArray(const T &val);


    /**
     * Initializer list constructor.
     * Each element in the list is copied to the data_ vector.
     * @note The size of the list must be equal to the number of entries in the
     * StaticMultiArray object, otherwise in Debug mode an exception will be raised.
     */
    StaticMultiArray(std::initializer_list<T> list);

    /**
     * Copy constructor. It performs a deep copy of the StaticMultiArray @p m_array.
     */
    StaticMultiArray(const StaticMultiArray<T,dim,rank> &m_array) = default;


    /**
     * Move constructor.
     */
    StaticMultiArray(StaticMultiArray<T,dim,rank> &&m_array) = default;

    /**
     * Destructor. It does nothing.
     */
    ~StaticMultiArray() = default;
    ///@}


    /** @name Assignment operators */
    ///@{

    /** Copy assignment operator. Performs a deep copy of the static multi-array @p m_array. */
    StaticMultiArray<T,dim,rank> &operator=(const StaticMultiArray<T,dim,rank> &m_array) = default;

    /** Move assignment operator.*/
    StaticMultiArray<T,dim,rank> &operator=(StaticMultiArray<T,dim,rank> &&m_array) = default;
    ///@}

private:

    /** Type of the base class. */
    using base_t = MultiArray<SafeSTLArray<T,n_entries>,rank>;

};

IGA_NAMESPACE_CLOSE

#include <igatools/utils/static_multi_array-inline.h>

#endif // #ifndef STATIC_MULTI_ARRAY_H_
