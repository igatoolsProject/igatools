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


#ifndef MULTI_ARRAY_UTILS_H_
#define MULTI_ARRAY_UTILS_H_

#include <igatools/base/config.h>
#include <igatools/utils/tensor_index.h>
#include <igatools/utils/tensor_size.h>
#include <igatools/utils/safe_stl_vector.h>

IGA_NAMESPACE_OPEN

/**
 * @brief Utils for managing multi-dimensional array indices.
 *
 * This class provides the tools to handle in a unified way
 * the tensor and flat index convention used in all the multidimensional array
 * containers provided by igatools.
 *
 * \section tensor_flat_index Convention for flat and tensor indices in multi-arrays
 * The library uses multi-dimensional arrays in many cases.
 * Sometimes is convenient to access its entries using the tensor index and
 * some other times using a flat index.
 * In this circumstances it is necessary to define a convention on the
 * conversion for flat to tensor index and viceversa.
 * In igatools we hide (encapsulate) the arbitrary choice that we follow
 * in the class MultiArrayUtils that provides the necessary method
 * to perform this transformations without knowing the internal details.
 *
 * For the curious or the developer interested in working on the class
 * the convention that we follow is:
 * "that the first index moves faster".
 * In this way each tensor index has a different weight for a
 * given tensor data, the weight is how many flatten units
 * a unit at a given index worth.
 * For example....(todo)
 * From here we see that the weight is computed from the size of
 * each direction vector
 * { 1, size[0], size[0]*size[1], size[0]*size[1]*size[2], ...}
 *
 * @ingroup multi_array_containers
 *
 * @author S.Pauletti
 */
template <int rank>
class MultiArrayUtils
{
public:
    /**
     * Given the tensor data component weights, and a
     * flat index it return the corresponding tensor index (see details in
     * the general documentation of CartesianProductArray.
     */
    static TensorIndex<rank>
    flat_to_tensor_index(const Index flat_index,
                         const TensorIndex<rank> &weight) noexcept ;

    /**
     * Given the tensor data component weights, and a
     * tensor index it return the corresponding linear index (see details in
     * the general documentation of CartesianProductArray.
     */
    static Index
    tensor_to_flat_index(const TensorIndex<rank> &tensor_index,
                         const TensorIndex<rank> &weight) noexcept ;

    /**
     * Given the size of the vectors in a tensor data, returns
     * the weight of each tensor index.
     * (see details in the general documentation of CartesianProductArray).
     */
    static TensorIndex<rank>
    compute_weight(const TensorSize<rank> &size) noexcept ;


    /**
     * Computes the size of a given tensor index
     * @todo remove this function because a TensorIndex is not a TensorSize
     */
    static Size
    size(const TensorIndex<rank> &extend) noexcept ;

    /**
     * This an auxiliary function needed in local to global,
     * what it does is given a raw tensor data of int type, a flat index and weight
     * it returns an array with the corresponding data in each entry
     */
    static TensorIndex<rank>
    get_tensor_entry(const SafeSTLArray< SafeSTLVector<Index>, rank> &data,
                     const Index flat_index, const TensorIndex<rank> &weight) noexcept ;
};




IGA_NAMESPACE_CLOSE

#include <igatools/utils/multi_array_utils-inline.h>


#endif // #ifndef MULTI_ARRAY_UTILS_H_
