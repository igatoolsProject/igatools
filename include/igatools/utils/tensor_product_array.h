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



#ifndef TENSOR_PRODUCT_ARRAY_H_
#define TENSOR_PRODUCT_ARRAY_H_

#include <igatools/utils/cartesian_product_array.h>

IGA_NAMESPACE_OPEN




/**
 * @brief Specialization of CartesianProductArray when the entries are a product
 * of Reals instead of a tuple of Real.
 *
 * This is a container for a special type of @p rank dimensional multiarray (or
 * rank rank tensor). More specifically, the type of multiarray that can be
 * generated as a tensor product of @p rank one-dimensional
 * vectors (that we call directions).
 *
 * @ingroup multi_array_containers
 * @author Martinelli, 2012, 2013, 2014
 *
 */
template<int rank>
class TensorProductArray :
    public CartesianProductArray<Real,rank>
{
public :

    /** @name Constructors. */
    ///@{
    TensorProductArray() = default;

    /**
     * The i-th vector of the the array is initialized to be
     * of size size[i], calling the default constructor of T for
     * each entry of the vectors.
     */
    explicit TensorProductArray(const TensorSize<rank> &size);

    /**
     * Same as CartesianProductArray(const array< int, rank > size)
     * but with all direction sizes equal to size.
     */
    explicit TensorProductArray(const Size size);


    /**
     * Copy constructor.
     */
    TensorProductArray(const TensorProductArray<rank> &data) = default;

    /**
     * Move constructor.
     */
    TensorProductArray(TensorProductArray<rank> &&data) = default;

    /**
     * Constructor from the CartesianProductArray object @p data.
     *
     * @note The use of this constructor is safe because TensorProductArray<rank> and
     * CartesianProductArray<Real,rank> differs only for the members functions and
     * not for the member variables.
     */
    TensorProductArray(const CartesianProductArray<Real,rank> &data);

    /**
     * Destructor
     */
    ~TensorProductArray() = default;
    ///@}


    /** @name Assignment operators */
    ///@{
    /**
     * Copy assignment operator.
     */
    TensorProductArray<rank> &operator=(const TensorProductArray<rank> &data) = default;

    /**
     * Move assignment operator.
     */
    TensorProductArray<rank> &operator=(TensorProductArray<rank> &&data) = default;
    ///@}



    /** @name Functions for performing dilation and translation of the values in the container */
    ///@{
    /** Dilation followed by a translation of the TensorProductArray data. */
    void dilate_translate(const std::array<Real,rank> &dilate,
                          const Point<rank> &translate) ;

    /** Dilation of the TensorProductArray data.*/
    void dilate(const std::array<Real,rank> &dilate) ;


    /** Translation of the TensorProductArray data.*/
    void translate(const Point<rank> &translate) ;
    ///@}

    /** @name Functions returning the tensor product entries */
    ///@{
    /**
     * Returns the entry identified by the tensor index given by @p tensor_id.
     */
    Real tensor_product(const TensorIndex<rank> &tensor_id) const;

    /**
     * Returns a flat vector with the component of the tensor generated
     * through the tensor product operation.
     */
    std::vector<Real> get_flat_tensor_product() const;
    ///@}


public :
};







IGA_NAMESPACE_CLOSE

// If we are in debug mode we do not inline to gain some compilation speed,
// but we do in Release mode
#ifdef NDEBUG
#include <igatools/utils/tensor_product_array-inline.h>
#endif


#endif // #ifndef TENSOR_PRODUCT_ARRAY_H_
