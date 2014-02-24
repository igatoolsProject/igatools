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


//TODO: maybe tensorproduct data and Product array can be just one class
//TODO: it should be called ProductStaticMultiArray

#ifndef __PRODUCT_ARRAY_H_
#define __PRODUCT_ARRAY_H_

#include <igatools/base/config.h>
#include <igatools/base/logstream.h>
#include <igatools/base/tensor.h>
#include <igatools/utils/tensor_sized_container.h>
//#include <type_traits>

IGA_NAMESPACE_OPEN



/**
 * @brief Dynamic sized, tensor product type of multi dimensional array.
 *
 * Special type of <tt>rank></tt>-dimensional multidimensional array container.
 * Each entry of the array can be though as indexed by
 * rank number of indices, but each entry is the ordered
 * rank-tuple generated as the  cartesian (or direct) product of
 * @p rank one-dimensional vectors.
 *
 * To understand better let's consider the following
 * typical example: a cartesian grid in two space dimensions,
 * in this case all two dimensional vertices of the grid can be generated
 * only by storing two vectors, one for the <b>x</b> and one for the <b>y</b> direction.
 *
 * \f{eqnarray*}{
 * \mathbf{x}&=&\{0,1,4\} \\
 * \mathbf{y}&=&\{0,2,3,5\}
 * \f}
 *
 * Then the 12 2D-vertices can be generated as the cartesian product:
 * \f{eqnarray*}{
 * \mathbf{P}_0 &=& (0,0) \\
 * \mathbf{P}_1 &=& (1,0) \\
 * \mathbf{P}_2 &=& (4,0) \\
 * \mathbf{P}_3 &=& (0,2) \\
 * \vdots \\
 * \mathbf{P}_{11} &=& (4,5)
 * \f}
 *
 * A related class that is derived from this on is the TensorProductArray, when
 * the type <tt>T</tt> is a scalar type.
 *
 *
 * @todo write code of how to use it
 *
 * @ingroup multi_array_containers
 * @author Martinelli, 2012, 2013, 2014
 * @author Pauletti, 2012,2013
 */
template< class T, int rank>
class ProductArray : public TensorSizedContainer<rank>
{
public:

    /** @name Constructors */
    ///@{
    /**
     * Construct a rank-dimensional empty array.
     */
    ProductArray() ;

    // TODO (pauletti, Dec 23, 2013): Document this
    ProductArray(std::initializer_list<std::initializer_list<T>> list);

    /**
     * Construct a rank-dimensional ProductArray where the
     * the i-th direction is initialized to be
     * of size size[i], calling the default constructor of T for
     * each entry of the vectors.
     */
    explicit ProductArray(const TensorSize<rank> size);

    /**
     * Constructor. Same as ProductArray(const array< int, rank > size)
     * but with all direction sizes equal to size.
     */
    explicit ProductArray(const Size size);


    /**
     * Constructor. Construct a rank-dimensional ProductArray where the
     * the i-th direction is initialized to be equal to @p data_directions[i]
     */
    explicit ProductArray(const std::array<std::vector<T>,rank> &data_directions) ;


    /**
     * Destructor.
     */
    ~ProductArray() = default;

    /**
     * Copy constructor.
     */
    ProductArray(const ProductArray<T,rank> &product_array) = default;

    /**
     * Move constructor.
     */
    ProductArray(ProductArray<T,rank> &&product_array) = default;

    ///@}


    /** @name Assignment operators */
    ///@{
    /**
     * Copy assignment operator.
     * @note Use with care as it may be an expensive operation.
     */
    ProductArray<T,rank> &
    operator=(const ProductArray<T,rank> &product_array) = default;


    /** Move assignment operator. */
    ProductArray<T,rank> &
    operator=(ProductArray<T,rank> &&product_array) = default;
    ///@}

private:
    /**
     * Type traits magic to return different types where T is an arithmetic type
     * (floating point or integer)
     * or not, used in the cartesian_product() functions.
     */
    using point_t = Conditional<std::is_arithmetic<T>::value,
          Conditional<std::is_floating_point<T>::value,
          Point<rank>,
          TensorIndex<rank>>,
          std::array<T,rank> >;

    /**
     * sub-ProductArray of the ProductArray.
     * @todo make it a template function
     * @warning There must be a bug in the compiler as it cannot put >0 instead !=0
     */
    using sub_product_t = Conditional< (rank > 0),
          ProductArray<T,rank-1>,
          ProductArray<T,0> >;

public:


    /** @name Functions for accessing (read/write/copy) the internal data */
    ///@{
    /**
     * Read/write access function to the <tt>j</tt>-th entry along the <tt>i</tt>-th direction.
     * @note In Debug mode there is a check on the indices <tt>(i,j)</tt> used as input parameters.
     */
    T &entry(const int i, const int j) ;

    /**
     * Copy input @p data to the data relative to the @p i-th direction.
     * @note The ProductArray object will be internally resized (if needed)
     * in order to contains all the entries in the input @p data.
     */
    void copy_data_direction(const int i, const std::vector<T> &data) ;

    /**
     * Get a const-reference to the vector data of the <tt>i</tt>-th direction.
     */
    const std::vector<T> &get_data_direction(const int i) const ;


    /**
     * Returns a sub ProductArray of the ProductArray.
     * @todo document more.
     */
    sub_product_t get_sub_product(const TensorIndex<rank-1> &index) const;

    /**
     * Given the product index, it returns the corresponding
     * cartesian product.
     * For example if index=[3,1], cartesian product returns
     * and array [data[0][3], data[1][1]].
     * In the case T is a floating point it returns a Point type.
     */
    point_t cartesian_product(const TensorIndex<rank> &index) const;

    /**
     * Returns a flat vector of the cartesian product of the
     * product array.
     * For example if the ProductArray is
     * \f[
     * \mathbf{x}=\{1,2\} \quad , \quad
     * \mathbf{y}=\{4,3\}
     * \f]
     * then this function returns
     * \f[
     * \{(1,4), (1,3), (2,4), (2, 3)\}
     * \f]
     */
    std::vector< point_t > get_flat_cartesian_product() const;

    ///@}


    /**
     * Resize the vectors containing the data, where
     * size[i] is the number of element for the i-th direction.
     */
    void resize(const TensorSize<rank> &size);

protected:
    /**
     * This member contain the data for each coordinate direction.
     * data_[i][j] refers to the j-th data element  along the
     * i-th coordinate direction.
     */
    std::array<std::vector<T>,rank> data_ ;

};


/**
 * Output operator
 * @relates ProductArray
 */
template< class T, int rank>
LogStream &operator<<(LogStream &out, const ProductArray< T, rank> &data) ;


IGA_NAMESPACE_CLOSE

// If we are in debug mode we do not inline to gain some compilation speed,
// but we do in Release mode
#ifdef NDEBUG
#include <igatools/utils/product_array-inline.h>
#endif


#endif
