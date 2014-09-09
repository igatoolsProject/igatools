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



#ifndef CARTESIAN_PRODUCT_ARRAY_H_
#define CARTESIAN_PRODUCT_ARRAY_H_

#include <igatools/base/config.h>
#include <igatools/base/logstream.h>
#include <igatools/base/tensor.h>
#include <igatools/utils/vector.h>
#include <igatools/utils/tensor_sized_container.h>
#include <igatools/utils/value_vector.h>


IGA_NAMESPACE_OPEN


/**
 * @brief Dynamic sized, tensor product type of multi dimensional array.
 *
 * Special type of <tt>rank</tt>-dimensional multidimensional array container.
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
class CartesianProductArray: public TensorSizedContainer<rank>
{
public:
    using EntryType = T;
    /**
     * Type for the <tt>rank-1</tt> CartesianProductArray
     */
    using SubProduct = Conditional< (rank > 0),
          CartesianProductArray<T,rank-1>,
          CartesianProductArray<T,0> >;



    /** @name Constructors */
    ///@{
    /**
     * Construct a rank-dimensional empty array.
     */
    CartesianProductArray() ;

    // TODO (pauletti, Dec 23, 2013): Document this
    CartesianProductArray(std::initializer_list<std::initializer_list<T>> list);

    /**
     * Construct a rank-dimensional CartesianProductArray where the
     * the i-th direction is initialized to be
     * of size size[i], calling the default constructor of T for
     * each entry of the vectors.
     */
    explicit CartesianProductArray(const TensorSize<rank> size);

    CartesianProductArray(const TensorSize<rank> size, T& val);

    /**
     * Constructor. Same as CartesianProductArray(const array< int, rank > size)
     * but with all direction sizes equal to @p size.
     */
    explicit CartesianProductArray(const Size size);


    /**
     * Constructor. Construct a rank-dimensional CartesianProductArray where the
     * the i-th direction is initialized to be equal to @p data_directions[i]
     */
    explicit CartesianProductArray(const std::array<vector<T>,rank> &data_directions) ;


    /**
     * Destructor.
     */
    ~CartesianProductArray() = default;

    /**
     * Copy constructor.
     */
    CartesianProductArray(const CartesianProductArray<T,rank> &product_array) = default;

    /**
     * Move constructor.
     */
    CartesianProductArray(CartesianProductArray<T,rank> &&product_array) = default;

    ///@}


    /** @name Assignment operators */
    ///@{
    /**
     * Copy assignment operator.
     * @note Use with care as it may be an expensive operation.
     */
    CartesianProductArray<T,rank> &
    operator=(const CartesianProductArray<T,rank> &product_array) = default;


    /** Move assignment operator. */
    CartesianProductArray<T,rank> &
    operator=(CartesianProductArray<T,rank> &&product_array) = default;
    ///@}

private:
    /**
     * Type traits magic to return different types where T is an arithmetic type
     * (floating point or integer)
     * or not, used in the cartesian_product() functions.
     */
    using point_t = Conditional<std::is_arithmetic<T>::value,
          Conditional<std::is_floating_point<T>::value,
          Points<rank>,
          TensorIndex<rank>>,
          std::array<T,rank> >;


public:


    /** @name Functions for accessing (read/write/copy) the internal data */
    ///@{
    /**
     * Read/write access function to the <tt>j</tt>-th entry along the <tt>i</tt>-th direction.
     * @note In Debug mode there is a check on the indices <tt>(i,j)</tt> used as input parameters.
     */
    T &entry(const int i, const int j) ;

    T const &entry(const int i, const int j) const;

    /**
     * Copy input @p data to the data relative to the @p i-th direction.
     * @note The CartesianProductArray object will be internally resized (if needed)
     * in order to contains all the entries in the input @p data.
     */
    void copy_data_direction(const int i, const vector<T> &data) ;

    /**
     * Get a const-reference to the vector data of the <tt>i</tt>-th direction.
     */
    const vector<T> &get_data_direction(const int i) const ;

    ///@}

    /** @name Functions returning rank-1 and rank+1 objects */
    ///@{
    using SubProductTensorIndex = Conditional<(rank>0),TensorIndex<rank-1>,TensorIndex<0>>;

    /**
     * Returns a rank-1 CartesianProductArray built
     * copying some part of the the data from the calling object.
     *
     * The data to be copied is selected using the @p index argument, where
     * <tt>index[j]</tt> means that the data along <tt>j</tt>-th direction of the new object
     * is the copy of the data along the <tt>index[j]</tt>-th direction of the old object.
     * @code
       //example
       CartesianProductArray<T,rank> old_obj;
       CartesianProductArray<T,rank-1> new_obj;
       //
       new_obj.data_[j] = old_obj.data_[index[j]];
       @endcode
     *
     */
    SubProduct get_sub_product(const SubProductTensorIndex &index) const;



    ///@}

    /** @name Functions returning the cartesian product entries */
    ///@{

    /**
     * Given the product index, it returns the corresponding
     * cartesian product.
     * For example if index=[3,1], cartesian product returns
     * and array [data[0][3], data[1][1]].
     * In the case T is a floating point it returns a Point type.
     */
    point_t cartesian_product(const TensorIndex<rank> &index) const;

    std::array<T const &, rank> cartesian_product_adresss(const TensorIndex<rank> &index) const;

    /**
     * Returns a flat vector of the cartesian product of the
     * product array.
     * For example if the CartesianProductArray is
     * \f[
     * \mathbf{x}=\{1,2\} \quad , \quad
     * \mathbf{y}=\{4,3\}
     * \f]
     * then this function returns
     * \f[
     * \{(1,4), (1,3), (2,4), (2, 3)\}
     * \f]
     */
    Conditional<std::is_floating_point<T>::value,ValueVector<point_t>,vector<point_t> >
    get_flat_cartesian_product() const;
    ///@}

    /**
     * Resize the vectors containing the data, where
     * size[i] is the number of element for the i-th direction.
     */
    void resize(const TensorSize<rank> &size);

    /**
     * Prints some internal informations. Mainly used for testing and debugging purposes.
     */
    void print_info(LogStream &out) const;

protected:
    /**
     * This member contain the data for each coordinate direction.
     * data_[i][j] refers to the j-th data element  along the
     * i-th coordinate direction.
     */
    std::array<vector<T>,rank> data_ ;
};


/**
 * Returns a CartesianProductArray of one higher rank built from the insertion
 * of a given @p new_vector at in the CartesianProductArray @p orig at the given direction @p index.
 */
template <class T, int rank>
CartesianProductArray<T, rank+1>
insert(const CartesianProductArray<T, rank> &orig,
       const int index,
       const vector<T> &new_vector);


IGA_NAMESPACE_CLOSE

// If we are in debug mode we do not inline to gain some compilation speed,
// but we do in Release mode
#ifdef NDEBUG
#include <igatools/utils/cartesian_product_array-inline.h>
#endif


#endif // #ifndef CARTESIAN_PRODUCT_ARRAY_H_
