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

#ifndef _TENSORS_H
#define _TENSORS_H

#include <igatools/base/config.h>
#include <igatools/base/exceptions.h>
#include <igatools/base/logstream.h>
#include <igatools/utils/tensor_index.h>

#include <type_traits>
#include <cstring>
#include <cmath>

IGA_NAMESPACE_OPEN

//TODO (Nov 9, 2013, antolin): move implementation of functions to tensor.cpp or
//                              in inline to tensor-inline.h.
// TODO (pauletti, Feb 20, 2014): some members are missing documentation

template<int dim_, int rank_, class tensor_type, class value_type>
class Tensor;

/**
 * This namespace collects structures to represent the different types of
 * tensor.
 * For more details of its meaning see the documentation on the
 * Tensor class.
 * @relates Tensor
 */
namespace tensor
{

struct covariant;

struct contravariant
{
    using co_type = covariant;
};

struct covariant
{
    using co_type = contravariant;
};

struct raw
{
    using co_type = raw;
};

} // end namespace tensor


/**
 * The double class we need for tensors.
 * @note From the theoretical point of view this class should not
 * be necessary, provided that we can inherit from double.
 * Unfortunately C++ does not allow to inherit from built-in types so
 * we have to wrap it using this class.
 * All functions are inlined so in principle there should be no difference
 * in performance.
 */
class Tdouble
{
public:
    static const bool is_tensor = false;

    /** Dimension of the vector space */
    static const int dim   = 1; //Could actually be any number

    /** Rank of the tensor */
    static const int rank  = 0;

    /** Flat size of the tensor, i.e.
     * total number of components of type value_type */
    static const int size = 1; //iga::constexpr_pow(dim_, rank_);

    using self_t = Tdouble;
    using co_tensor_t = self_t;
    using value_t = Real;

    /** @name Constructors */
    ///@{
    /** Default constructor */
    Tdouble(const Real val = 0.);

    /**
     * Advance constructor to optimization when zero
     * init is not necessary
     */
    Tdouble(const bool non_init);

    /** Copy constructor */
    Tdouble(const Tdouble &td) = default;

    /** Move constructor */
    Tdouble(Tdouble &&td) = default;

    /** Destructor */
    ~Tdouble() = default;
    ///@}

    /**
     * @name Assignment operators
     */
    ///@{

    /**
     * Copy assignment operator.
     */
    Tdouble &operator=(const Tdouble &td) = default;


    /**
     * Move assignment operator.
     */
    Tdouble &operator=(Tdouble &&td) = default;


    /**
     * Assignment operator from a Real @p val1.
     */
    Tdouble &operator=(const value_t &val1);
    ///@}


    // TODO (pauletti, Feb 25, 2014): document and implementation in other file
    operator value_t &() noexcept
    {
        return val_;
    }

    // TODO (pauletti, Feb 25, 2014): document and implementation in other file
    operator value_t const &() const noexcept
    {
        return val_;
    }

    /**
     * Read write access operator
     */
    value_t &operator[](const int i) noexcept;


    /**
     * Read access operator
     */
    const value_t &operator[](const int i) const noexcept;

    /**
     * Read write access operator using a tensor index
     */
    value_t &operator()(const TensorIndex<0> &i) noexcept;

    /**
     * Read access operator using a tensor index
     */
    const value_t &operator()(const TensorIndex<0> &i) const noexcept;

    /**
     * Read write access operator using the flat index
     */
    value_t &operator()(const int i) noexcept;

    /**
     * Read access operator using the flat index
     */
    const value_t &operator()(const int i) const noexcept;

    Tdouble &operator+=(const Real td) noexcept;

    Tdouble &operator-=(const Real td) noexcept;

    Tdouble &operator*=(const Real td) noexcept;

    Tdouble &operator/=(const Real td) noexcept;

    Real norm() const noexcept;

    Real norm_square() const noexcept;

    TensorIndex<0>
    flat_to_tensor_index(const int flat_index) const noexcept;

    /**
     * The last index moves faster.
     */
    int tensor_to_flat_index(const TensorIndex<0> &tensor_index) const noexcept;


private:
    value_t val_;

};



/**
 * @name Tensor type aliases
 * @relates Tensor
 */
///@{

/**
 * This alias represents the co-tensor of the tensor T
 * @relates Tensor
 */
template<typename T>
using CoTensor =
    Conditional<
    T::is_tensor,
    Tensor<T::dim, T::rank, typename T::tensor_t::co_type, typename T::value_t::co_tensor_t>,
    Tdouble>;




/**
 * This alias represents the transpose of the tensor T
 * @relates Tensor
 */
template<typename T>
using Transpose =
    Conditional<
    T::is_tensor,
    Conditional<
    T::value_t::rank==0,
    T,
    Tensor<T::value_t::dim,T::value_t::rank,typename T::value_t::tensor_t,
    Tensor<T::dim,T::rank,typename T::tensor_t,typename T::value_t::value_t> >
    >,
    Tdouble >;



/**
 * This alias represents the sub-tensor part of tensor of type T,
 * i.e. a tensor with same dim and tensor_type of T
 * but with rank equal to T::rank - 1.
 * @note If T has rank <= 1 then the sub-tensor type is T::value_t
 * @relates Tensor
 */
template<typename T>
using SubTensor =
    Conditional<
    (T::rank<=1),
    typename T::value_t,
    Tensor<T::dim,T::rank-1,typename T::tensor_t,typename T::value_t>>;


/**
 * This alias represents the tensor A(x), i.e. the action of the
 * tensor A (of type T) on x
 * @relates Tensor
 */
template<typename T>
using ActionTensor = Conditional<
                     T::is_tensor,
                     SubTensor<T>,
                     Real>;
///@}

/**
 *  @brief Mathematical tensor of given rank and type.
 *
 *  The Tensor class is used to represent a mathematical tensor of
 *  given rank and type.
 *  This class is a core computational support for the library.
 *  It is used for storing and computations with values and derivatives
 *  of all orders of basis function, mapping and general functions.
 *
 *
 *  Mathematically, given a \f$ d \f$ -dimensional vector space \f$ V \f$,
 *  an (m,n)-tensor over  \f$ V \f$ is a scalar multilinear function
 *  of @p m one-forms and @p n vectors, i.e.
 *  \f[ A: \underbrace{V^* \times \dots \times V^*}_{m \mbox{ times }} \times
 *         \underbrace{V   \times \dots \times V  }_{n \mbox{ times }} \to \mathbb R \f]
 *  where \f$ V^* \f$ denotes the dual of \f$ V \f$.
 * It is usually refer to as an @p m times @e contravariant, @p n times @e covariant tensor.
 * The \e rank of \f$ A \f$ is m+n.
 *
 * For example, the bilinear form \f$ T: V \times V^* \to \mathbb R \f$
 * is a (1,1)-tensor or rank 1 covariant, rank 1 contravariant tensor.
 * In the language of the Tensor class it is
 * represented with the following code
 * \code
 * Tensor<dim, 1, tensor::covariant, Tensor<dim, 1, tensor::contravariant> > T;
 * \endcode
 *
 * Similarly, the bilinear function \f$ S: V^* \times V \to \mathbb R \f$ would be
 * \code
 * Tensor<dim, 1, tensor::contravariant, Tensor<dim, 1, tensor::covariant> > S;
 * \endcode
 *
 * Once a basis for the vectors is selected and the corresponding dual basis for
 * the covectors
 * considered, we can talk about the components of a tensor.
 * Assume \f$ \mathbf{e}_1, \dots, \mathbf{e}_{dim} \f$ and \f$ \mathbf{e}^1, \dots, \mathbf{e}^{dim} \f$ are vector and co-vector
 * bases respectively then the tensor can be written in terms of its components \f$ T^{i}_{\phantom{i} j} \f$ as :
 * \f[
 * T = \sum_{i=1}^{d} \sum_{j=1}^{d} T^{i}_{\phantom{i} j} \mathbf{e}_i \otimes \mathbf{e}^j
 * \f]
 * where we have used the upper index for the contravariant components and the lower index for the covariant components.
 *
 * In this case the bilinear function \f$ T(\mathbf{e}_i,\mathbf{e}^j) \f$ is given
 * by T[i][j] and the contravariant tensor \f$ T(\mathbf{e}_i) : V^* \to \mathbb R\f$
 * is given by T[i].
 *
 * Now if we consider two vector spaces \f$ V \f$ and \f$ W \f$ of dimension m and
 * n respectively, the bilinear form \f$ R: V \times W^* \to \mathbb R \f$,
 * which is nothing else than a linear transformation
 * from \f$ V \f$ into \f$ W \f$, is in the Tensor class
 * \code
 * Tensor<m, 1, tensor::covariant, Tensor<n, 1, tensor::contravariant> > R;
 * \endcode
 *
 *
 * \section deriv Using tensors for derivatives
 * The main use of the Tensor class in igatools are derivatives of any
 * order (even the zero order).
 *
 * Consider the function \f$ F: V \to W \f$, then the derivative of
 * \f$ F \f$ at \f$ p \f$ is a linear function from \f$ V \f$ to \f$ W \f$, i.e.  \f$ DF(p) \in L(V,W) \f$
 * or \f$ DF(p): V \times W^* \to \mathbb R \f$,
 * and similarly \f$ D^2 F(p): V \times V \times W^* \to \mathbb R \f$.
 *
 * \code
 * Tensor<m, 1, tensor::covariant, Tensor<n, 1, tensor::contravariant > > DF;
 * Tensor<m, 2, tensor::covariant, Tensor<n, 1, tensor::contravariant > > D2F;
 * Tensor<m, N, tensor::covariant, Tensor<n, 1, tensor::contravariant > > DNF;
 * \endcode
 *
 * \section pathological Pathological cases
 * - dim_ == 0 is an empty object
 * - rank_ == 0 is a "scalar" or value_type
 * \section rational Rational for the Tensor class design
 * - Tensors are defined recursively on the rank_
 * - The memory required by a Tensor object is internally statically
 *   allocated using the informations passed by the template parameters
 * - A 0-rank_ tensor is a value_type for all practical purposes
 * - If an operation is applied to a rank_ 0 tensor
 *   it should be directly transfered to value_type.
 *
 * Access operators:
 * The tensor class provides two three different types of access operators
 * - sub tensor access []
 * - entry flat index access ()
 * - entry tensor index access ()
 *
 * @note Usig metaprogramming techniques we have implemented the following differences
 * between the general tensor T and a rank 1 tensor:
 * - if rank != 1, SubTensor<Tensor<dim,rank,tensor_type,value_type>> is Tensor<dim,rank-1,tensor_type,value_type>
 * - if rank == 1, SubTensor<Tensor<dim,rank,tensor_type,value_type>> is value_type
 *
 *
 * @author Martinelli 2012, 2013
 * @author Cavallini 2012
 * @author Pauletti 2012, 2013
 *
 */
template<int dim_, int rank_, class tensor_type, class value_type>
class Tensor
{
public:
    // TODO (pauletti, Feb 25, 2014): Document this
    static const bool is_tensor = true;

    /** Dimension of the vector space */
    static const int dim   = dim_;

    /** Rank of the tensor */
    static const int rank  = rank_;

    /** Flat size of the tensor, i.e. total number of components of type value_type */
    static const int size = iga::constexpr_pow(dim_, rank_);


    using tensor_t = tensor_type;

    using value_t = value_type;

    using self_t = Tensor<dim,rank,tensor_t,value_t>;

    using co_tensor_t = CoTensor<self_t>;

    /**
     * @name Constructors
     */
    ///@{

    /**
     * Default constructor. Sets all entries to 0.
     */
    Tensor() = default;

    /** Copy constructor */
    Tensor(const Tensor<dim_, rank_, tensor_type, value_type> &tensor) = default;

    /** Move constructor */
    Tensor(Tensor<dim_, rank_, tensor_type, value_type> &&tensor) = default;

    /** Constructor from an initializer-list of value_type */
    Tensor(std::initializer_list<value_type> list);


    /** Destructor */
    ~Tensor() = default;
    ///@}


    /**
     * @name Assignment operators
     */
    ///@{

    /** Copy assignment operator */
    Tensor<dim_, rank_, tensor_type, value_type> &operator=(
        const Tensor<dim_, rank_, tensor_type, value_type> &tensor) = default;



    /** Move assignment operator */
    Tensor<dim_, rank_, tensor_type, value_type> &operator=(
        Tensor<dim_, rank_, tensor_type, value_type> &&tensor) = default;


    /** Initializer-list assignment */
    Tensor<dim_, rank_, tensor_type, value_type>
    &operator=(std::initializer_list<value_type>);


    /**
     * Assignment operator using a tensor entry value.
     * After the use of this operator the tensor entries will be set to the value specified by
     * @p entry_val.
     */
    Tensor<dim_, rank_, tensor_type, value_type> &operator=(
        const value_type &entry_val);

    /**
     * Assigning using a Real. For safety of meaning only assigning 0 is allowed.
     */
    Tensor<dim_, rank_, tensor_type, value_type> &operator=(const Real value);
    ///@}


    /**
     * @name Flat- and Tensor-index access operators
     */
    ///@{

    /** Read write access operator using a tensor index */
    value_type &operator()(const TensorIndex<rank_>  &i);

    /** Read access operator using a tensor index */
    const value_type &operator()(const TensorIndex<rank_>  &i) const;

    /** Read write access operator using the flat index */
    value_type &operator()(const int i);

    /**
     * Read access operator using the flat index
     */
    const value_type &operator()(const int i) const;

    ///@}

    /**
     * @name Sub-tensor access operators
     */
    ///@{
    /**
     * Read write access operator
     */
    SubTensor<self_t> &operator[](const int  i);

    /**
     * Read access operator
     */
    const SubTensor<self_t> &operator[](const int i) const;
    ///@}


    /**
     * @name Basic mathematical operations
     */
    ///@{
    /**
     * Increment the tensor
     */
    Tensor<dim_, rank_, tensor_type, value_type> &
    operator+=(const Tensor<dim_, rank_, tensor_type, value_type> &tensor);

    /**
     * Decrement the tensor
     */
    Tensor<dim_, rank_, tensor_type, value_type> &
    operator-=(const Tensor<dim_, rank_, tensor_type, value_type> &tensor);

    /**
     * Multiply the tensor by a scalar value
     */
    Tensor<dim_, rank_, tensor_type, value_type> &
    operator*=(const Real value);

    /**
     * Divide tensor by a scalar value
     */
    Tensor<dim_, rank_, tensor_type, value_type> &
    operator/=(const Real value);
    ///@}

    /**
     * @name Norm evaluation
     */
    ///@{
    /**
     * Return the Frobenius-norm of a tensor, i.e. the square root of the
     * sum of squares of all entries.
     */
    Real norm() const noexcept;

    /**
     * Return the square of the Frobenius-norm of a tensor, i.e.
     * the square root of the sum of squares of all entries.
     *
     * This function mainly exists because it makes computing the
     * norm simpler recursively, but may also be useful in other
     * contexts.
     */
    Real norm_square() const noexcept;
    ///@}

public:

    /** @name Dealing with the indices
     * @todo maybe these functions could be outside
     * TODO: maybe these functions could be outside
     */
    ///@{
    TensorIndex<rank_>
    flat_to_tensor_index(const int flat_index) const noexcept;

    /**
     * The last index moves faster.
     */
    int tensor_to_flat_index(const TensorIndex<rank_> &tensor_index) const noexcept;
    ///@}

private :
    SubTensor<self_t> tensor_[dim_== 0? 1: dim_];

};





/**
 * Specialization of the tensor class to handle Derivatives.
 * It fits the derivative of a function ... at a given point.
 * @todo complete the documentation.
 */
template<int dim, int range, int rank, int order>
using Derivatives =
    Conditional<
    order==0,
    Tensor<1,1,tensor::covariant,
    Tensor<range,rank,tensor::contravariant,Tdouble>>,
    Tensor<dim,order,tensor::covariant,
    Tensor<range,rank,tensor::contravariant,Tdouble>>
    >;


template<int dim, int range, int rank>
using Values = Tensor<range, rank, tensor::contravariant, Tdouble>;


/**
 * The <tt>Point</tt> class provides for a point or vector in a space with
 * arbitrary dimension <tt>dim</tt>.
 *
 * In igatools we use <tt>Point</tt> for:
 * - the type to evaluate a function
 * - the type for a vertex of the reference patch
 * - the return type for the values of mapping
 * @note technically  <tt>Point</tt> is not a class but
 * a template alias (or templated typedef).
 * @note In the dim==1 case the tensor is of rank 0, and in
 * all other cases of rank 1, this is so that we can use the
 * point as a Real in the dim ==1 case.
 *
 */
template< int dim >
using Point = Tensor<dim, 1, tensor::contravariant, Tdouble>;

/**
 *  Returns the tensor product of two rank one tensors.
 *  Mathematically denoted by \f$ value\_type = a \otimes b \f$,
 *  and defined by \f$ T(u) = a (b \cdot u) \f$.
 *
 * @relates Tensor
 */
template<class V1, class V2>
Tensor<V2::dim, 1, typename V2::tensor_t::co_type,Tensor<V1::dim, 1, typename V1::tensor_t,Tdouble> >
tensor_product(const V1 &a, const V2 &b)
{
    Tensor<V2::dim, 1, typename V2::tensor_t::co_type,Tensor<V1::dim, 1, typename V1::tensor_t,Tdouble>> R;

    for (int u = 0; u < V2::dim; ++u)
        R[u] = b[u] * a;

    return R;
}



/** @name tensor_arith_oper Tensor arithmetic operators */
/** @} */
/**
 * Adds two tensors of the same type.
 * @note When possible use *= instead as it
 * does not require the creation of new tensor.
 *
 * @relates Tensor
 */
template<class T>
EnableIf<T::is_tensor, T>
operator+(const T &A, const T &B) noexcept;

/**
 * Subtract two tensors.
 * @note When possible use *= instead as it
 * does not require the creation of new tensor.
 * @relates Tensor
 */
template<class T>
EnableIf<T::is_tensor,T>
operator-(const T &A, const T &B) noexcept;

/**
 * Multiply a tensor by a scalar.
 * @note When possible use *= instead as it
 * does not require the creation of new tensor.
 *
 * @relates Tensor
 */
template<class T>
EnableIf<T::is_tensor,T>
operator*(const T &A, const Real scalar) noexcept;

/**
 * Multiply a tensor by a scalar.
 * @note When possible use *= instead as it
 * does not require the creation of new tensor.
 *
 * @relates Tensor
 */
template<class T>
EnableIf<T::is_tensor,T>
operator*(const Real scalar, const T &A) noexcept;

/**
 * Divide a tensor by a scalar.
 * @note When possible use /= instead as it
 * does not require the creation of new tensor.
 *
 * @relates Tensor
 */
template<class T>
EnableIf<T::is_tensor,T>
operator/(const T &A, const Real scalar) noexcept;
/** @} */

template<class T>
EnableIf<T::is_tensor,T>
action(const T &A, const Tdouble &x)
{
    T B(A);
    B *= x[0];
    return B;
}




/**
 * Computes the action of tensor A on x,
 * returning the tensor A(x).
 */
template <class A_t, class V_t>
EnableIf<!std::is_same<V_t,Tdouble>::value,
         ActionTensor<A_t> >
         action(const A_t &A, const V_t &x)
{
    Assert(bool(std::is_same<typename A_t::tensor_t::co_type, typename V_t::tensor_t>::value),
           ExcMessage("Wrong tensor types in action"));
    ActionTensor<A_t> R;
    for (int i = 0; i < A_t::dim; ++i)
        R += action(A[i],x[i]);

    return R;
}


/**
 * Composition of two tensors for which is well defined.
 * R = S compose by T
 * this is like a contract 1 index */
template< class T1, class T2 >
Tensor<T2::dim, T2::rank, typename T2::tensor_t, SubTensor<T1> >
compose(const T1 &S, const T2 &T)
{
    Tensor<T2::dim, T2::rank, typename T2::tensor_t, SubTensor<T1> > R;
    for (int i = 0; i < T2::dim; ++i)
        R[i] = action(S,T[i]);

    return R;
}



/**
 * Computes the transpose of a tensor.
 * A   : V x W -> R, then
 * A^t : W x V -> R, and
 * A(v,w)=A^t(w,v).
 *
 */
template < class T >
EnableIf<T::is_tensor,Transpose<T>>
                                 transpose(const T &A)
{
    Transpose<T> B;

    const int size1 = T::size;
    const int size2 = T::value_t::size;

    for (int i=0; i<size1; ++i)
    {
        auto index_a = A.flat_to_tensor_index(i);
        for (int j=0; j<size2; ++j)
        {
            auto index_b = A(i).flat_to_tensor_index(j);
            B(index_b)(index_a) = A(index_a)(index_b);
        }
    }
    return B;
}

/**
 * Returns the symmetric tensor of A, i.e.
 * 0.5*(A+A^t)
 *
 * @relates Tensor
 */
template < class T >
EnableIf<T::is_tensor,T>
symmetric_tensor(const T &A)
{
    Assert(T::dim == T::value_t::dim, ExcMessage("Only for square tensors."));

    T B(A);
    for (int i=0; i<T::size; ++i)
        for (int j=0; j<T::size; ++j)
            B[i][j] += A[j][i];

    B *= 0.5;

    return B;
}

/**
* Returns the co-type tensor
*
* @relates Tensor
*/
template < class T >
EnableIf<T::is_tensor,CoTensor<T> >
co_tensor(const T &A)
{
    // we copy the memory of A in coA in order to avoid the
    // aliasing due to the use of reinterpret_cast
    CoTensor<T> coA;
    memcpy(&coA, &A, sizeof(A));
    return coA;
}


/**
 * Scalar product of two number.
 * It is the last step in the recursively defined scalar product for tensors.
 *
 * @relates Tensor
 */
inline
Real
scalar_product(const Tdouble a, const Tdouble b) noexcept
{
    return a[0]*b[0];
}

/**
 * Scalar product of two tensors of the same type.
 *
 * @relates Tensor
 */
template< class T >
inline
EnableIf<T::is_tensor,Real>
scalar_product(const T &t1, const T &t2) noexcept
{
    Real result = 0.;
    for (int i = 0; i < T::dim; ++i)
        result += scalar_product(t1[i],t2[i]);

    Assert(!std::isnan(result),ExcNotANumber());
    Assert(!std::isinf(result),ExcNumberNotFinite());

    return result;
}

/**
 * Contracts a tensor of rank k+1 with one of rank k
 * returning a tensor of rank 1.
 *
 * @relates Tensor
 */
template< class T1, class T2 >
inline
Tensor<T1::dim,1,typename T1::tensor_t,Tdouble>
contract_1(const T1 &t1, const T2 &t2) noexcept
{
    Tensor<T1::dim,1,typename T1::tensor_t,Tdouble> result;
    for (int i = 0; i < T1::dim; ++i)
        result[i] = scalar_product(t1[i],t2);

    return result;
}

template< class tensor_type1, class tensor_type2>
inline
Real
inverse_(const Tensor<0,1,tensor_type1,Tensor<0,1,tensor_type2,Tdouble>> &t,
         Tensor<0,1,typename tensor_type2::co_type,Tensor<0,1,typename tensor_type1::co_type,Tdouble>> &t_inv)
{
    return 1.;
}

template< class tensor_type1, class tensor_type2>
inline
Real
inverse_(const Tensor<1,1,tensor_type1,Tensor<1,1,tensor_type2,Tdouble>> &t,
         Tensor<1,1,typename tensor_type2::co_type,Tensor<1,1,typename tensor_type1::co_type,Tdouble>> &t_inv)
{
    const Real det = t[0][0];
    Assert(det != Real(0.0), ExcDivideByZero());

    const Real InvDet = 1.0 / det;

    t_inv[0][0] =  InvDet;

    return det;
}


template< class tensor_type1, class tensor_type2>
inline
Real
inverse_(const Tensor<2,1,tensor_type1,Tensor<2,1,tensor_type2,Tdouble>> &t,
         Tensor<2,1,typename tensor_type2::co_type,Tensor<2,1,typename tensor_type1::co_type,Tdouble>> &t_inv)
{
    const Real det = t[0][0] * t[1][1] - t[0][1] * t[1][0];
    Assert(det != Real(0.0), ExcDivideByZero());

    const Real InvDet = 1.0 / det;


    t_inv[0][0] = t[1][1] * InvDet;
    t_inv[0][1] = t[0][1] * (-InvDet);
    t_inv[1][0] = t[1][0] * (-InvDet);
    t_inv[1][1] = t[0][0] * InvDet;

    return det;
}




template< class type1, class type2>
inline
Real
inverse_(const Tensor<3,1,type1,Tensor<3,1,type2,Tdouble>> &t,
         Tensor<3,1,typename type2::co_type,Tensor<3,1,typename type1::co_type,Tdouble>> &t_inv)
{
    const Real t4 = t[0][0]*t[1][1];
    const Real t6 = t[0][0]*t[1][2];
    const Real t8 = t[0][1]*t[1][0];
    const Real t00 = t[0][2]*t[1][0];
    const Real t01 = t[0][1]*t[2][0];
    const Real t04 = t[0][2]*t[2][0];
    const Real det = (t4*t[2][2]-t6*t[2][1]-t8*t[2][2]+
                      t00*t[2][1]+t01*t[1][2]-t04*t[1][1]);
    Assert(det != Real(0.0), ExcDivideByZero());

    const Real t07 = 1.0/det;
    t_inv[0][0] = (t[1][1]*t[2][2]-t[1][2]*t[2][1])*t07;
    t_inv[0][1] = (t[0][2]*t[2][1]-t[0][1]*t[2][2])*t07;
    t_inv[0][2] = (t[0][1]*t[1][2]-t[0][2]*t[1][1])*t07;
    t_inv[1][0] = (t[1][2]*t[2][0]-t[1][0]*t[2][2])*t07;
    t_inv[1][1] = (t[0][0]*t[2][2]-t04)*t07;
    t_inv[1][2] = (t00-t6)*t07;
    t_inv[2][0] = (t[1][0]*t[2][1]-t[1][1]*t[2][0])*t07;
    t_inv[2][1] = (t01-t[0][0]*t[2][1])*t07;
    t_inv[2][2] = (t4-t8)*t07;


    return det;
}



template< int dim, int range >
inline
EnableIf<(dim==range),Real>
inverse(const Derivatives<dim, range, 1, 1> &DF,
        Derivatives<range, dim, 1, 1> &DF_inv)
{
    return inverse_(DF, DF_inv);
}


template< int dim, int range >
inline
EnableIf<(dim<range),Real>
inverse(const Derivatives<dim,range,1,1> &DF,
        Derivatives<range,dim,1,1> &DF_inv)
{
    const auto DF_t = co_tensor(transpose(DF));
    const auto G = compose(DF_t, DF);

    Derivatives<dim,dim,1,1> G_inv;
    const Real det = inverse_(G, G_inv);

    DF_inv = compose(G_inv, DF_t);

    return sqrt(det);
}


template< class tensor_type1, class tensor_type2>
inline
Real
det_(const Tensor<0,1,tensor_type1,Tensor<0,1,tensor_type2,Tdouble>> &t)
{
    return 1.;
}

template< class tensor_type1, class tensor_type2>
inline
Real
det_(const Tensor< 1, 1, tensor_type1, Tensor< 1, 1, tensor_type2, Tdouble > >  &t)
{
    return t[0][0];
}


template< class tensor_type1, class tensor_type2>
inline
Real
det_(const Tensor< 2, 1, tensor_type1, Tensor< 2, 1, tensor_type2, Tdouble > > &t)
{
    const Real det = t[0][0] * t[1][1] - t[0][1] * t[1][0];
    return det;
}




template< class type1, class type2>
inline
Real
det_(const Tensor< 3, 1, type1, Tensor< 3, 1, type2, Tdouble > > &t)
{
    const Real t4 = t[0][0]*t[1][1];
    const Real t6 = t[0][0]*t[1][2];
    const Real t8 = t[0][1]*t[1][0];
    const Real t00 = t[0][2]*t[1][0];
    const Real t01 = t[0][1]*t[2][0];
    const Real t04 = t[0][2]*t[2][0];
    const Real det = (t4*t[2][2]-t6*t[2][1]-t8*t[2][2]+
                      t00*t[2][1]+t01*t[1][2]-t04*t[1][1]);


    return det;
}



/**
 * Determinant of a square tensor.
 *
 * Specialization for dim_==1.
 *
 * @relates Tensor
 */
template< int dim, int range >
inline
EnableIf<(dim==range),Real>
determinant(const Derivatives<dim,range,1,1> &DF)
{
    return det_(DF);
}

/**
 * Determinant of a rectangular rank 2 tensor.
 *
 * Specialization for dim_==1.
 *
 * @relates Tensor
 */
template< int dim, int range >
inline
EnableIf<(dim<range),Real>
determinant(const Derivatives<dim, range,1,1> &DF)
{
    // G = DF_t o DF
    const auto DF_t = co_tensor(transpose(DF));
    const auto G = compose(DF_t, DF);
    const Real det = det_(G);


//    const Real det = det_(compose(co_tensor(transpose(DF)), DF));

    return sqrt(det);
}

/**
 * Trace of a rank_ 2 square tensor.
 *
 * @relates Tensor
 */
template< class T >
Real
trace(const T &A)
{
    Assert((T::dim == T::value_t::dim) || (T::rank==2),
           ExcMessage("Only for square tensors."));
    Real res = 0.;
    for (int i = 0; i < T::dim; ++i)
        res += A[i][i];
    return res;
}


/**
 * Output operator for tensor.
 *
 * @relates Tensor
*/
template< int dim_, int rank_, class tensor_type, class T >
LogStream &
operator<<(LogStream &out, const Tensor<dim_,rank_,tensor_type,T > &tensor)
{
    out << "[ ";
    for (int i=0; i<dim_; i++)
    {
        out << tensor[i] << " ";
    }
    out << "] ";

    return (out);
}

/**
 * Output operator for Tdouble.
 *
 * @relates Tensor
*/
inline
LogStream &
operator<<(LogStream &out, const Tdouble &tensor)
{
    out << tensor[0];
    return (out);
}


/**
 * @brief The <tt>TMatrix</tt> class provides for a statically allocated \p m x \p n matrix.
 *
 * @note technically  <tt>TMatrix</tt> is not a class but
 * a template alias (or templated typedef).
 */
template< int m, int n = m >
using TMatrix = Tensor< m, 1, tensor::raw, Tensor< n, 1, tensor::raw, Tdouble > >;

IGA_NAMESPACE_CLOSE

#include <igatools/base/tensor-inline.h>


#endif
