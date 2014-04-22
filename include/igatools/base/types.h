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

#ifndef __IGA_TYPES_H_
#define __IGA_TYPES_H_

#include <igatools/base/config.h>
#include <igatools/base/exceptions.h>
#include <array>
#include <vector>
#include <cmath>
#include <cstdint>
#include <iostream>

IGA_NAMESPACE_OPEN

/**
 * This file contains the declaration of types used in the library
 * as well as some utility functions associated with these types.
 */

/**
 * Real type.
 */
#ifdef REAL_IS_LONG_DOUBLE
using Real = long double;
#else
using Real = double;
#endif

/**
 * Unsigned integer.
 */
typedef unsigned int uint;

/**
 * Long integer.
 */
typedef unsigned long long ulongint;

/**
 * This type is used in several places of the library for
 * indexing a product like structure, to simplify developers life
 * if it ever becomes a class on its own we use int_array
 */
template<int dim>
using int_array = std::array<int, dim >;

/**
 * This type is used for array of vectors such as:
 * - raw knot vectors
 * - raw multiplicity vectors
 */
template< class T, int dim >
using vector_array = std::array< std::vector<T>, dim >;

/**
 * This function is used to initialize  an array at
 * construction time to a uniform value.
 */
template <class T, int dim>
inline
std::array<T,dim> filled_array(const T v)
{
    std::array<T,dim> res;
    res.fill(v);
    return res;
}

/**
 * Bounding Box, a dim-dimensional rectangular
 * box described by the intervals.
 * eg. BBox<2> box {{0,0.5},{1,2}}.
 */
template<int dim>
using BBox = std::array<std::array<Real, 2>, dim>;

/**
 * Id used for the boundary indicators,
 * Negative ids are reserved for gluing.
 */
using boundary_id = int;

/**
 * Id use to identify patches, could be interpret as
 * a material id.
 */
using patch_id = int;


/**
 * Id to assign to elements inside a patch.
 */
using element_id = int;


/**
 * Type for the flat indices of igatools container.
 */
using Index = int;

/**
 * Type for the size of igatools containers
 */
using Size = int;

/**
 * Bit field flags for specifying which element values will be needed
 * through the iterator on the grid like containers.
 * These are used in the init_value of the element iterators
 * as a uniform way of handling the caches
 * for efficient computations.
 * @todo put link to iterators_accessors.dox
 */
enum class ValueFlags : std::int64_t
{
    /** Fill nothing */
    none           =    0,

    ///@name Applicable to elements of all grid-like containers
    ///@{
    /** Quadrature points on the element */
    point          =    1L << 1,

    /** Differential of measure (length, area, volume) at the
     * quadrature points */
    measure        =    1L << 2,

    /** Differential of measure (length, area, volume) at the
     * quadrature points times the unit quadrature weights */
    w_measure      =    1L << 3,

    /** Quadrature points on a face of the element */
    face_point     =    1L << 4,

    /** Differential of face measure (length, area, volume) at the
     * quadrature points */
    face_measure   =    1L << 5,

    /** Differential of face measure (length, area, volume) at the
     * quadrature points times the unit quadrature weights */
    face_w_measure =    1L << 6,

    /** Face outer (element boundary) normal at the quadrature points */
    face_normal    =    1L << 7,

    ///@}

    ///@name Reference element related values
    ///@{
    /** reference element measure (one single value) */
//    ref_elem_measure       =    1L << 8,

    /** reference element face measure (one single value) */
//    ref_elem_face_measure  =    1L << 9,

    /** reference element coordinate lengths */
//    ref_elem_coord_length  =    1L << 10,

    /** reference element face outer normal (one single value) */
//    ref_elem_face_normal   =    1L << 11,

    ///@}

    ///@name Applicable to mapping, push-forward and physical space
    ///@{
    /** compute the element normal spaces */
    elem_normal            =    1L << 12,

    /** compute the element curvatures */
    elem_curvature         =    1L << 13,
    ///@}

    ///@name Mapping related values
    ///@{
    /** compute the values of the map */
    map_value              =    point,

    /** compute the gradients of the map */
    map_gradient           =    1L << 14,

    /** compute the gradients of the map */
    map_hessian            =    1L << 15,

    /** compute the gradients of the map */
    map_inv_gradient       =    1L << 16,

    /** compute the gradients of the map */
    map_inv_hessian        =    1L << 17,

    /** compute the values of the map */
    map_face_value         =    face_point,

    /** compute the gradients of the map */
    map_face_gradient      =    1L << 18,

    /** compute the gradients of the map */
    map_face_hessian       =    1L << 19,

    /** compute the gradients of the map */
    map_face_inv_gradient  =    1L << 20,

    /** compute the gradients of the map */
    map_face_inv_hessian   =    1L << 21,

    ///@}


    ///@name Transformation (pushforward) related
    ///@{
    /** transform the values of basis functions */
    tran_value    =    1L << 22,

    /** transform the gradients of basis functions */
    tran_gradient =    1L << 23,

    /** transform the second derivatives of basis functions */
    tran_hessian  =    1L << 24,

    ///@}


    ///@name Basis function (space) related
    ///@{
    /** compute values of basis functions */
    value         =    1L << 25,

    /** compute the gradients of basis functions */
    gradient      =    1L << 26,

    /** compute the second derivatives of basis functions */
    hessian       =    1L << 27,

    /** compute the second derivatives of basis functions */
    divergence    =    1L << 28,

    /** compute values of basis functions */
    face_value    =    1L << 29,

    /** compute the gradients of basis functions */
    face_gradient =    1L << 30,

    /** compute the second derivatives of basis functions */
    face_hessian  =    1L << 31,

    /** compute the second derivatives of basis functions */
    face_divergence    =    1L << 32

                            ///@}
};


/**
 * Type of transformation for values and derivatives of basis functions.
 * If \f$\vec F\f$ is the mapping, the transform functions are defined
 * for each class as documented below.
 */
enum class Transformation : int
{
    /** \f[ \phi = \hat\phi \circ \vec{F}^{-1}\f] */
    h_grad = 1,

    /** \f[\vec u = D\vec{F}^{-T}(\hat{\vec{u}} \circ \vec{F}^{-1})\f] */
    h_div = 2,

    /** \f[\vec v=\frac{D\vec{F}}{\det D\vec{F}}(\hat{\vec{v}}\circ\vec{F}^{-1})\f] */
    h_curl = 3 ,

    /** \f[ \varphi = \frac{\hat\varphi \circ \vec{F}^{-1}}{\det D\vec{F}} \f] */
    l_2 = 4
};

/**
 *  Bit field flags for specifying which kind of error-norm must be evaluated.
 */
enum class Norm
{
    /** use of the L2 norm */
    L2 = 1 << 0,

    /** use of the H1 norm */
    H1 = 1 << 1,

    /** use of the H1 semi-norm */
    H1_semi = 1 << 2

};


/**
 * Type for specifying which linear algebra package to use.
 */
enum class LinearAlgebraPackage
{
    /** Use the internal (igatools) linear algebra implementation.*/
    internal = 0,

    /** Use the Trilinos linear algebra implementation.*/
    trilinos = 1,

    /** Use the PETSc linear algebra implementation.*/
    petsc = 2
};

/**
 * Bit field flags specifying the verbosity level of the output.
 */
enum class VerbosityLevel : int
{
    /** Normal level. */
    normal = 1 << 0,

    /** The output contains information about the memory address of the objects involved. */
    debug = 1 << 1
};

/**
 * Generic convert function that can be used to convert any enum class to
 * its underlying integral type.
 *
 * Example:
   @code
   auto value = to_integral(enum_class::enum_element);

   auto redValue = to_integral(Color::Red); //where Color is an enum class! *
   @endcode
 * @author M.Martinelli
 * @date 2013
 */
template<typename E>
constexpr auto
to_integral(E e) -> typename std::underlying_type<E>::type
{
    return static_cast<typename std::underlying_type<E>::type>(e);
}

/**
 * Return the number of elements in any enum class.
 * @author M.Martinelli
 * @date 2013
 */
template<typename E>
constexpr auto
get_enum_size() -> typename std::underlying_type<E>::type
{
    return to_integral(E::ENUM_SIZE) ;
}


/** Bitwise or operator to use with the Flags of igatools */
template<class Flag>
constexpr
inline Flag
operator|(const Flag a, const Flag b)
{
    return (static_cast< Flag >(static_cast< int >(a) | static_cast< int >(b))) ;
}

/** Bitwise and operator to use with the Flags of igatools */
template<class Flag>
inline Flag
operator&(const Flag a, const Flag b)
{
    return (static_cast< Flag >(static_cast< int >(a) & static_cast< int >(b))) ;
}

/** Bitwise xor operator to use with the Flags of igatools */
template<class Flag>
inline Flag
operator^(const Flag a, const Flag b)
{
    return (static_cast< Flag >(static_cast< int >(a) ^ static_cast< int >(b))) ;
}

/** Operator to add a flag of a given flag */
template<class Flag>
inline void
operator|=(Flag &a, const Flag b)
{
    a = a|b ;
}

/** Checks if a given flag a contains the given flag b */
template<class Flag>
inline bool contains(const Flag &a, const Flag b)
{
    return ((a & b) ==  b);
}

//TODO(pauletti, Feb 19, 2014): the item below should be documented
template<class Flag>
int bitcount(Flag a)
{
    //Loop the value while there are still bits
    //Remove the end bit
    int count = 0;
    while (a != 0)
    {
        a = a & static_cast< Flag >((static_cast< int >(a) - 1));
        count++;
    }

    return count;
}


//TODO(pauletti, Feb 19, 2014): the item below should be documented
inline
std::ostream &operator<< (std::ostream &stream, const ValueFlags &flag)
{
    return (stream << static_cast< int >(flag));
}



/**
 * Const expression power
 */
constexpr int constexpr_pow(int a, int b)
{
    return (b>=1) ? a * pow(a,b-1) : 1 ;
}



/**
 * Const expression factorial
 */
constexpr int constexpr_factorial(int a)
{
    return (a>1) ? a * constexpr_factorial(a-1) : 1 ;
}



/**
 * Const expression binomial coefficient
 */
constexpr Real constexpr_binomial_coefficient(int a, int b)
{
    return Real(constexpr_factorial(a)) /
           Real(constexpr_factorial(b) * constexpr_factorial(a-b));
}



/**
 * Constants to use for the state of the iterators
 */
namespace IteratorState
{
static const int pass_the_end = -1;
static const int invalid = -2;
}



/**
 * Return the type T if C is true, the type F if C is false
 */
template<bool C, typename T, typename F>
using Conditional = typename std::conditional<C,T,F>::type;


template<bool B, typename T = void>
using EnableIf = typename std::enable_if<B,T>::type;


IGA_NAMESPACE_CLOSE

#endif /* __IGA_TYPES_H_ */
