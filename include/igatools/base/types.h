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

#ifndef __IGA_TYPES_H_
#define __IGA_TYPES_H_

#include <igatools/base/config.h>
#include <igatools/base/exceptions.h>
#include <array>
#include <cmath>
#include <cstdint>
#include <iostream>


#include <boost/tti/has_member_function.hpp>
#include <boost/mpl/int.hpp>
#include <boost/mpl/vector.hpp>
#include <boost/variant.hpp>

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

    /**@name Applicable to elements of all grid-like containers */
    ///@{

    /** Quadrature points on the element */
    point          =    1L << 1,

    /** Differential of measure (length, area, volume) at the
     * quadrature points */
    measure        =    1L << 2,

    /** Differential of measure (length, area, volume) at the
     * quadrature points times the unit quadrature weights */
    w_measure      =    1L << 3,

    /** normal space */
    boundary_normal  =    1L << 4,

    /** normal space */
    outer_normal  =    1L << 5,

    /** curvatures */
    curvature     =    1L << 6,
    ///@}


    /** compute the values of the map */
    value              =  1L << 9,

    /** compute the gradients of the map */
    gradient           =    1L << 14,

    /** compute the hessians of the map */
    hessian            =    1L << 15,

    /** compute the (pseudo) inverse gradients of the map */
    inv_gradient       =    1L << 16,

    /** compute the (pseudo) inverse hessians of the map */
    inv_hessian        =    1L << 17,



#if 0
    ///@name Transformation (pushforward) related
    ///@{
    /** transform the values of basis functions */
    tran_value    =    1L << 22,

    /** transform the gradients of basis functions */
    tran_gradient =    1L << 23,

    /** transform the second derivatives of basis functions */
    tran_hessian  =    1L << 24,
    ///@}
#endif

    /** compute the divergences of basis functions */
    divergence    =    1L << 28
};


const std::array<ValueFlags, max_der> DerivativeFlags =
{ValueFlags::value, ValueFlags::gradient,  ValueFlags::hessian};

enum class TransformationFlags : std::int64_t
{
    /** transform nothing */
    tran_none    =    1L << 0,

    /** transform the values of basis functions */
    tran_value    =    1L << 1,

    /** transform the gradients of basis functions */
    tran_gradient =    1L << 2,

    /** transform the second derivatives of basis functions */
    tran_hessian  =    1L << 3
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
    h_curl = 2,

    /** \f[\vec v=\frac{D\vec{F}}{\det D\vec{F}}(\hat{\vec{v}}\circ\vec{F}^{-1})\f] */
    h_div = 3 ,

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
enum class LAPack : int
{
    /** Use the internal (igatools) linear algebra implementation.*/
    internal = 0,

    /** Use the Trilinos Tpetra linear algebra implementation.*/
    trilinos_tpetra = 1,

    /** Use the Trilinos Epetra linear algebra implementation.*/
    trilinos_epetra = 2,

    /** Use the PETSc linear algebra implementation.*/
    petsc = 3
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


/** Bitwise OR operator to use with the Flags of igatools */
template<class Flag>
constexpr
inline Flag
operator|(const Flag a, const Flag b)
{
    return (static_cast< Flag >(static_cast< int >(a) | static_cast< int >(b))) ;
}

/** Bitwise AND operator to use with the Flags of igatools */
template<class Flag>
inline Flag
operator&(const Flag a, const Flag b)
{
    return (static_cast< Flag >(static_cast< int >(a) & static_cast< int >(b))) ;
}

/** Bitwise XOR operator to use with the Flags of igatools */
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
    using underlying_type = typename std::underlying_type<Flag>::type;

    //Loop the value while there are still bits
    //Remove the end bit
    int count = 0;
    while (static_cast<underlying_type>(a) != 0)
    {
        a = a & static_cast<Flag>((static_cast<underlying_type>(a) - 1));
        count++;
    }

    return count;
}


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

/**
 * Macro used to generate the metafunction
 * <tt>has_member_function_print_info</tt>
 */
BOOST_TTI_HAS_MEMBER_FUNCTION(print_info)



/**
 * Type for specifying the kind of boundary condition on a face of a space-component.
 */
enum class BoundaryConditionType : int
{
    /** Dirichlet boundary condition. */
    Dirichlet = 0,

    /** Dirichlet homogeneous boundary condition. */
    DirichletHomogeneous = 1,

    /** Neumann boundary condition. */
    Neumann = 2,

    /** Neumann homogeneous boundary condition. */
    NeumannHomogeneous = 3,

    /** Robin boundary condition. */
    Robin = 4,

    /**
     * The boundary part of an interface, i.e. has a relation with another space-component face.
     * @sa InterfaceType
     * @sa MultiPatchSpace
     */
    Interface = 5
};

/**
 * Type for specifying the kind of interface between two faces of a space-component.
 * @sa MultiPatchSpace
 */
enum class InterfaceType : int
{
    /** No conditions on the interface.*/
    none = 0,

    /**
     * The interface is defined as strong C0 gluing between one side of each patch.
     * The dofs of the patches <em>will be not renumbered</em> in order to ensure this interface condition.
     */
    C0_strong = 1,

    /**
     * The interface is defined as strong C0 gluing between one side of each patch.
     * The dofs of the patches <em>will be renumbered</em> in order to ensure this interface condition.
     */
    C0_strong_renumbering = 2,

    /**
     * This interface is defined using Mortar gluing.
     *
     * @todo Complete the documentation.
     */
    Mortar = 3
};


/**
 * Bit field flags specifying the type of a linear constraint.
 */
enum class LinearConstraintType : int
{
    /** Lagrange multiplier. */
    lagrange = 1 << 0,

    /** Penalty. */
    penalty = 1 << 1,

    /** Augmented Lagrange multiplier. */
    augmented_lagrange = 1 << 2,

    /** Any of the above.*/
    any = lagrange | penalty | augmented_lagrange
};


enum class CopyPolicy : int
{
    /** Use the shallow copy policy: the pointer are copied. */
    shallow = 1,

    /** Use the deep copy policy: the pointer are allocated using the copy constructor. */
    deep = 2
};


/**
 * Enumerator for different kind of element properties.
 */
enum class ElementProperty : int
{
    /**
     * No-property indicator (useful for example if you want to iterate over all the elements of
     * the CartesianGrid, without taking into account the properties of the elements.).
     */
    none = 0,   //!< active

    /**
     * Active elements indicator (used for example in hierarchical spaces).
     */
    active = 1,   //!< active

    /**
    * Marked elements indicators.
    */
    marked = 2,   //!< marked

    /**
    * Influence elements indicators  (used for example in hierarchical spaces).
    */
    influence = 3,//!< influence

    /**
     * Number of different element properties allowed.
     */
    ENUM_SIZE = 3,//!< ENUM_SIZE
};



// TODO (pauletti, Nov 14, 2014): delete after gcc implements correct std::max
constexpr int max(int a, int b)
{
    return a>b ? a : b;
}


//---------------------------------------------------------------------------------------
template<template<int> class Q, int start, std::size_t N>
struct seq;

template<template<int> class Q, int start>
struct seq<Q, start, start>
{
    using type = boost::mpl::vector<Q<start>>;
};

template<template<int> class Q, int start, std::size_t N>
struct seq
{
    using v1 = typename seq<Q, start, N-1>::type;
    using type = typename boost::mpl::push_back<v1, Q<N>>::type;
};


template <template<int> class T,int dim>
using SubElemVariants = typename boost::make_variant_over<typename seq<T, iga::max(0, dim-num_sub_elem), dim>::type>::type;


template <int dim>
using Topology = boost::mpl::int_<dim>;


template <int dim>
using TopologyVariants = SubElemVariants<Topology,dim>;
//---------------------------------------------------------------------------------------




IGA_NAMESPACE_CLOSE

#endif /* __IGA_TYPES_H_ */
