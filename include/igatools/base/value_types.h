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

#ifndef VALUE_TYPES_H_
#define VALUE_TYPES_H_

#include <igatools/base/config.h>

#include <boost/mpl/map.hpp>
#include <boost/mpl/int.hpp>


#include <string>

IGA_NAMESPACE_OPEN

//---------------------------------------------------------------------------------------
namespace grid_element
{
enum class Flags
{
  /** Fill nothing */
  none           =    0,

  /** Quadrature points on the element */
  point          =    1L << 1,

  /** Quadrature weigths on the element */
  weight    =    1L << 2
};

/**
 * Alias used to define the container for the points in the cache.
 */
class _Point
{
public:
  static const std::string name;
  static const auto flag = Flags::point;
};

/**
 * Alias used to define the container for the quadrature weights in the cache.
 */
class _Weight
{
public:
  static const std::string name;
  static const auto flag = Flags::weight;
};

} // end namespace grid_element
//---------------------------------------------------------------------------------------



namespace domain_element
{
enum class Flags
{
  /** Fill nothing */
  none           =    0,

  //public element information c
  /** Quadrature points on the element */
  point          =    1L << 1,

  /** Quadrature weigths on the element */
  measure        =    1L << 2,

  /** Quadrature weigths on the element */
  w_measure      =    1L << 3,

  // internal cache flags
  gradient       =    1L << 4
};

/**
 * Alias used to define the container for the points in the cache.
 */
class _Point
{
public:
  static const std::string name;
  static const auto flag = Flags::point;
};

/**
 * Alias used to define the container for the quadrature weights in the cache.
 */
class _Measure
{
public:
  static const std::string name;
  static const auto flag = Flags::measure;
};

struct _Gradient
{
  static const std::string name;
  static const auto flag = Flags::gradient;
};


} // end namespace grid_element
//------------------------------------------------------------------------------


//---------------------------------------------------------------------------------------
/**
 * @name Admissible flags for the classes that are derived from FunctionElement.
 */
namespace function_element
{
enum class Flags
{
  /** Fill nothing */
  none           =    0,

  /** Function values */
  value          =    1L << 1,

  /** Function gradients */
  gradient       =    1L << 2
};

struct _Value
{
  static const std::string name;
  static const auto flag = Flags::value;
};

struct _Gradient
{
  static const std::string name;
  static const auto flag = Flags::gradient;
};

} // end namespace function_element
//---------------------------------------------------------------------------------------





//---------------------------------------------------------------------------------------
/**
 * @name Admissible flags for the classes that are derived from SpaceElement.
 */
namespace space_element
{
enum class Flags
{
  /** Fill nothing */
  none           =    0,          //!< none

  /** Basis functions value */
  value          =    1L << 1,    //!< value

  /** Basis functions gradient */
  gradient          =    1L << 2, //!< gradient

  /** Basis functions hessian */
  hessian          =    1L << 3,  //!< hessian

  /** Basis functions divergence */
  divergence          =    1L << 4//!< divergence
};

/**
 * Alias used to define/select the container for the basis function values in the cache.
 */
class _Value
{
public:
  static const std::string name;
  static const auto flag = Flags::value;
};

/**
 * Alias used to define/select the container for the basis function gradients in the cache.
 */
class _Gradient
{
public:
  static const std::string name;
  static const auto flag = Flags::gradient;
};

/**
 * Alias used to define/select the container for the basis function hessians in the cache.
 */
class _Hessian
{
public:
  static const std::string name;
  static const auto flag = Flags::hessian;
};

/**
 * Alias used to define the container for the basis function divergences in the cache.
 */
class _Divergence
{
public:
  static const std::string name;
  static const auto flag = Flags::divergence;
};

} // end namespace space_element
//---------------------------------------------------------------------------------------




IGA_NAMESPACE_CLOSE

#endif //#ifndef  VALUE_TYPES_H_
