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
#include <igatools/utils/safe_stl_array.h>

#include <boost/mpl/map.hpp>
#include <boost/mpl/int.hpp>

#include<map>
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

enum class CacheFlags
{
  /** Fill nothing */
  none           =    0,

  //public element information c
  /** Quadrature points on the element */
  point          =    1L << 1,

  /** Quadrature weigths on the element */
  weight         =    1L << 2
};

struct activate
{
  using FlagsToCache = std::map<Flags, CacheFlags>;
  static FlagsToCache grid;
};

/**
 * Alias used to define the container for the points in the cache.
 */
class _Point
{
public:
  static const std::string name;
  static const auto flag = CacheFlags::point;
};

/**
 * Alias used to define the container for the quadrature weights in the cache.
 */
class _Weight
{
public:
  static const std::string name;
  static const auto flag = CacheFlags::weight;
};

} // end namespace grid_element
//---------------------------------------------------------------------------------------


//------------------------------------------------------------------------------
namespace grid_function_element
{
/** Quantities that can be requested from a grid_function element */
enum class Flags
{
  none           =    0,

  D0       =    1L << 1,

  D1       =    1L << 3,

  D2       =    1L << 4,

  D3       =    1L << 5
};


static const SafeSTLArray<Flags, 3> all_flags = {Flags::D0, Flags::D1, Flags::D2};

/** Auxiliary quantities stored in a local cache */
enum class CacheFlags
{
  none           =    0,

  D0       =    1L << 1,

  D1       =    1L << 3,

  D2       =    1L << 4,

  D3       =    1L << 5
};


struct activate
{
  using FlagsToCache = std::map<Flags, CacheFlags>;
  static FlagsToCache grid_function;

  using FlagsToGrid = std::map<Flags, grid_element::Flags>;
  static FlagsToGrid grid;
};

template <int order>
struct _D
{
  static const std::string name;
  static const CacheFlags flag;
};

}




namespace domain_element
{
/** Quantities that can be requested from a domain element */
enum class Flags
{
  /** Fill nothing */
  none           =    0,

  point          =    1L << 1,

  measure        =    1L << 2,

  w_measure      =    1L << 3,

  ext_normal     =    1L << 4,

  jacobian       =    1L << 5,

  inv_jacobian   =    1L << 6,

  hessian        =    1L << 7,

  inv_hessian    =    1L << 8,

  boundary_normal =   1L << 9
};

static const SafeSTLArray<Flags,9> all_flags =
{
  Flags::point,
  Flags::measure,
  Flags::w_measure,
  Flags::ext_normal,
  Flags::jacobian,
  Flags::inv_jacobian,
  Flags::hessian,
  Flags::inv_hessian,
  Flags::boundary_normal
};

/** Auxiliary quantities stored in a local cache */
enum class CacheFlags
{
  none           =    0,

  measure        =    1L << 1,

  inv_jacobian   =    1L << 2,

  inv_hessian    =    1L << 3,

  boundary_normal =   1L << 4,

  ext_normal     =    1L << 5
};


struct activate
{
  using FlagsToCache = std::map<Flags, CacheFlags>;
  static FlagsToCache domain;

  using FlagsToGridFunc = std::map<Flags, grid_function_element::Flags>;
  static FlagsToGridFunc grid_func;

  using FlagsToGrid = std::map<Flags, grid_element::Flags>;
  static FlagsToGrid grid;
};



struct _Measure
{
  static const std::string name;
  static const auto flag = CacheFlags::measure;
};

struct _InvJacobian
{
  static const std::string name;
  static const auto flag = CacheFlags::inv_jacobian;
};

struct _InvHessian
{
  static const std::string name;
  static const auto flag = CacheFlags::inv_hessian;
};

struct _BoundaryNormal
{
  static const std::string name;
  static const auto flag = CacheFlags::boundary_normal;
};


struct _ExtNormal
{
  static const std::string name;
  static const auto flag = CacheFlags::ext_normal;
};


}

//---------------------------------------------------------------------------------------
/**
 * @name Admissible flags for the classes that are derived from FunctionElement.
 */
namespace function_element
{
/** Quantities that can be requested from a function element */
enum class Flags
{
  /** Fill nothing */
  none           =    0,

  value          =    1L << 1,

  gradient       =    1L << 2,

  D2             =    1L << 3

};


/** Auxiliary quantities stored in a local cache */
enum class CacheFlags
{
  none           =    0,

  value          =    1L << 1,

  gradient       =    1L << 2,

  D2             =    1L << 3
};


struct activate
{
  using FlagsToCache = std::map<Flags, CacheFlags>;
  static FlagsToCache function;

  using FlagsToDomain = std::map<Flags, domain_element::Flags>;
  static FlagsToDomain domain;
};

struct _Value
{
  static const std::string name;
  static const auto flag = CacheFlags::value;
};

struct _Gradient
{
  static const std::string name;
  static const auto flag = CacheFlags::gradient;
};

struct _D2
{
  static const std::string name;
  static const auto flag = CacheFlags::D2;
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
  gradient       =    1L << 2, //!< gradient

  /** Basis functions hessian */
  hessian        =    1L << 3,  //!< hessian

  /** Basis functions divergence */
  divergence     =    1L << 4, //!< divergence

  point          = 1L << 5, //!< point in the domain

  w_measure      = 1L << 6 //!< quad. weight * domain measure element
};

#if 0
static const SafeSTLArray<Flags,6> all_flags =
{
  Flags::value,
  Flags::gradient,
  Flags::hessian,
  Flags::divergence,
  Flags::point,
  Flags::w_measure
};


enum class CacheFlags
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


struct activate
{
  using SpaceFlagsToSpaceCache = std::map<Flags, CacheFlags>;
  static SpaceFlagsToSpaceCache phys_space_flags_to_phys_space_cache;

  using SpaceFlagsToDomainFlags = std::map<Flags,domain_element::Flags>;
  static SpaceFlagsToDomainFlags phys_space_flags_to_domain_flags;

  using SpaceCacheToDomainFlags = std::map<CacheFlags,domain_element::Flags>;
  static SpaceCacheToDomainFlags phys_space_cache_to_domain_flags;
};
#endif

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
