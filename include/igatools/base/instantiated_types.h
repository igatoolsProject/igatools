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

#ifndef __INSTANTIATED_TYPES_H_
#define __INSTANTIATED_TYPES_H_

#include <igatools/base/config.h>

#include <boost/fusion/container/map.hpp>
#include <boost/fusion/container/vector.hpp>

#include <igatools/utils/safe_stl_map.h>

IGA_NAMESPACE_OPEN

template <int dim> class Grid;
template <int dim, int range> class GridFunction;
template <int dim, int codim> class Domain;
template <int dim, int range, int rank> class SplineSpace;
template <int dim, int range, int rank> class BSpline;
template <int dim, int range, int rank> class NURBS;
template <int dim, int range, int rank> class ReferenceSpaceBasis;
template <int dim, int range, int rank, int codim> class PhysicalSpaceBasis;
template <int dim, int codim, int range, int rank> class Function;

/**
 * @brief This is a helper struct for defining containers for all
 * the igatools instantiated types.
 *
 * This is a helper struct for defining containers for all
 * the igatools instantiated types.
 *
 * The types are stored as shared pointers into <tt>boost::fusion</tt>
 * vectors.
 *
 * @author P. Antolin
 * @date 2015
 */
struct InstantiatedTypes
{
private:

  /** Alias type for the shared pointer of a grid. */
  template <int dim>
  using GridPtr = std::shared_ptr<Grid<dim>>;

  /** Alias type for the map of indices and grid shared pointers. */
  template <int dim>
  using MapGridPtr = SafeSTLMap<Index, GridPtr<dim>>;

  /** Alias type for the shared pointer of a spline space. */
  template <int dim, int range, int rank>
  using SplineSpacePtr = std::shared_ptr<SplineSpace<dim, range, rank>>;

  /** Alias type for the map of indices and spline space shared pointers. */
  template <int dim, int range, int rank>
  using MapSplineSpacePtr = SafeSTLMap<Index, SplineSpacePtr<dim, range, rank>>;

  /** Alias for the reference space basis. */
  template <int dim, int range, int rank>
  using RefSpace = ReferenceSpaceBasis<dim, range, rank>;

  /** Alias type for the shared pointer of a reference space basis. */
  template <int dim, int range, int rank>
  using RefSpacePtr = std::shared_ptr<RefSpace<dim, range, rank>>;

  /** Alias type for the map of indices and reference space basis shared pointers. */
  template <int dim, int range, int rank>
  using MapRefSpacePtr = SafeSTLMap<Index, RefSpacePtr<dim, range, rank>>;

  /** Alias type for the shared pointer of a grid function. */
  template <int dim, int range>
  using GridFuncPtr = std::shared_ptr<GridFunction<dim, range>>;

  /** Alias type for the map of indices and grid function shared pointers. */
  template <int dim, int range>
  using MapGridFuncPtr = SafeSTLMap<Index, GridFuncPtr<dim, range>>;

  /** Alias type for the shared pointer of a domain. */
  template <int dim, int codim>
  using DomainPtr = std::shared_ptr<Domain<dim, codim>>;

  /** Alias type for the map of indices and domain shared pointers. */
  template <int dim, int codim>
  using MapDomainPtr = SafeSTLMap<Index, DomainPtr<dim, codim>>;

  /** Alias for physical space basis. */
  template <int dim, int range, int rank, int codim>
  using PhysSpace = PhysicalSpaceBasis<dim, range, rank, codim>;

  /** Alias type for the shared pointer of a physical space basis. */
  template <int dim, int range, int rank, int codim>
  using PhysSpacePtr = std::shared_ptr<PhysSpace<dim, range, rank, codim>>;

  /** Alias type for the map of indices and physical space basis shared pointers. */
  template <int dim, int range, int rank, int codim>
  using MapPhysSpacePtr = SafeSTLMap<Index, PhysSpacePtr<dim, range, rank, codim>>;

  /** Alias type for the shared pointer of a function. */
  template <int dim, int codim, int range, int rank>
  using FunctionPtr = std::shared_ptr<Function<dim, codim, range, rank>>;

  /** Alias type for the map of indices and physical space basis shared pointers. */
  template <int dim, int codim, int range, int rank>
  using MapFunctionPtr = SafeSTLMap<Index, FunctionPtr<dim, codim, range, rank>>;
#include <map>

public:
  /** All grid instantiations. */
  typedef boost::fusion::vector<std::shared_ptr<Grid<0>>,
                                std::shared_ptr<Grid<1>>,
                                std::shared_ptr<Grid<2>>,
                                std::shared_ptr<Grid<3>>> GridPtrs;

  /** All spline space instantiations. */
  typedef boost::fusion::vector<
          std::shared_ptr<SplineSpace<0,0,1>>,
          std::shared_ptr<SplineSpace<0,1,1>>,
          std::shared_ptr<SplineSpace<0,2,1>>,
          std::shared_ptr<SplineSpace<0,3,1>>,
          std::shared_ptr<SplineSpace<1,1,1>>,
          std::shared_ptr<SplineSpace<1,2,1>>,
          std::shared_ptr<SplineSpace<1,3,1>>,
          std::shared_ptr<SplineSpace<2,1,1>>,
          std::shared_ptr<SplineSpace<2,2,1>>,
          std::shared_ptr<SplineSpace<2,3,1>>,
          std::shared_ptr<SplineSpace<3,1,1>>,
          std::shared_ptr<SplineSpace<3,3,1>>> SplineSpacePtrs;

  /** All grid function instantiations. */
  typedef boost::fusion::vector<
          std::shared_ptr<GridFunction<0,0>>,
          std::shared_ptr<GridFunction<0,1>>,
          std::shared_ptr<GridFunction<0,2>>,
          std::shared_ptr<GridFunction<0,3>>,
          std::shared_ptr<GridFunction<1,1>>,
          std::shared_ptr<GridFunction<1,2>>,
          std::shared_ptr<GridFunction<1,3>>,
          std::shared_ptr<GridFunction<2,1>>,
          std::shared_ptr<GridFunction<2,2>>,
          std::shared_ptr<GridFunction<2,3>>,
          std::shared_ptr<GridFunction<3,1>>,
          std::shared_ptr<GridFunction<3,3>>> GridFunctionPtrs;

  /** All domain instantiations. */
  typedef boost::fusion::vector<
          std::shared_ptr<Domain<0,0>>,
          std::shared_ptr<Domain<0,1>>,
          std::shared_ptr<Domain<0,2>>,
          std::shared_ptr<Domain<0,3>>,
          std::shared_ptr<Domain<1,0>>,
          std::shared_ptr<Domain<1,1>>,
          std::shared_ptr<Domain<1,2>>,
          std::shared_ptr<Domain<2,0>>,
          std::shared_ptr<Domain<2,1>>,
          std::shared_ptr<Domain<3,0>>> DomainPtrs;

  /** All physical space basis instantiations. */
  typedef boost::fusion::vector<
          std::shared_ptr<PhysSpace<0,0,1,0>>,
          std::shared_ptr<PhysSpace<0,1,1,2>>,
          std::shared_ptr<PhysSpace<0,1,1,1>>,
          std::shared_ptr<PhysSpace<0,1,1,3>>,
          std::shared_ptr<PhysSpace<0,3,1,2>>,
          std::shared_ptr<PhysSpace<0,3,1,1>>,
          std::shared_ptr<PhysSpace<0,2,1,2>>,
          std::shared_ptr<PhysSpace<0,2,1,1>>,
          std::shared_ptr<PhysSpace<0,3,1,3>>,
          std::shared_ptr<PhysSpace<1,1,1,0>>,
          std::shared_ptr<PhysSpace<1,1,1,1>>,
          std::shared_ptr<PhysSpace<1,2,1,1>>,
          std::shared_ptr<PhysSpace<1,3,1,1>>,
          std::shared_ptr<PhysSpace<1,1,1,2>>,
          std::shared_ptr<PhysSpace<1,3,1,2>>,
          std::shared_ptr<PhysSpace<2,1,1,0>>,
          std::shared_ptr<PhysSpace<2,2,1,0>>,
          std::shared_ptr<PhysSpace<2,1,1,1>>,
          std::shared_ptr<PhysSpace<2,3,1,1>>,
          std::shared_ptr<PhysSpace<3,1,1,0>>,
          std::shared_ptr<PhysSpace<3,3,1,0>>,
          std::shared_ptr<PhysSpace<1,2,1,0>>,
          std::shared_ptr<PhysSpace<1,3,1,0>>,
          std::shared_ptr<PhysSpace<2,3,1,0>>,
          std::shared_ptr<PhysSpace<0,1,1,0>>,
          std::shared_ptr<PhysSpace<0,2,1,0>>,
          std::shared_ptr<PhysSpace<0,3,1,0>>> PhysSpacePtrs;

  /** All function instantiations. */
  typedef boost::fusion::vector<
          std::shared_ptr<Function<0,0,0,1>>,
          std::shared_ptr<Function<0,2,1,1>>,
          std::shared_ptr<Function<0,1,1,1>>,
          std::shared_ptr<Function<0,3,1,1>>,
          std::shared_ptr<Function<0,2,3,1>>,
          std::shared_ptr<Function<0,1,3,1>>,
          std::shared_ptr<Function<0,2,2,1>>,
          std::shared_ptr<Function<0,1,2,1>>,
          std::shared_ptr<Function<0,3,3,1>>,
          std::shared_ptr<Function<1,0,1,1>>,
          std::shared_ptr<Function<1,1,1,1>>,
          std::shared_ptr<Function<1,1,2,1>>,
          std::shared_ptr<Function<1,1,3,1>>,
          std::shared_ptr<Function<1,2,1,1>>,
          std::shared_ptr<Function<1,2,3,1>>,
          std::shared_ptr<Function<2,0,1,1>>,
          std::shared_ptr<Function<2,0,2,1>>,
          std::shared_ptr<Function<2,1,1,1>>,
          std::shared_ptr<Function<2,1,3,1>>,
          std::shared_ptr<Function<3,0,1,1>>,
          std::shared_ptr<Function<3,0,3,1>>,
          std::shared_ptr<Function<1,0,2,1>>,
          std::shared_ptr<Function<1,0,3,1>>,
          std::shared_ptr<Function<2,0,3,1>>,
          std::shared_ptr<Function<0,0,1,1>>,
          std::shared_ptr<Function<0,0,2,1>>,
          std::shared_ptr<Function<0,0,3,1>>> FunctionPtrs;

  /** Fusion map relating types and map of indices and shared pointers of the types. */
  typedef boost::fusion::map<
  boost::fusion::pair<Grid<0>, MapGridPtr<0>>,
  boost::fusion::pair<Grid<1>, MapGridPtr<1>>,
  boost::fusion::pair<Grid<2>, MapGridPtr<2>>,
  boost::fusion::pair<Grid<3>, MapGridPtr<3>>,
  boost::fusion::pair<SplineSpace<0,0,1>, MapSplineSpacePtr<0,0,1>>,
  boost::fusion::pair<SplineSpace<0,1,1>, MapSplineSpacePtr<0,1,1>>,
  boost::fusion::pair<SplineSpace<0,2,1>, MapSplineSpacePtr<0,2,1>>,
  boost::fusion::pair<SplineSpace<0,3,1>, MapSplineSpacePtr<0,3,1>>,
  boost::fusion::pair<SplineSpace<1,1,1>, MapSplineSpacePtr<1,1,1>>,
  boost::fusion::pair<SplineSpace<1,2,1>, MapSplineSpacePtr<1,2,1>>,
  boost::fusion::pair<SplineSpace<1,3,1>, MapSplineSpacePtr<1,3,1>>,
  boost::fusion::pair<SplineSpace<2,1,1>, MapSplineSpacePtr<2,1,1>>,
  boost::fusion::pair<SplineSpace<2,2,1>, MapSplineSpacePtr<2,2,1>>,
  boost::fusion::pair<SplineSpace<2,3,1>, MapSplineSpacePtr<2,3,1>>,
  boost::fusion::pair<SplineSpace<3,1,1>, MapSplineSpacePtr<3,1,1>>,
  boost::fusion::pair<SplineSpace<3,3,1>, MapSplineSpacePtr<3,3,1>>,
  boost::fusion::pair<RefSpace<0,0,1>, MapRefSpacePtr<0,0,1>>,
  boost::fusion::pair<RefSpace<0,1,1>, MapRefSpacePtr<0,1,1>>,
  boost::fusion::pair<RefSpace<0,2,1>, MapRefSpacePtr<0,2,1>>,
  boost::fusion::pair<RefSpace<0,3,1>, MapRefSpacePtr<0,3,1>>,
  boost::fusion::pair<RefSpace<1,1,1>, MapRefSpacePtr<1,1,1>>,
  boost::fusion::pair<RefSpace<1,2,1>, MapRefSpacePtr<1,2,1>>,
  boost::fusion::pair<RefSpace<1,3,1>, MapRefSpacePtr<1,3,1>>,
  boost::fusion::pair<RefSpace<2,1,1>, MapRefSpacePtr<2,1,1>>,
  boost::fusion::pair<RefSpace<2,2,1>, MapRefSpacePtr<2,2,1>>,
  boost::fusion::pair<RefSpace<2,3,1>, MapRefSpacePtr<2,3,1>>,
  boost::fusion::pair<RefSpace<3,1,1>, MapRefSpacePtr<3,1,1>>,
  boost::fusion::pair<RefSpace<3,3,1>, MapRefSpacePtr<3,3,1>>,
  boost::fusion::pair<GridFunction<0,0>, MapGridFuncPtr<0,0>>,
  boost::fusion::pair<GridFunction<0,1>, MapGridFuncPtr<0,1>>,
  boost::fusion::pair<GridFunction<0,2>, MapGridFuncPtr<0,2>>,
  boost::fusion::pair<GridFunction<0,3>, MapGridFuncPtr<0,3>>,
  boost::fusion::pair<GridFunction<1,1>, MapGridFuncPtr<1,1>>,
  boost::fusion::pair<GridFunction<1,2>, MapGridFuncPtr<1,2>>,
  boost::fusion::pair<GridFunction<1,3>, MapGridFuncPtr<1,3>>,
  boost::fusion::pair<GridFunction<2,1>, MapGridFuncPtr<2,1>>,
  boost::fusion::pair<GridFunction<2,2>, MapGridFuncPtr<2,2>>,
  boost::fusion::pair<GridFunction<2,3>, MapGridFuncPtr<2,3>>,
  boost::fusion::pair<GridFunction<3,1>, MapGridFuncPtr<3,1>>,
  boost::fusion::pair<GridFunction<3,3>, MapGridFuncPtr<3,3>>,
  boost::fusion::pair<Domain<0,0>, MapDomainPtr<0,0>>,
  boost::fusion::pair<Domain<0,1>, MapDomainPtr<0,1>>,
  boost::fusion::pair<Domain<0,2>, MapDomainPtr<0,2>>,
  boost::fusion::pair<Domain<0,3>, MapDomainPtr<0,3>>,
  boost::fusion::pair<Domain<1,0>, MapDomainPtr<1,0>>,
  boost::fusion::pair<Domain<1,1>, MapDomainPtr<1,1>>,
  boost::fusion::pair<Domain<1,2>, MapDomainPtr<1,2>>,
  boost::fusion::pair<Domain<2,0>, MapDomainPtr<2,0>>,
  boost::fusion::pair<Domain<2,1>, MapDomainPtr<2,1>>,
  boost::fusion::pair<Domain<3,0>, MapDomainPtr<3,0>>,
  boost::fusion::pair<PhysSpace<0,0,1,0>, MapPhysSpacePtr<0,0,1,0>>,
  boost::fusion::pair<PhysSpace<0,1,1,2>, MapPhysSpacePtr<0,1,1,2>>,
  boost::fusion::pair<PhysSpace<0,1,1,1>, MapPhysSpacePtr<0,1,1,1>>,
  boost::fusion::pair<PhysSpace<0,1,1,3>, MapPhysSpacePtr<0,1,1,3>>,
  boost::fusion::pair<PhysSpace<0,3,1,2>, MapPhysSpacePtr<0,3,1,2>>,
  boost::fusion::pair<PhysSpace<0,3,1,1>, MapPhysSpacePtr<0,3,1,1>>,
  boost::fusion::pair<PhysSpace<0,2,1,2>, MapPhysSpacePtr<0,2,1,2>>,
  boost::fusion::pair<PhysSpace<0,2,1,1>, MapPhysSpacePtr<0,2,1,1>>,
  boost::fusion::pair<PhysSpace<0,3,1,3>, MapPhysSpacePtr<0,3,1,3>>,
  boost::fusion::pair<PhysSpace<1,1,1,0>, MapPhysSpacePtr<1,1,1,0>>,
  boost::fusion::pair<PhysSpace<1,1,1,1>, MapPhysSpacePtr<1,1,1,1>>,
  boost::fusion::pair<PhysSpace<1,2,1,1>, MapPhysSpacePtr<1,2,1,1>>,
  boost::fusion::pair<PhysSpace<1,3,1,1>, MapPhysSpacePtr<1,3,1,1>>,
  boost::fusion::pair<PhysSpace<1,1,1,2>, MapPhysSpacePtr<1,1,1,2>>,
  boost::fusion::pair<PhysSpace<1,3,1,2>, MapPhysSpacePtr<1,3,1,2>>,
  boost::fusion::pair<PhysSpace<2,1,1,0>, MapPhysSpacePtr<2,1,1,0>>,
  boost::fusion::pair<PhysSpace<2,2,1,0>, MapPhysSpacePtr<2,2,1,0>>,
  boost::fusion::pair<PhysSpace<2,1,1,1>, MapPhysSpacePtr<2,1,1,1>>,
  boost::fusion::pair<PhysSpace<2,3,1,1>, MapPhysSpacePtr<2,3,1,1>>,
  boost::fusion::pair<PhysSpace<3,1,1,0>, MapPhysSpacePtr<3,1,1,0>>,
  boost::fusion::pair<PhysSpace<3,3,1,0>, MapPhysSpacePtr<3,3,1,0>>,
  boost::fusion::pair<PhysSpace<1,2,1,0>, MapPhysSpacePtr<1,2,1,0>>,
  boost::fusion::pair<PhysSpace<1,3,1,0>, MapPhysSpacePtr<1,3,1,0>>,
  boost::fusion::pair<PhysSpace<2,3,1,0>, MapPhysSpacePtr<2,3,1,0>>,
  boost::fusion::pair<PhysSpace<0,1,1,0>, MapPhysSpacePtr<0,1,1,0>>,
  boost::fusion::pair<PhysSpace<0,2,1,0>, MapPhysSpacePtr<0,2,1,0>>,
  boost::fusion::pair<PhysSpace<0,3,1,0>, MapPhysSpacePtr<0,3,1,0>>,
  boost::fusion::pair<Function<0,0,0,1>, MapFunctionPtr<0,0,0,1>>,
  boost::fusion::pair<Function<0,2,1,1>, MapFunctionPtr<0,2,1,1>>,
  boost::fusion::pair<Function<0,1,1,1>, MapFunctionPtr<0,1,1,1>>,
  boost::fusion::pair<Function<0,3,1,1>, MapFunctionPtr<0,3,1,1>>,
  boost::fusion::pair<Function<0,2,3,1>, MapFunctionPtr<0,2,3,1>>,
  boost::fusion::pair<Function<0,1,3,1>, MapFunctionPtr<0,1,3,1>>,
  boost::fusion::pair<Function<0,2,2,1>, MapFunctionPtr<0,2,2,1>>,
  boost::fusion::pair<Function<0,1,2,1>, MapFunctionPtr<0,1,2,1>>,
  boost::fusion::pair<Function<0,3,3,1>, MapFunctionPtr<0,3,3,1>>,
  boost::fusion::pair<Function<1,0,1,1>, MapFunctionPtr<1,0,1,1>>,
  boost::fusion::pair<Function<1,1,1,1>, MapFunctionPtr<1,1,1,1>>,
  boost::fusion::pair<Function<1,1,2,1>, MapFunctionPtr<1,1,2,1>>,
  boost::fusion::pair<Function<1,1,3,1>, MapFunctionPtr<1,1,3,1>>,
  boost::fusion::pair<Function<1,2,1,1>, MapFunctionPtr<1,2,1,1>>,
  boost::fusion::pair<Function<1,2,3,1>, MapFunctionPtr<1,2,3,1>>,
  boost::fusion::pair<Function<2,0,1,1>, MapFunctionPtr<2,0,1,1>>,
  boost::fusion::pair<Function<2,0,2,1>, MapFunctionPtr<2,0,2,1>>,
  boost::fusion::pair<Function<2,1,1,1>, MapFunctionPtr<2,1,1,1>>,
  boost::fusion::pair<Function<2,1,3,1>, MapFunctionPtr<2,1,3,1>>,
  boost::fusion::pair<Function<3,0,1,1>, MapFunctionPtr<3,0,1,1>>,
  boost::fusion::pair<Function<3,0,3,1>, MapFunctionPtr<3,0,3,1>>,
  boost::fusion::pair<Function<1,0,2,1>, MapFunctionPtr<1,0,2,1>>,
  boost::fusion::pair<Function<1,0,3,1>, MapFunctionPtr<1,0,3,1>>,
  boost::fusion::pair<Function<2,0,3,1>, MapFunctionPtr<2,0,3,1>>,
  boost::fusion::pair<Function<0,0,1,1>, MapFunctionPtr<0,0,1,1>>,
  boost::fusion::pair<Function<0,0,2,1>, MapFunctionPtr<0,0,2,1>>,
  boost::fusion::pair<Function<0,0,3,1>, MapFunctionPtr<0,0,3,1>>> ObjectTypes;

};

IGA_NAMESPACE_CLOSE

#endif /*__INSTANTIATED_TYPES_H_ */
