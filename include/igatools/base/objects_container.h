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

#ifndef __OBJECTS_CONTAINER_H_
#define __OBJECTS_CONTAINER_H_

#include <igatools/base/config.h>

#include <igatools/utils/safe_stl_vector.h>

#include <boost/fusion/container/map.hpp>
#include <boost/fusion/container/vector.hpp>

IGA_NAMESPACE_OPEN

template <int dim> class Grid;
template <int dim, int range, int rank> class SplineSpace;
template <int dim, int range, int rank> class ReferenceSpaceBasis;
template <int dim, int range, int rank> class BSpline;
template <int dim, int range, int rank> class NURBS;
template <int dim, int range> class GridFunction;
template <int dim, int codim> class Domain;
template <int dim, int range, int rank, int codim> class PhysicalSpaceBasis;
template <int dim, int codim, int range, int rank> class Function;

/**
 * @brief To be documented.
 *
 * @author P. Antolin 2015
 *
 * @todo document more
 */
class ObjectsContainer
{
private:

  /** Type for current class. */
  using self_t = ObjectsContainer;

  /** Type for the shared pointer of a grid. */
  template <int dim>
  using GridPtr_t = std::shared_ptr<Grid<dim>>;

  /** Type for the map of indices and grid shared pointers. */
  template <int dim>
  using MapGridPtr_t = std::map<Index, GridPtr_t<dim>>;

  /** Type for the shared pointer of a spline space. */
  template <int dim, int range, int rank>
  using SplineSpacePtr_t = std::shared_ptr<SplineSpace<dim, range, rank>>;

  /** Type for the map of indices and spline space shared pointers. */
  template <int dim, int range, int rank>
  using MapSplineSpacePtr_t = std::map<Index, SplineSpacePtr_t<dim, range, rank>>;

  /** Alias for the reference space basis. */
  template <int dim, int range, int rank>
  using RefSpace_t = ReferenceSpaceBasis<dim, range, rank>;

  /** Type for the shared pointer of a reference space basis. */
  template <int dim, int range, int rank>
  using RefSpacePtr_t = std::shared_ptr<RefSpace_t<dim, range, rank>>;

  /** Type for the map of indices and reference space basis shared pointers. */
  template <int dim, int range, int rank>
  using MapRefSpacePtr_t = std::map<Index, RefSpacePtr_t<dim, range, rank>>;

  /** Type for the shared pointer of a grid function. */
  template <int dim, int range>
  using GridFuncPtr_t = std::shared_ptr<GridFunction<dim, range>>;

  /** Type for the map of indices and grid function shared pointers. */
  template <int dim, int range>
  using MapGridFuncPtr_t = std::map<Index, GridFuncPtr_t<dim, range>>;

  /** Type for the shared pointer of a domain. */
  template <int dim, int codim>
  using DomainPtr_t = std::shared_ptr<Domain<dim, codim>>;

  /** Type for the map of indices and domain shared pointers. */
  template <int dim, int codim>
  using MapDomainPtr_t = std::map<Index, DomainPtr_t<dim, codim>>;

  /** Alias for physical space basis. */
  template <int dim, int range, int rank, int codim>
  using PhysSpace_t = PhysicalSpaceBasis<dim, range, rank, codim>;

  /** Type for the shared pointer of a physical space basis. */
  template <int dim, int range, int rank, int codim>
  using PhysSpacePtr_t = std::shared_ptr<PhysSpace_t<dim, range, rank, codim>>;

  /** Type for the map of indices and physical space basis shared pointers. */
  template <int dim, int range, int rank, int codim>
  using MapPhysSpacePtr_t = std::map<Index, PhysSpacePtr_t<dim, range, rank, codim>>;

  /** Type for the shared pointer of a function. */
  template <int dim, int codim, int range, int rank>
  using FunctionPtr_t = std::shared_ptr<Function<dim, codim, range, rank>>;

  /** Type for the map of indices and physical space basis shared pointers. */
  template <int dim, int codim, int range, int rank>
  using MapFunctionPtr_t = std::map<Index, FunctionPtr_t<dim, codim, range, rank>>;

  typedef boost::fusion::map<
  boost::fusion::pair<Grid<0>, MapGridPtr_t<0>>,
  boost::fusion::pair<Grid<1>, MapGridPtr_t<1>>,
  boost::fusion::pair<Grid<2>, MapGridPtr_t<2>>,
  boost::fusion::pair<Grid<3>, MapGridPtr_t<3>>,
  boost::fusion::pair<SplineSpace<0,0,1>, MapSplineSpacePtr_t<0,0,1>>,
  boost::fusion::pair<SplineSpace<0,1,1>, MapSplineSpacePtr_t<0,1,1>>,
  boost::fusion::pair<SplineSpace<0,2,1>, MapSplineSpacePtr_t<0,2,1>>,
  boost::fusion::pair<SplineSpace<0,3,1>, MapSplineSpacePtr_t<0,3,1>>,
  boost::fusion::pair<SplineSpace<1,1,1>, MapSplineSpacePtr_t<1,1,1>>,
  boost::fusion::pair<SplineSpace<1,2,1>, MapSplineSpacePtr_t<1,2,1>>,
  boost::fusion::pair<SplineSpace<1,3,1>, MapSplineSpacePtr_t<1,3,1>>,
  boost::fusion::pair<SplineSpace<2,1,1>, MapSplineSpacePtr_t<2,1,1>>,
  boost::fusion::pair<SplineSpace<2,2,1>, MapSplineSpacePtr_t<2,2,1>>,
  boost::fusion::pair<SplineSpace<2,3,1>, MapSplineSpacePtr_t<2,3,1>>,
  boost::fusion::pair<SplineSpace<3,1,1>, MapSplineSpacePtr_t<3,1,1>>,
  boost::fusion::pair<SplineSpace<3,3,1>, MapSplineSpacePtr_t<3,3,1>>,
  boost::fusion::pair<RefSpace_t<0,0,1>, MapRefSpacePtr_t<0,0,1>>,
  boost::fusion::pair<RefSpace_t<0,1,1>, MapRefSpacePtr_t<0,1,1>>,
  boost::fusion::pair<RefSpace_t<0,2,1>, MapRefSpacePtr_t<0,2,1>>,
  boost::fusion::pair<RefSpace_t<0,3,1>, MapRefSpacePtr_t<0,3,1>>,
  boost::fusion::pair<RefSpace_t<1,1,1>, MapRefSpacePtr_t<1,1,1>>,
  boost::fusion::pair<RefSpace_t<1,2,1>, MapRefSpacePtr_t<1,2,1>>,
  boost::fusion::pair<RefSpace_t<1,3,1>, MapRefSpacePtr_t<1,3,1>>,
  boost::fusion::pair<RefSpace_t<2,1,1>, MapRefSpacePtr_t<2,1,1>>,
  boost::fusion::pair<RefSpace_t<2,2,1>, MapRefSpacePtr_t<2,2,1>>,
  boost::fusion::pair<RefSpace_t<2,3,1>, MapRefSpacePtr_t<2,3,1>>,
  boost::fusion::pair<RefSpace_t<3,1,1>, MapRefSpacePtr_t<3,1,1>>,
  boost::fusion::pair<RefSpace_t<3,3,1>, MapRefSpacePtr_t<3,3,1>>,
  boost::fusion::pair<GridFunction<0,0>, MapGridFuncPtr_t<0,0>>,
  boost::fusion::pair<GridFunction<0,1>, MapGridFuncPtr_t<0,1>>,
  boost::fusion::pair<GridFunction<0,2>, MapGridFuncPtr_t<0,2>>,
  boost::fusion::pair<GridFunction<0,3>, MapGridFuncPtr_t<0,3>>,
  boost::fusion::pair<GridFunction<1,1>, MapGridFuncPtr_t<1,1>>,
  boost::fusion::pair<GridFunction<1,2>, MapGridFuncPtr_t<1,2>>,
  boost::fusion::pair<GridFunction<1,3>, MapGridFuncPtr_t<1,3>>,
  boost::fusion::pair<GridFunction<2,1>, MapGridFuncPtr_t<2,1>>,
  boost::fusion::pair<GridFunction<2,2>, MapGridFuncPtr_t<2,2>>,
  boost::fusion::pair<GridFunction<2,3>, MapGridFuncPtr_t<2,3>>,
  boost::fusion::pair<GridFunction<3,1>, MapGridFuncPtr_t<3,1>>,
  boost::fusion::pair<GridFunction<3,3>, MapGridFuncPtr_t<3,3>>,
  boost::fusion::pair<Domain<0,0>, MapDomainPtr_t<0,0>>,
  boost::fusion::pair<Domain<0,1>, MapDomainPtr_t<0,1>>,
  boost::fusion::pair<Domain<0,2>, MapDomainPtr_t<0,2>>,
  boost::fusion::pair<Domain<0,3>, MapDomainPtr_t<0,3>>,
  boost::fusion::pair<Domain<1,0>, MapDomainPtr_t<1,0>>,
  boost::fusion::pair<Domain<1,1>, MapDomainPtr_t<1,1>>,
  boost::fusion::pair<Domain<1,2>, MapDomainPtr_t<1,2>>,
  boost::fusion::pair<Domain<2,0>, MapDomainPtr_t<2,0>>,
  boost::fusion::pair<Domain<2,1>, MapDomainPtr_t<2,1>>,
  boost::fusion::pair<Domain<3,0>, MapDomainPtr_t<3,0>>,
  boost::fusion::pair<PhysSpace_t<0,0,1,0>, MapPhysSpacePtr_t<0,0,1,0>>,
  boost::fusion::pair<PhysSpace_t<0,1,1,2>, MapPhysSpacePtr_t<0,1,1,2>>,
  boost::fusion::pair<PhysSpace_t<0,1,1,1>, MapPhysSpacePtr_t<0,1,1,1>>,
  boost::fusion::pair<PhysSpace_t<0,1,1,3>, MapPhysSpacePtr_t<0,1,1,3>>,
  boost::fusion::pair<PhysSpace_t<0,3,1,2>, MapPhysSpacePtr_t<0,3,1,2>>,
  boost::fusion::pair<PhysSpace_t<0,3,1,1>, MapPhysSpacePtr_t<0,3,1,1>>,
  boost::fusion::pair<PhysSpace_t<0,2,1,2>, MapPhysSpacePtr_t<0,2,1,2>>,
  boost::fusion::pair<PhysSpace_t<0,2,1,1>, MapPhysSpacePtr_t<0,2,1,1>>,
  boost::fusion::pair<PhysSpace_t<0,3,1,3>, MapPhysSpacePtr_t<0,3,1,3>>,
  boost::fusion::pair<PhysSpace_t<1,1,1,0>, MapPhysSpacePtr_t<1,1,1,0>>,
  boost::fusion::pair<PhysSpace_t<1,1,1,1>, MapPhysSpacePtr_t<1,1,1,1>>,
  boost::fusion::pair<PhysSpace_t<1,2,1,1>, MapPhysSpacePtr_t<1,2,1,1>>,
  boost::fusion::pair<PhysSpace_t<1,3,1,1>, MapPhysSpacePtr_t<1,3,1,1>>,
  boost::fusion::pair<PhysSpace_t<1,1,1,2>, MapPhysSpacePtr_t<1,1,1,2>>,
  boost::fusion::pair<PhysSpace_t<1,3,1,2>, MapPhysSpacePtr_t<1,3,1,2>>,
  boost::fusion::pair<PhysSpace_t<2,1,1,0>, MapPhysSpacePtr_t<2,1,1,0>>,
  boost::fusion::pair<PhysSpace_t<2,2,1,0>, MapPhysSpacePtr_t<2,2,1,0>>,
  boost::fusion::pair<PhysSpace_t<2,1,1,1>, MapPhysSpacePtr_t<2,1,1,1>>,
  boost::fusion::pair<PhysSpace_t<2,3,1,1>, MapPhysSpacePtr_t<2,3,1,1>>,
  boost::fusion::pair<PhysSpace_t<3,1,1,0>, MapPhysSpacePtr_t<3,1,1,0>>,
  boost::fusion::pair<PhysSpace_t<3,3,1,0>, MapPhysSpacePtr_t<3,3,1,0>>,
  boost::fusion::pair<PhysSpace_t<1,2,1,0>, MapPhysSpacePtr_t<1,2,1,0>>,
  boost::fusion::pair<PhysSpace_t<1,3,1,0>, MapPhysSpacePtr_t<1,3,1,0>>,
  boost::fusion::pair<PhysSpace_t<2,3,1,0>, MapPhysSpacePtr_t<2,3,1,0>>,
  boost::fusion::pair<PhysSpace_t<0,1,1,0>, MapPhysSpacePtr_t<0,1,1,0>>,
  boost::fusion::pair<PhysSpace_t<0,2,1,0>, MapPhysSpacePtr_t<0,2,1,0>>,
  boost::fusion::pair<PhysSpace_t<0,3,1,0>, MapPhysSpacePtr_t<0,3,1,0>>,
  boost::fusion::pair<Function<0,0,0,1>, MapFunctionPtr_t<0,0,0,1>>,
  boost::fusion::pair<Function<0,2,1,1>, MapFunctionPtr_t<0,2,1,1>>,
  boost::fusion::pair<Function<0,1,1,1>, MapFunctionPtr_t<0,1,1,1>>,
  boost::fusion::pair<Function<0,3,1,1>, MapFunctionPtr_t<0,3,1,1>>,
  boost::fusion::pair<Function<0,2,3,1>, MapFunctionPtr_t<0,2,3,1>>,
  boost::fusion::pair<Function<0,1,3,1>, MapFunctionPtr_t<0,1,3,1>>,
  boost::fusion::pair<Function<0,2,2,1>, MapFunctionPtr_t<0,2,2,1>>,
  boost::fusion::pair<Function<0,1,2,1>, MapFunctionPtr_t<0,1,2,1>>,
  boost::fusion::pair<Function<0,3,3,1>, MapFunctionPtr_t<0,3,3,1>>,
  boost::fusion::pair<Function<1,0,1,1>, MapFunctionPtr_t<1,0,1,1>>,
  boost::fusion::pair<Function<1,1,1,1>, MapFunctionPtr_t<1,1,1,1>>,
  boost::fusion::pair<Function<1,1,2,1>, MapFunctionPtr_t<1,1,2,1>>,
  boost::fusion::pair<Function<1,1,3,1>, MapFunctionPtr_t<1,1,3,1>>,
  boost::fusion::pair<Function<1,2,1,1>, MapFunctionPtr_t<1,2,1,1>>,
  boost::fusion::pair<Function<1,2,3,1>, MapFunctionPtr_t<1,2,3,1>>,
  boost::fusion::pair<Function<2,0,1,1>, MapFunctionPtr_t<2,0,1,1>>,
  boost::fusion::pair<Function<2,0,2,1>, MapFunctionPtr_t<2,0,2,1>>,
  boost::fusion::pair<Function<2,1,1,1>, MapFunctionPtr_t<2,1,1,1>>,
  boost::fusion::pair<Function<2,1,3,1>, MapFunctionPtr_t<2,1,3,1>>,
  boost::fusion::pair<Function<3,0,1,1>, MapFunctionPtr_t<3,0,1,1>>,
  boost::fusion::pair<Function<3,0,3,1>, MapFunctionPtr_t<3,0,3,1>>,
  boost::fusion::pair<Function<1,0,2,1>, MapFunctionPtr_t<1,0,2,1>>,
  boost::fusion::pair<Function<1,0,3,1>, MapFunctionPtr_t<1,0,3,1>>,
  boost::fusion::pair<Function<2,0,3,1>, MapFunctionPtr_t<2,0,3,1>>,
  boost::fusion::pair<Function<0,0,1,1>, MapFunctionPtr_t<0,0,1,1>>,
  boost::fusion::pair<Function<0,0,2,1>, MapFunctionPtr_t<0,0,2,1>>,
  boost::fusion::pair<Function<0,0,3,1>, MapFunctionPtr_t<0,0,3,1>>
  > ObjectTypes_t;

public:
  /**
   * Valid Grid instantiations.
   */
  typedef boost::fusion::vector<std::shared_ptr<Grid<0>>,
                                std::shared_ptr<Grid<1>>,
                                std::shared_ptr<Grid<2>>,
                                std::shared_ptr<Grid<3>>> ValidGridPtrs;

  /**
   * Valid Grid instantiations.
   */
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
          std::shared_ptr<SplineSpace<3,3,1>>> ValidSplineSpacePtrs;

  /**
   * Valid Domain instantiations.
   */
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
          std::shared_ptr<Domain<3,0>>> ValidDomainPtrs;

  /**
   * Valid Physical Spaces instantiations.
   */
  typedef boost::fusion::vector<
          std::shared_ptr<PhysSpace_t<0,0,1,0>>,
          std::shared_ptr<PhysSpace_t<0,1,1,2>>,
          std::shared_ptr<PhysSpace_t<0,1,1,1>>,
          std::shared_ptr<PhysSpace_t<0,1,1,3>>,
          std::shared_ptr<PhysSpace_t<0,3,1,2>>,
          std::shared_ptr<PhysSpace_t<0,3,1,1>>,
          std::shared_ptr<PhysSpace_t<0,2,1,2>>,
          std::shared_ptr<PhysSpace_t<0,2,1,1>>,
          std::shared_ptr<PhysSpace_t<0,3,1,3>>,
          std::shared_ptr<PhysSpace_t<1,1,1,0>>,
          std::shared_ptr<PhysSpace_t<1,1,1,1>>,
          std::shared_ptr<PhysSpace_t<1,2,1,1>>,
          std::shared_ptr<PhysSpace_t<1,3,1,1>>,
          std::shared_ptr<PhysSpace_t<1,1,1,2>>,
          std::shared_ptr<PhysSpace_t<1,3,1,2>>,
          std::shared_ptr<PhysSpace_t<2,1,1,0>>,
          std::shared_ptr<PhysSpace_t<2,2,1,0>>,
          std::shared_ptr<PhysSpace_t<2,1,1,1>>,
          std::shared_ptr<PhysSpace_t<2,3,1,1>>,
          std::shared_ptr<PhysSpace_t<3,1,1,0>>,
          std::shared_ptr<PhysSpace_t<3,3,1,0>>,
          std::shared_ptr<PhysSpace_t<1,2,1,0>>,
          std::shared_ptr<PhysSpace_t<1,3,1,0>>,
          std::shared_ptr<PhysSpace_t<2,3,1,0>>,
          std::shared_ptr<PhysSpace_t<0,1,1,0>>,
          std::shared_ptr<PhysSpace_t<0,2,1,0>>,
          std::shared_ptr<PhysSpace_t<0,3,1,0>>> ValidPhysSpacePtrs;

  /**
   * Valid Function instantiations.
   */
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
          std::shared_ptr<Function<0,0,3,1>>> ValidFunctionPtrs;


private:


  /** @name Constructors*/
  ///@{
  /**
   * To document
   */
  ObjectsContainer();

  /**
   * Copy constructor.
   */
  ObjectsContainer(const self_t &container) = default;

  /**  Move constructor */
  ObjectsContainer(self_t &&container) = default;

public:
  /** Destructor */
  ~ObjectsContainer() = default;
  ///@}

public:
  /**
   * @name Creators
   * @note The functions here return:
   * - a <b>non-const</b> Container object wrapped by a std::shared_ptr
   * in the case of <b>create()</b> functions;
   * - a <b>const</b> Container object wrapped by a std::shared_ptr
   * in the case of <b>const_create()</b> functions.
   */
  ///@{
  /**
   * Creates a objects container (non-const).
   */
  static std::shared_ptr<self_t> create();

  /**
   * Creates a objects container (const).
   */
  static std::shared_ptr<const self_t> const_create();
  ///@}

private:
  /**
   * @name Assignment operators
   */
  ///@{

  /**
   * Copy assignment operator.
   */
  self_t &operator=(const self_t &container) = default;

  /**
   * Move assignment operator.
   */
  self_t &operator=(self_t &&container) = default;
  ///@}

public:
  /**
   * @todo To be documented.
   */
  template <int dim>
  void insert_grid (const GridPtr_t<dim> grid, const Index &id);

  /**
   * @todo To be documented.
   */
  template <int dim>
  GridPtr_t<dim> get_grid (const Index &id) const;

  /**
   * @todo To be documented.
   */
  template <int dim>
  bool is_grid (const Index &id) const;

  /**
   * @todo To be documented.
   */
  template <int dim, int range, int rank>
  void insert_spline_space (const SplineSpacePtr_t<dim, range, rank> spline_space,
                         const Index &id);

  /**
   * @todo To be documented.
   */
  template <int dim, int range, int rank>
  SplineSpacePtr_t<dim, range, rank> get_spline_space (const Index &id) const;

  /**
   * @todo To be documented.
   */
  template <int dim, int range, int rank>
  bool is_spline_space (const Index &id) const;

  /**
   * @todo To be documented.
   */
  template <int dim, int range, int rank>
  void insert_ref_space (const RefSpacePtr_t<dim, range, rank> ref_space,
                      const Index &id);

  /**
   * @todo To be documented.
   */
  template <int dim, int range, int rank>
  RefSpacePtr_t<dim, range, rank> get_ref_space (const Index &id) const;

  /**
   * @todo To be documented.
   */
  template <int dim, int range, int rank>
  bool is_ref_space (const Index &id) const;

  /**
   * @todo To be documented.
   */
  template <int dim, int range, int rank>
  std::shared_ptr<BSpline<dim, range, rank>> get_bspline (const Index &id) const;

  /**
   * @todo To be documented.
   */
  template <int dim, int range, int rank>
  std::shared_ptr<NURBS<dim, range, rank>> get_nurbs (const Index &id) const;

  /**
   * @todo To be documented.
   */
  template <int dim, int range, int rank>
  bool is_bspline (const Index &id) const;

  /**
   * @todo To be documented.
   */
  template <int dim, int range, int rank>
  bool is_nurbs (const Index &id) const;

  /**
   * @todo To be documented.
   */
  template <int dim, int range>
  void insert_grid_function (const GridFuncPtr_t<dim, range> grid_func,
                          const Index &id);

  /**
   * @todo To be documented.
   */
  template <int dim, int range>
  GridFuncPtr_t<dim, range> get_grid_function (const Index &id) const;

  /**
   * @todo To be documented.
   */
  template <int dim, int range>
  bool is_grid_function (const Index &id) const;

  /**
   * @todo To be documented.
   */
  template <int dim, int codim>
  void insert_domain (const DomainPtr_t<dim, codim> domain,
                   const Index &id);

  /**
   * @todo To be documented.
   */
  template <int dim, int codim>
  DomainPtr_t<dim, codim> get_domain (const Index &id) const;

  /**
   * @todo To be documented.
   */
  template <int dim, int codim>
  bool is_domain (const Index &id) const;

  /**
   * @todo To be documented.
   */
  template <int dim, int range, int rank, int codim>
  void insert_phys_space_basis (const PhysSpacePtr_t<dim, range, rank, codim> space,
                             const Index &id);

  /**
   * @todo To be documented.
   */
  template <int dim, int range, int rank, int codim>
  PhysSpacePtr_t<dim, range, rank, codim> get_phys_space_basis (const Index &id) const;

  /**
   * @todo To be documented.
   */
  template <int dim, int range, int rank, int codim>
  bool is_phys_space_basis (const Index &id) const;

  /**
   * @todo To be documented.
   */
  template <int dim, int codim, int range, int rank>
  void insert_function (const FunctionPtr_t<dim, codim, range, rank> function,
                     const Index &id);

  /**
   * @todo To be documented.
   */
  template <int dim, int codim, int range, int rank>
  FunctionPtr_t<dim, codim, range, rank> get_function (const Index &id) const;

  /**
   * @todo To be documented.
   */
  template <int dim, int codim, int range, int rank>
  bool is_function (const Index &id) const;

  /**
   * @todo To be documented.
   */
  template <class T>
  void insert_object (const std::shared_ptr<T> object, const Index &id);

  /**
   * @todo To be documented.
   */
  template <class T>
  std::shared_ptr<T> get_object (const Index &id) const;

  /**
   * @todo To be documented.
   */
  template <class T>
  bool is_object (const Index &id) const;

  bool is_id_present (const Index &id) const;

private:
  /**
   * Container for the objects.
   */
  ObjectTypes_t objects_;

};

IGA_NAMESPACE_CLOSE

#endif /*OBJECTS_CONTAINER_H_ */
