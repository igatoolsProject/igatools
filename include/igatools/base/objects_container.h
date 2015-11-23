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

template <int dim> class Grid;
template <int dim, int range, int rank> class SplineSpace;
template <int dim, int range, int rank> class ReferenceSpaceBasis;

IGA_NAMESPACE_OPEN

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

  typedef boost::fusion::map<
  boost::fusion::pair<Grid<0>, MapGridPtr_t<0>>,
  boost::fusion::pair<Grid<1>, MapGridPtr_t<1>>,
  boost::fusion::pair<Grid<2>, MapGridPtr_t<2>>,
  boost::fusion::pair<Grid<3>, MapGridPtr_t<3>>,
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
  boost::fusion::pair<RefSpace_t<3,3,1>, MapRefSpacePtr_t<3,3,1>>> ObjectTypes_t;


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
  void add_grid (const GridPtr_t<dim> grid, const Index &id);

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
  void add_spline_space (const SplineSpacePtr_t<dim, range, rank> spline_space,
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
  void add_ref_space (const RefSpacePtr_t<dim, range, rank> ref_space,
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

private:
  /**
   * Container for the objects.
   */
  ObjectTypes_t objects_;

};

IGA_NAMESPACE_CLOSE

#endif /*OBJECTS_CONTAINER_H_ */
