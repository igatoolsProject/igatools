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

#ifndef __GRID_ELEMENT_H_
#define __GRID_ELEMENT_H_

#include <igatools/base/config.h>
#include <igatools/utils/tensor_range.h>
#include <igatools/base/quadrature.h>
#include <igatools/geometry/grid.h>
#include <igatools/geometry/grid_handler.h>
#include <igatools/utils/value_vector.h>
#include <igatools/basis_functions/values_cache.h>

IGA_NAMESPACE_OPEN

/**
 * @brief This is the base class for all kind of elements.
 *
 * It a pure virtual class and contains the interfaces for the functions
 * that the derived classes must implement
 * (basically the advance/increment operator <tt>++</tt>).
 * Moreover it contains the element-property identifier and the associated getter.
 *
 * @ingroup elements
 *
 * @author martinelli, 2105
 *
 */
class Element
{
public:

  Element(const PropId &property);

  virtual ~Element() = default;

  virtual void operator++() = 0;


  /**
   * @name Functions for managing/querying the element properties.
   */
  ///@{
  /**
   * Tests if a certain element @p property is TRUE.
   */
  bool has_property(const PropId &property) const;

  /**
   * Returns the property of the element.
   */
  const PropId &get_property() const;
  ///@}

private:

  PropId property_;

};


/**
 * @brief Element accessor for the Grid.
 *
 *
 *
 * See module (and the submodules) on \ref elements for a general overview.
 * @ingroup elements
 *
 * ### Quantities handled by the cache
 * - grid_element::_Point i.e. evaluation points mapped in a element of the parametric domain
 * - grid_element::_Weight i.e. quadrature weights associated to each evaluation point,
 * multiplied by the element <tt>dim</tt>-dimensional measure.
 *
 * @author pauletti, 2012, 2013, 2014, 2015
 * @author martinelli, 2012, 2013, 2014, 2105
 *
 */
template <int dim>
class GridElement : public Element
{
private:
  using self_t = GridElement<dim>;

public:
  /** Type required by the GridIterator templated iterator */
  using ContainerType = const Grid<dim>;
  using IndexType = typename ContainerType::IndexType;
  using List = typename ContainerType::List;
  using ListIt = typename ContainerType::ListIt;

  using Point = Points<dim>;

  using CacheFlags = grid_element::CacheFlags;
  using Flags = grid_element::Flags;

  /** @name Constructors */
  ///@{
protected:

  /**
   * Default constructor. Not allowed to be used.
   */
  GridElement() = delete;

public:
  /**
   * Construct an accessor pointing to the element with
   * flat index @p elem_index of the Grid @p grid.
   */
  GridElement(const std::shared_ptr<const Grid<dim>> &grid,
              const ListIt &index,
              const PropId &prop = ElementProperties::active);

  /**
   * Copy constructor. Not allowed to be used.
   */
  GridElement(const self_t &elem) = delete;

  /**
   * Move constructor.
   */
  GridElement(self_t &&elem) = default;

  /**
   * Destructor.
   */
  ~GridElement() = default;
  ///@}


  /** @name Assignment operators */
  ///@{
  /**
   * Copy assignment operator. Not allowed to be used.
   */
  self_t &operator=(const self_t &element) = delete;

  /**
   * Move assignment operator.
   */
  self_t &operator=(self_t &&elem) = default;
  ///@}





  /**
   * @name Functions/operators for moving the element in the Grid.
   *
   * @note They should be called only by the GridIterator.
   */
  ///@{
  void operator++() override;

  /**
   * Move the element to the one specified by <tt>elem_id</tt>.
   *
   * In Debug mode an assertion will be raised
   * if the GridElement specified by <tt>elem_id</tt> has not the same property of the
   * calling GridElement.
   *
   * @warning Use this function only if you know what you are doing
   */
  void move_to(const IndexType &elem_id);
  ///@}

  /**
   * @name Comparison operators
   * @note In order to be meaningful, the comparison must be performed on elements defined on
   * the <b>same Grid</b>
   * (in the sense that the pointer to the grid held by the elements involved in the comparison
   * must point to the same Grid object).
   */
  ///@{
  /**
   * True if the elements have the same index.
   *  @note In debug mode, it is also check they both refer to
   *  the same cartesian grid. No check is done on the cache.
   */
  bool operator==(const self_t &elem) const;

  /**
   * True if the elements have different index.
   *  @note In debug mode, it is also check they both refer to
   *  the same cartesian grid. No check is done on the cache.
   */
  bool operator!=(const self_t &elem) const;

  /**
   * True if the flat-index of the element on the left is smaller than
   * the flat-index of the element on the right.
   *  @note In debug mode, it is also check they both refer to
   *  the same cartesian grid. No check is done on the cache.
   */
  bool operator<(const self_t &elem) const;

  /**
   * True if the flat-index of the element on the left is bigger than
   * the flat-index of the element on the right.
   *  @note In debug mode, it is also check they both refer to
   *  the same cartesian grid. No check is done on the cache.
   */
  bool operator>(const self_t &elem) const;
  ///@}

  ///@name Query information that does not require the use of the cache
  ///@{

  const IndexType &get_index() const;

  /** Return the Grid from which the element belongs.*/
  std::shared_ptr<const Grid<dim>> get_grid() const;


  /**
     * Test if the element has a boundary face.
     */
  template<int k = (dim > 0) ? (dim-1) : 0 >
  bool is_boundary() const;

  /**
   * Test if the face identified by @p face_id on the current element is on the
   * boundary of the cartesian grid.
   */
  template<int k = (dim > 0) ? (dim-1) : 0>
  bool is_boundary(const Index sub_elem_id) const;

  /**
   * Return the @p i-th vertex
   */
  Point vertex(const int i) const;


  /**
   * Returns the lengths of the coordinate sides of the <tt>sdim</tt>-dimensional element
   * of the Grid.
   * For example in 3 dimensions
   * \code{.cpp}
     const int s_id = 1;
     auto length = elem.template get_side_lengths<2>(s_id);
     // length[0] is the length of along the first active direction of the face 1
     // length[1] is the length of along the second active direction of the face 1
     \endcode
   */
  template<int sdim>
  const Points<sdim> get_side_lengths(const int s_id) const;

  /**
   * Prints internal information about the GridElementAccessor.
   * Its main use is for testing and debugging.
   */
  void print_info(LogStream &out) const;

  /**
   * Returns true if two elements belongs from the same Grid.
   */
  bool same_grid_of(const self_t &elem) const;

  ///@}

  ///@name Query information that requires the use of the cache
  ///@{
  template <int sdim>
  Real get_measure(const int s_id) const;


  /**
   * Returns the quadrature weights corresponding to the <tt>sdim</tt>
   * dimensional s_id-th sub-element.
   *
   * @note The returned weights are the quadrature unit weights multiplied by the
   * <tt>sdim>-dimensional element measure.
   */
  template <int sdim>
  ValueVector<Real> get_weights(const int s_id) const;

  /**
   * Returns the quadrature points corresponding to the <tt>sdim</tt>
   * dimensional s_id-th sub-element.
   */
  template <int sdim>
  ValueVector<Point> get_points(const int s_id = 0) const;

  ValueVector<Point> get_element_points() const;


  /**
   * Returns the unitary quadrature scheme corresponding to the <tt>sdim</tt>-dimensional
   * s_id-th sub-element.
   */
  template <int sdim>
  std::shared_ptr<const Quadrature<sdim>> get_quad() const
  {
    return quad_list_.template get_quad<sdim>();
  }


  void print_cache_info(LogStream &out) const;

  ///@}


private:
  template <class Accessor> friend class GridIteratorBase;
  friend class GridHandler<dim>;

protected:

  /** Cartesian grid from which the element belongs.*/
  std::shared_ptr<const Grid<dim>> grid_;

private:

  /** Index in the property list of the current element */
  typename List::iterator index_it_;

  /**
   * @name Types and data needed for the definition/use of the cache.
   */
  ///@{

public:


  using _Point = grid_element::_Point;
  using _Weight = grid_element::_Weight;

private:
  using CType = boost::fusion::map<
                boost::fusion::pair< _Point,DataWithFlagStatus<ValueVector<Points<dim>>>>,
                boost::fusion::pair<_Weight,DataWithFlagStatus<ValueVector<Real>>>
                >;


  using ValuesCache = PointValuesCache<dim,CType>;

  using CacheType = AllSubElementsCache<ValuesCache>;

  /** List of quadrature pointers for all sub elements */
  QuadList<dim> quad_list_;

  /** The local cache. */
  CacheType all_sub_elems_cache_;


  template <class ValueType, int sdim>
  const auto &
  get_values_from_cache(const int s_id) const
  {
    const auto &cache = all_sub_elems_cache_.template get_sub_elem_cache<sdim>(s_id);
    return cache.template get_data<ValueType>();
  }
  ///@}


};



IGA_NAMESPACE_CLOSE

#endif
