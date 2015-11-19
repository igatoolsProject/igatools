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

#ifndef __CARTESIAN_GRID_H_
#define __CARTESIAN_GRID_H_

#include <igatools/base/config.h>
#include <igatools/base/tensor.h>
#include <igatools/base/logstream.h>
#include <igatools/utils/safe_stl_array.h>
#include <igatools/utils/safe_stl_map.h>
#include <igatools/utils/cartesian_product_array.h>
#include <igatools/utils/dynamic_multi_array.h>
#include <igatools/utils/element_index.h>
#include <igatools/geometry/unit_element.h>
#include <igatools/base/array_utils.h>
#include <igatools/geometry/bbox.h>
#include <igatools/geometry/grid_iterator.h>
#include <igatools/base/properties_id_container.h>

#include<memory>
#ifdef MESH_REFINEMENT
#include <boost/signals2.hpp>
#endif

IGA_NAMESPACE_OPEN

template <int> class Grid;
template <int> class GridElement;
template <int> class GridHandler;


/**
 * @brief Grid in <tt>dim</tt>-dimensional space with cartesian-product structure.
 *
 * The data defining the grid is a <tt>dim</tt>-dimensional array of coordinates,
 * representing the positions of the grid wires.
 *
 * Two consecutive coordinates (along the same direction) defines an <em>interval</em>.
 *
 * Then, the tensor-product of the intervals along the coordinate directions define
 * the elements that are tiling the domain covered by the Grid;
 *
 * The element type for the Grid is GridElement.
 *
 * The elements can be iterated with the GridIterator object.
 * The iterators for traversing the the entire set of elements of the grid
 * are build with the begin(), end(), cbegin(), cend() functions called with no
 * arguments.
 *
 * ### Element properties
 * The elements can have associated a certain list of <em>element properties</em>
 * (identified by one PropId).
 * There is no pre-defined list of element properties: any property can be
 * defined and added to the list
 * of properties stored within the Grid. This choice is made because
 * the properties are usually specific for a given problem (e.g. <em>active</em>
 * elements in hierarchical grids or <em>marked</em> elements for a-posteriori
 * error estimators).
 *
 * For example, let suppose we want to build a bi-dimensional grid (over the
 * unit square \f$ [0,1]\times[0,1] \f$) made of 4 equally-sized intervals in
 * each coordinate direction, and then assign the property <tt>"active"</tt> to
 * and the <tt>"marked"</tt> in an alternated (i.e. checkerboard) manner.
 * @code{.cpp}
   auto grid = Grid<2>::create(5); // here we create a 4x4 grid

   const PropId property_active = "active"; // here we choose the string identifying the first property
   const PropId property_marked = "marked"; // here we choose the string identifying the second property

   grid->add_elements_property(property_active); // here we add the first property to the elements property database managed by the grid
   grid->add_elements_property(property_marked); // here we add the second property to the elements property database managed by the grid

   //the next loop set the elements with the right property
   int flat_id = 0;
   for (const auto elem : *grid)
   {
       const auto elem_id = elem.get_index();
       if (flat_id % 2 == 0)
           grid->set_property_status_elem(property_active,elem_id,true);
       else
           grid->set_property_status_elem(property_marked,elem_id,true);

       ++flat_id;
   }
   @endcode
 *
 * The list of elements with a given property can be obtained from the Grid
 * using the function get_elements_with_property(const PropId &property).
 * For example if we want the IDs of the elements having the property <tt>"marked"</tt> we can write
 * @code{.cpp}
   const auto & elems_id_marked = grid->get_elements_with_property("marked");
   // now elems_id_marked is a const reference to a std::set<Index> containing the values [1 3 5 7 9 11 13 15]
   @endcode
 *
 * The elements with the same property can be iterated with the GridIterator object.
 * The iterators for elements with a given property can be obtained with the
 * begin(), end(), cbegin(), cend() functions invoked
 * with the property name as input argument. For example, to iterate over the elements that share the property
 * <tt>"active"</tt> (provided that property <tt>"active"</tt> is defined for some elements in the grid)
 * you can use
 * @code{.cpp}
   PropId property_active = "active";
   auto elem_active = grid->begin(property_active);
   auto  end_active = grid->end(property_active);
   for ( ; elem_active != end_active ; ++elem_active)
   {
   // do something with the active elements
   }
   @endcode
 *
 * In Debug mode, an assertion will be raised if the the property string used as input argument
 * in the begin() / end() functions is not defined in the Grid.
 *
 * ### Getting a Grid by an XML structure.
 * A Grid object can be obtained by a Boost XML structure using the
 * function
 * get_grid_from_xml().
 * An example of valid XML structure for a 2-dimensional Grid is given
 * by the following
 * XML block
 * @code{.xml}
   <Grid Dim="2">
     <Knots Direction="0" Size="2">
       0.0 1.0
     </Knots>
     <Knots Direction="1" Size="2">
        0.0 1.0
     </Knots>
   </Grid>
 * @endcode
 * where the <tt>Dim</tt> attribute in the <tt>"Grid"</tt> tag is the
 * dimensionality of the Grid
 * and the <tt>Size</tt> attribute in each <tt>"Knots"</tt> tag is the number
 *  of unique knots values (also called
 * "breakpoints") along each coordinate direction.
 *
 * @author M. Martinelli 2012, 2013, 2014, 2015
 * @author pauletti 2012, 2013, 2014, 2015
 *
 * @see get_grid_from_xml()
 * @ingroup h_refinement
 * @ingroup containers
 * @ingroup serializable
 * @todo document more
 */
template<int dim_>
class Grid :
  public std::enable_shared_from_this<Grid<dim_>>
{
private:
  /** Type for current class. */
  using self_t = Grid<dim_>;
  using base_t = TensorSizedContainer<dim_>;

public:
  static const int dim = dim_;

  /** Type for the element accessor. */
  using ElementAccessor = GridElement<dim_>;

  /** Type for the iterator over the elements of the grid (const version).  */
  using ElementIterator = GridIterator<ElementAccessor>;


  using ElementHandler = GridHandler<dim_>;


  using Point = Points<dim_>;

  using IndexType = ElementIndex<dim_>;

  using PropertyList = PropertiesIdContainer<IndexType>;
  using List = typename PropertyList::List;
  using ListIt = typename PropertyList::List::iterator;


  /** Type for the vector of knot vectors */
  using KnotCoordinates = CartesianProductArray<Real, dim_>;

  /** @name Constructors*/
  ///@{
public:
  /**
   * Construct a uniform cartesian grid of the unit <tt>dim</tt>-dimensional
   * hypercube \f$[0,1]^{dim}\f$, with @p n knots (equally spaced) in each dimension.
   *
   * @note This constructor is public because it can be used as default constructor and
   * for the serialization mechanism, the defult constructor must be public.
   */
  explicit Grid(const Size n = 2);

protected:

  /**
   * Construct a uniform cartesian grid of the unit <tt>dim</tt>-dimensional
   * hypercube \f$[0,1]^{dim}\f$,
   * with <tt>n[0],..,n[dim-1</tt>] knots in each dimension
   * respectively.
   */
  explicit Grid(const TensorSize<dim_> &n_knots);

  /**
   * @todo Document me
   */
  explicit Grid(const BBox<dim_> &bbox,
                const Size n_knots);

  /**
   * @todo Document me
   */
  explicit Grid(const BBox<dim_> &bbox,
                const TensorSize<dim_> &n_knots);
  /**
   * Construct a cartesian grid where the knot coordinate in each
   * direction is provided as CartesianProductArray object.
   *
   * The knot coordinate in each direction must be sorted and without
   * repetition.
   * @note In Debug mode, a check for this precondition (up to machine precision)
   * is perform and if not satistified an exception is raised.
   */
  explicit
  Grid(const KnotCoordinates &knots);

  /**
   * Construct a cartesian grid where the knot coordinate in each
   * direction is provided as SafeSTLArray of SafeSTLVector<Real>.
   *
   * The knot coordinate in each direction must be sorted and without
   * repetition.
   * @note In Debug mode, a check for this precondition (up to machine precision)
   * is perform and if not satisfied an exception is raised.
   */
  explicit
  Grid(const SafeSTLArray<SafeSTLVector<Real>,dim_> &knot_coords);

public:
  /**
   * Copy constructor.
   * Perform a deep copy of the member variables except the
   * signal_refine_ variable, that is not copied at all.
   */
  Grid(const self_t &grid);

  /**  Move constructor */
  Grid(self_t &&grid) = default;

  /** Destructor */
  ~Grid() = default;
  ///@}

  /**
   * @name Creators
   * @note The functions here return:
   * - a <b>non-const</b> Grid object wrapped by a std::shared_ptr
   * in the case of <b>create()</b> functions;
   * - a <b>const</b> Grid object wrapped by a std::shared_ptr
   * in the case of <b>const_create()</b> functions.
   */
  ///@{
  /**
   * Creates a uniform (non-const) cartesian grid of the unit <tt>dim</tt>-dimensional
   * hypercube \f$[0,1]^{dim}\f$, with @p n knots (equally spaced) in each
   * dimension.
   */
  static std::shared_ptr<self_t> create(const Index n_knots = 2);

  /**
   * Creates a uniform (const) cartesian grid of the unit <tt>dim</tt>-dimensional
   * hypercube \f$[0,1]^{dim}\f$, with @p n knots (equally spaced) in each
   * dimension.
   */
  static std::shared_ptr<const self_t> const_create(const Index n_knots = 2);

  /**
   * Creates a uniform (non-const) cartesian grid of the unit <tt>dim</tt>-dimensional
   * hypercube \f$[0,1]^{dim}\f$,
   * with <tt>n_knots[0],..,n_knots[dim-1</tt>] knots in each dimension
   * respectively.
   */
  static std::shared_ptr<self_t>
  create(const TensorSize<dim_> &n_knots);

  /**
   * Creates a uniform (const) cartesian grid of the unit <tt>dim</tt>-dimensional
   * hypercube \f$[0,1]^{dim}\f$,
   * with <tt>n_knots[0],..,n_knots[dim-1</tt>] knots in each dimension
   * respectively.
   */
  static std::shared_ptr<const self_t>
  const_create(const TensorSize<dim_> &n_knots);

  static std::shared_ptr<self_t>
  create(const BBox<dim_> &bbox, const Size n_knots);

  static std::shared_ptr<const self_t>
  const_create(const BBox<dim_> &bbox, const Size n_knots);

  /**
   * Construct a (non-const) cartesian grid where the knot coordinate in each
   * direction is provided as CartesianProductArray object.
   *
   * @pre The knot coordinate in each direction must be sorted and without
   * repetition.
   *
   * @note In Debug mode, a check for this precondition
   * (up to machine precision)
   * is perform and if not satisfied an exception is raised.
   */
  static std::shared_ptr<self_t>
  create(const KnotCoordinates &knots);

  /**
   * Construct a (const) cartesian grid where the knot coordinate in each
   * direction is provided as CartesianProductArray object.
   *
   * @pre The knot coordinate in each direction must be sorted and without
   * repetition.
   *
   * @note In Debug mode, a check for this precondition
   * (up to machine precision)
   * is perform and if not satisfied an exception is raised.
   */
  static std::shared_ptr<const self_t>
  const_create(const KnotCoordinates &knots);

  /**
   * Construct a (non-const) cartesian grid where the knot coordinate in each
   * direction is provided as SafeSTLArray of SafeSTLVector<Real>.
   *
   * @pre The knot coordinate in each direction must be sorted and without
   * repetition.
   *
   * @note In Debug mode, a check for this precondition
   * (up to machine precision)
   * is perform and if not satisfied an exception is raised.
   */
  static std::shared_ptr<self_t>
  create(const SafeSTLArray<SafeSTLVector<Real>,dim_> &knots);

  /**
   * Construct a (const) cartesian grid where the knot coordinate in each
   * direction is provided as SafeSTLArray of SafeSTLVector<Real>.
   *
   * @pre The knot coordinate in each direction must be sorted and without
   * repetition.
   *
   * @note In Debug mode, a check for this precondition
   * (up to machine precision)
   * is perform and if not satisfied an exception is raised.
   */
  static std::shared_ptr<const self_t>
  const_create(const SafeSTLArray<SafeSTLVector<Real>,dim_> &knots);

  /**
   * @todo document me
   */
  static std::shared_ptr<self_t>
  create(const BBox<dim_> &bbox, const TensorSize<dim_> &n_knots);

  /**
   * @todo document me
   */
  static std::shared_ptr<const self_t>
  const_create(const BBox<dim_> &bbox, const TensorSize<dim_> &n_knots);

  /**
   * Creates a (non-const) cartesian grid object (wrapped by a shared_ptr) using
   * the copy constructor.
   *
   * @note A grid built in this way is totally uncoupled from the grid used as argument
   * of this function. For example, a refinement of a grid does not affect the other gird.
   */
  static std::shared_ptr<self_t> create(const self_t &grid);

  /**
   * Creates a (const) cartesian grid object (wrapped by a shared_ptr) using
   * the copy constructor.
   *
   * @note A grid built in this way is totally uncoupled from the grid used as argument
   * of this function. For example, a refinement of a grid does not affect the other gird.
   */
  static std::shared_ptr<const self_t> const_create(const self_t &grid);
  ///@}

  /**
   * @name Assignment operators
   */
  ///@{

  /**
   * Copy assignment operator.
   */
  self_t &operator=(const self_t &grid) = default;

  /**
   * Move assignment operator.
   */
  self_t &operator=(self_t &&grid) = default;
  ///@}

  /**
   * @name Getting grid information
   */
  ///@{
  /**
   * Returns the number of elements with the <tt>property</tt> specified
   * as input argument.
   */
  Size get_num_elements(const PropId &prop = ElementProperties::active) const;

  /** Total number of elements, including active and non-active */
  Size get_num_all_elems() const;

  /**
   * Total number of one dimensional intervals along each
   * coordinate direction.
   */
  TensorSize<dim_> get_num_intervals() const;

  /**
   * Query the number of knot values along each coordinate direction
   * represented by the Grid.
   */
  TensorSize<dim_> get_num_knots_dim() const;

  /**
   * Returns the knot coordinates along the direction @p i.
   */
  SafeSTLVector<Real> const &get_knot_coordinates(const int i) const;

  /**
   * Returns the knot coordinates along all the directions (const version).
   */
  KnotCoordinates const &get_knot_coordinates() const;

  /**
   * Computes the interval lengths along each direction.
   */
  CartesianProductArray<Real, dim_> get_element_lengths() const;

  /**
   * Returns the smallest <tt>dim_</tt>-dimensional bounding box enclosing the
   * domain represented by the Grid object.
   */
  BBox<dim_> get_bounding_box() const;
  ///@}


  std::unique_ptr<ElementHandler> create_cache_handler() const;

  /**
   * @name Iterating of grid elements
   */
  ///@{
  /**
   * This function returns a element (const) iterator to the first element of the patch.
   */
  ElementIterator begin(const PropId &prop = ElementProperties::active) const;

  /**
   * This function returns a element (const) iterator to one-pass the end of patch.
   */
  ElementIterator end(const PropId &prop = ElementProperties::active) const;

  /**
   * This function returns a element (const) iterator to the first element of the patch.
   */
  ElementIterator cbegin(const PropId &prop = ElementProperties::active) const;

  /**
   * This function returns a element (const) iterator to one-pass the end of patch.
   */
  ElementIterator cend(const PropId &prop = ElementProperties::active) const;
  ///@}

  /**
   * @name Dealing with boundary information
   */
  ///@{
  /**
   * Get the patch @p face boundary id.
   */
  boundary_id get_boundary_id(const int face) const;

  /**
   * Set the patch @p face to have the boundary @p id.
   */
  void set_boundary_id(const int face, const boundary_id id);

  template<int sdim>
  using BoundaryNormal = SafeSTLArray<Points<dim_>, dim_-sdim>;

  /**
   * Returns the outward pointing
   * unit normal vector space to the element of sub dim_ k.
   */
  template<int sdim>
  BoundaryNormal<sdim> get_boundary_normals(const int s_id) const;

  template<int sdim>
  using SubGridMap =
    SafeSTLMap<typename Grid<sdim>::IndexType, IndexType>;

  /**
   * Construct a sub grid of dimension k conforming to
   * the grid sub element sub_elem_id and a map from the elements of
   * the sub_element grid to the corresponding element of the current
   * grid.
   */
  template<int sdim>
  std::shared_ptr<Grid<sdim> >
  get_sub_grid(const int sub_elem_id, SubGridMap<sdim> &elem_map) const;
  ///@}

  /**
   * \brief Returns the TensorIndex corresponding to the element with flat index @p elem_flat_id.
   */
  TensorIndex<dim_> flat_to_tensor_element_id(const int elem_flat_id) const;

  /**
   * \brief Returns the flat index corresponding to the element with TensorIndex @p elem_tensor_id.
   *
   * @note In the case of Grid<0>, this function returns 0 despite the input argument @p elem_tensor_id.
   */
  int tensor_to_flat_element_id(const TensorIndex<dim_> &elem_tensor_id) const;

  /**
   * Given a vector of points, this function return a map with
   * in which the keys are the element id and the associated value is
   * a list containing the point id lying on the associated element.
   *
   * @warning If the point is lying exactly on knot line(s),
   * then the point can have intersection with multiple elements, but the function
   * returns only the first element that owns the point.
   */
  std::map<IndexType, SafeSTLVector<int> >
  find_elements_id_of_points(const ValueVector<Points<dim_>> &points) const;

#if 0
  /**
   * Given a vector of points, this function return the id of the
   * elements that intersects with the point.
   *
   * @note If the point is lying exactly on knot line(s),
   * then the point can have intersection with multiple elements.
   */
  SafeSTLVector<IndexType>
  find_elements_id_of_points(const ValueVector<Points<dim_>> &points) const;
#endif

  /**
   * Returns TRUE if the point is exactly on an internal knots line (at least).
   */
  bool test_if_point_on_internal_knots_line(const Points<dim_> &point) const;


  /**
   * Prints debug information of the Grid to a LogStream.
   */
  void print_info(LogStream &out) const;

  /**
   * Comparison operator. Returns true if the knot coordinates of two grid
   * are equal.
   */
  bool operator==(const Grid<dim_> &grid) const;


  /**
    * This function returns TRUE if the current grid object is a <em>refinement</em> of
    *  @p grid_to_compare_with,
    * i.e. if the knots of the current grid object are present in @p grid_to_compare_with.
    * @note The functions returns TRUE also if the knots in the current grid object
    * are equal to the knots in @p grid_to_compare_with.
    */
  bool same_knots_or_refinement_of(const Grid<dim_> &grid_to_compare_with) const;

#ifdef MESH_REFINEMENT

private:

  /** Type for the insert_knots signal. */
  using signal_insert_knots_t =
    boost::signals2::signal<
    void (const SafeSTLArray<SafeSTLVector<Real>,dim_> &new_knots,const Grid<dim_> &old_grid)>;

public:

  /** Slot type for the insert_knots signal. */
  using SignalInsertKnotsSlot = typename signal_insert_knots_t::slot_type;

  /** @name Functions for performing grid refinement */
  ///@{
  /**
   * Perform a uniform refinement of the grid along the @p direction_id
   * direction,
   * dividing each interval into @p n_subdivisions intervals.
   * @param[in] direction_id Direction along which the refinement is
   * performed.
   * @param[in] n_subdivisions Number of subdivision in which each interval
   * in the grid
   * along the specified direction is divided. This value must be >= 2.
   */
  void refine_direction(const int direction_id, const Size n_subdivisions);

  /**
   * Refine the cartesian grid and the objects connected to it (if any),
   * e.g. maps, spaces, etc.
   * along the directions specified by the true values in the entries of the
   * array of bools @p refinement_directions,
   * and with a number of subdivisions for each interval
   * (along each direction)
   * specified by @p n_subdivisions.
   *
   * @note If the i-th direction is not active for the refinement
   * (i.e. <tt>refinement_directions[i] == false</tt>),
   *  then the corresponding value <tt>n_subdivisions[i]</tt> will be ignored.
   */
  void refine_directions(
    const SafeSTLArray<bool,dim_> &refinement_directions,
    const SafeSTLArray<Size,dim_> &n_subdivisions);

  /**
   * Refine the cartesian grid and the objects connected to it (if any),
   * e.g. maps, spaces, etc.
   *
   * Each interval in the unrefined grid is uniformly divided
   * in @p n_subdivisions sub-intervals.
   */
  void refine(const Size n_subdivisions = 2);


  /**
   *  Connect a slot (i.e. a function pointer) to the refinement signals
   *  which will be
   *  emitted whenever a insert_knots() function is called by an object holding
   *  a Grid member.
   */
  boost::signals2::connection
  connect_insert_knots(const SignalInsertKnotsSlot &subscriber);


  /**
   * Insert the @p knots_to_insert to the grid and to the object that are using the grid.
   * @note The @p knots_to_insert may contain multiple knot values in each direction.
   */
  void insert_knots(SafeSTLArray<SafeSTLVector<Real>,dim_> &knots_to_insert);

  /**
   * Returns the grid before the last refinement. If no refinement is
   * performed, this function returns a null pointer.
   */
  std::shared_ptr<const self_t > get_grid_pre_refinement() const;
  ///@}
#endif // MESH_REFINEMENT

private:
  /**
   *  Knot coordinates along each coordinate direction.
   *  For each direction the knot coordinates are sorted in an increasing
   *  order and does not contain duplicate values.
   */
  KnotCoordinates knot_coordinates_;

  /** Boundary ids, one id per face */
  SafeSTLArray<boundary_id, UnitElement<dim_>::template num_elem<dim_-1>()> boundary_id_;

  /**
   * Properties assigned to the elements.
   * Elements are referred to by their flat ids.
   */
  PropertyList elem_properties_;

  /**
   * Unique identifier associated to each object instance.
   */
  Index object_id_;

  TensorSizedContainer<dim_> elems_size_;

  template <int sdim>
  class SubGridData
  {
  private:
    using SubGrid = Grid<sdim>;

    //TODO (martinelli, November 11, 2015): the number of subgrid must be the number of combination of dim elems taken sdim at time
    SafeSTLArray<std::shared_ptr<SubGrid>,dim_> sub_grids_;

    using SubGridElemMap =
      SafeSTLMap<typename SubGrid::IndexType, IndexType>;

    SafeSTLArray<SubGridElemMap,UnitElement<dim>::template num_elem<sdim>()> sub_grid_element_maps_;
  };

  SubGridData<(dim_>0)? dim_-1 : 0> sub_grid_data_;

public:
  /**
   * Create an element (defined on this grid) with a given index and the given property
   */
  std::unique_ptr<ElementAccessor>
  create_element(const ListIt &index, const PropId &property) const;


  /**
   * Returns the unique identifier associated to each object instance.
   */
  Index get_object_id() const;

  /**
   * @name Functions related to the management/query of the element properties.
   */
  ///@{
  /**
  * Adds a new <tt>property</tt> definition for the elements in the
  * Grid.
  *
  * @note If the <tt>property</tt> is already present, an assertion will
  *  be raised (in Debug mode).
  */
  void add_property(const PropId &property);


  /**
   * Returns a const-reference to the list of element ids with the property specified by @p prop.
   */
  const List &get_elements_with_property(const PropId &prop) const;

  /**
   * Returns a reference to the list of element ids with the property specified by @p prop.
   */
  List &get_elements_with_property(const PropId &prop);


  /**
   * Sets the @p status of the given @p property for the element with id @p elem_id.
   */
  void set_property_status_elem(const PropId &property,
                                const IndexType &elem_id,
                                const bool status);

  /**
   * Returns true if the element identified by <tt>elem_id</tt> has
   * the property <tt>prop</tt>.
   */
  bool element_has_property(const IndexType &elem_id,
                            const PropId &prop) const;

#if 0
  /**
   * Sets the @p status of the given @p property for the entire set of elements in the grid.
   */
  void set_all_elements_property_status(const PropId &property,
                                        const bool status);
#endif

  ///@}


private:
  /**
   * Returns the flat ids of the sub-elements corresponding to the element
   * with index @p elem_id, referred to a Grid built as a refinement
   * of the current one using@p n_sub_elems for each element.
   */
  SafeSTLVector<Index>
  get_sub_elements_id(const TensorSize<dim_> &n_sub_elems,
                      const Index elem_id) const;




#ifdef MESH_REFINEMENT
  /**
   * This class member is the grid before the last refinement. If no
   * refinement is performed, this is a null pointer.
   */
  std::shared_ptr<self_t> grid_pre_refinement_ = nullptr;

  /**
   * Signals for the insert_knots() invocations. It can be viewed as a FIFO list of
   * function pointers.
   */
  signal_insert_knots_t insert_knots_signals_;
#endif
  friend class GridElement<dim_>;
//  friend class NonConstGridElement<dim_>;
//  friend class ConstGridElement<dim_>;

#ifdef SERIALIZATION
private:
  /**
   * @name Functions needed for serialization
   * @see <a href="http://uscilab.github.io/cereal/index.html">Cereal serialization library</a>
   */
  ///@{
  friend class cereal::access;

  template<class Archive>
  void serialize(Archive &ar);

  ///@}
#endif // SERIALIZATION


};

IGA_NAMESPACE_CLOSE

//#include <igatools/geometry/grid.serial>

#endif /* CARTESIAN_GRID_H_ */
