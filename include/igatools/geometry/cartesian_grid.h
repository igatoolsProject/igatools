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

#ifndef CARTESIAN_GRID_H_
#define CARTESIAN_GRID_H_

#include <igatools/base/config.h>
#include <igatools/base/tensor.h>
#include <igatools/base/logstream.h>
#include <igatools/utils/cartesian_product_array.h>
#include <igatools/utils/dynamic_multi_array.h>
#include <igatools/geometry/unit_element.h>
#include <igatools/base/array_utils.h>
#include <igatools/geometry/cartesian_grid_iterator.h>

// TODO (pauletti, Oct 9, 2014): should we use iga array
#include <array>
#include <memory>
#include <map>
#include <set>

#include <boost/signals2.hpp>

IGA_NAMESPACE_OPEN

template <int> class CartesianGridElement;
template <int> class GridElementHandler;
/**
 * @brief Grid in <tt>dim</tt>-dimensional space with cartesian-product structure.
 *
 * The data defining the grid is a <tt>dim</tt>-dimensional array of coordinates,
 * representing the positions of the grid wires.
 *
 * Two consecutive coordinates (along the same direction) defines an <em>interval</em>.
 *
 * Then, the tensor-product of the intervals along the coordinate directions define
 * the elements that are tiling the domain covered by the CartesianGrid;
 *
 * The elements can have associated a certain list of ElementProperty,
 * and then the list of elements with a given property can be extracted from the CartesianGrid.
 *
 * The element type for the CartesianGrid is CartesianGridElement.
 *
 * The elements can be iterated with the CartesianGridIterator object.
 *
 * ### Getting a CartesianGrid by an XML structure.
 * A CartesianGrid object can be obtained by a Boost XML structure using the
 * function
 * get_cartesian_grid_from_xml().
 * An example of valid XML structure for a 2-dimensional CartesianGrid is given
 * by the following
 * XML block
 * @code{.xml}
   <CartesianGrid Dim="2">
     <Knots Direction="0" Size="2">
       0.0 1.0
     </Knots>
     <Knots Direction="1" Size="2">
        0.0 1.0
     </Knots>
   </CartesianGrid>
 * @endcode
 * where the <tt>Dim</tt> attribute in the <tt>"CartesianGrid"</tt> tag is the
 * dimensionality of the CartesianGrid
 * and the <tt>Size</tt> attribute in each <tt>"Knots"</tt> tag is the number
 *  of unique knots values (also called
 * "breakpoints") along each coordinate direction.
 *
 * @author M. Martinelli 2012, 2013, 2014
 * @author pauletti 2012, 2013, 2014
 *
 * @see get_cartesian_grid_from_xml()
 * @ingroup h_refinement
 * @ingroup containers
 * @todo document more
 */
template<int dim_>
class CartesianGrid :
    protected TensorSizedContainer<dim_>,
    public std::enable_shared_from_this<CartesianGrid<dim_>>
{
private:
    /** Type for current class. */
    using self_t = CartesianGrid<dim_>;

public:

    /**
     * Alias for the (static) class holding the topological information.
     */
    using Topology = UnitElement<dim_>;

    using Point = Points<dim_>;

    static const int dim = dim_;

    /** Type for the element accessor. */
    using ElementAccessor = CartesianGridElement<dim_>;

    /** Type for the iterator over the elements of the grid (non-const version).  */
    using ElementIterator = CartesianGridIterator<ElementAccessor>;

    /** Type for the iterator over the elements of the grid (const version).  */
    using ElementConstIterator = CartesianGridConstIterator<ElementAccessor>;

    using ElementHandler = GridElementHandler<dim_>;

    /** Type for the vector of knot vectors */
    using KnotCoordinates = CartesianProductArray<Real, dim_>;

    /**
     * Types of grid for future optimization
     */
    enum class Kind
    {
        uniform, direction_uniform, non_uniform
    };



    /** @name Constructors*/
    ///@{
protected:
    /**
     * Construct a uniform cartesian grid of the unit <tt>dim</tt>-dimensional
     * hypercube \f$[0,1]^{dim}\f$, with @p n knots (equally spaced) in each dimension.
     */
    explicit CartesianGrid(const Size n = 2);

    /**
     * Construct a uniform cartesian grid of the unit <tt>dim</tt>-dimensional
     * hypercube \f$[0,1]^{dim}\f$,
     * with <tt>n[0],..,n[dim-1</tt>] knots in each dimension
     * respectively.
     */
    explicit CartesianGrid(const TensorSize<dim_> &n_knots, const Kind kind);

    /**
     * @todo Document me
     */
    explicit CartesianGrid(const BBox<dim_> &end_points,
                           const Size n_knots,
                           const Kind kind);

    /**
     * @todo Document me
     */
    explicit CartesianGrid(const BBox<dim_> &end_points,
                           const TensorSize<dim_> &n_knots,
                           const Kind kind);

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
    CartesianGrid(const KnotCoordinates &knot_coordinates,
                  const Kind kind);

    /**
     * Construct a cartesian grid where the knot coordinate in each
     * direction is provided as std::array of vector<Real>.
     *
     * The knot coordinate in each direction must be sorted and without
     * repetition.
     * @note In Debug mode, a check for this precondition (up to machine precision)
     * is perform and if not satistified an exception is raised.
     */
    explicit
    CartesianGrid(const std::array<vector<Real>,dim_> &knot_coordinates);

public:
    /**
     * Copy constructor.
     * Perform a deep copy of the member variables except the
     * signal_refine_ variable, that is not copied at all.
     */
    CartesianGrid(const self_t &grid);

    /**  Move constructor */
    CartesianGrid(self_t &&grid) = default;

    /** Destructor */
    ~CartesianGrid() = default;
    ///@}

    /**
     * @name Creators
     * @note The functions here return a CartesianGrid<dim_> object wrapped by
     * a std::shared_ptr
     */
    ///@{

    /**
     * Creates a uniform cartesian grid of the unit <tt>dim</tt>-dimensional
     * hypercube \f$[0,1]^{dim}\f$, with @p n knots (equally spaced) in each dimension.
     */
    static std::shared_ptr<self_t> create(const Index n = 2);


    /**
     * Creates a uniform cartesian grid of the unit <tt>dim</tt>-dimensional
     * hypercube \f$[0,1]^{dim}\f$,
     * with <tt>n[0],..,n[dim-1</tt>] knots in each dimension
     * respectively.
     */
    static std::shared_ptr<self_t> create(const TensorSize<dim_> &n);


    /**
     * Construct a cartesian grid where the knot coordinate in each
     * direction is provided as CartesianProductArray object.
     *
     * The knot coordinate in each direction must be sorted and without
     * repetition.
     * @note In Debug mode, a check for this precondition
     * (up to machine precision)
     * is perform and if not satisfied an exception is raised.
     */
    static std::shared_ptr<self_t>
    create(const KnotCoordinates &knot_coordinates);

    /**
     * Construct a cartesian grid where the knot coordinate in each
     * direction is provided as std::array of vector<Real>.
     *
     * The knot coordinate in each direction must be sorted and without
     * repetition.
     * @note In Debug mode, a check for this precondition
     * (up to machine precision)
     * is perform and if not satisfied an exception is raised.
     */
    static std::shared_ptr<self_t>
    create(const std::array<vector<Real>,dim_> &knot_coordinates);


    static std::shared_ptr<self_t>
    create(const BBox<dim_> &end_points, const Size n_knots);

    /**
     * @todo document me
     */
    static std::shared_ptr<self_t>
    create(const BBox<dim_> &end_points, const TensorSize<dim_> &n);

    /**
     * Creates a CastesianGrid object (wrapped by a shared_ptr) using
     * the copy constructor.
     *
     * @note A grid built in this way is totally uncoupled from the grid used as argument
     * of this function. For example, a refinement of a grid does not affect the other gird.
     */
    static std::shared_ptr<self_t>
    create(const self_t &grid);

    ///@}


    /**
     * Create an element (defined on this grid) with a given flat_index.
     */
    std::shared_ptr<ElementAccessor> create_element(const Index flat_index) const;

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


    ///@name Getting grid information
    ///@{
    /**
     * Returns the number of elements with the <tt>property</tt> specified as input argument.
     */
    Size get_num_elements_same_property(const ElementProperty &property) const;

    /**
     * Total number of active elements.
     */
    Size get_num_active_elems() const;

    /**
     * Total number of influence elements.
     */
    Size get_num_influence_elems() const;

    /** Total number of elements, including active and non-active */
    Size get_num_all_elems() const;

    /**
     * Total number of one dimensional intervals along each
     * coordinate direction.
     */
    TensorSize<dim_> get_num_intervals() const;

    /**
     * Query the number of knot values along each coordinate direction
     * represented by the CartesianGrid.
     */
    TensorSize<dim_> get_num_knots_dim() const;

    /**
     * Returns the knot coordinates along the direction @p i.
     */
    vector<Real> const &get_knot_coordinates(const int i) const;

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
     * domain represented by the CartesianGrid object.
     */
    BBox<dim_> get_bounding_box() const;
    ///@}


    ///@name Iterating of grid elements
    ///@{
    /**
     * This function returns a element iterator to the first element of the patch.
     */
    ElementIterator begin();

    /**
     * This function returns a element iterator to one-pass the end of patch.
     */
    ElementIterator end();

    /**
     * This function returns a element (const) iterator to the first element of the patch.
     */
    ElementConstIterator begin() const;

    /**
     * This function returns a element (const) iterator to one-pass the end of patch.
     */
    ElementConstIterator end() const;

    /**
     * This function returns a element (const) iterator to the first element of the patch.
     */
    ElementConstIterator cbegin() const;

    /**
     * This function returns a element (const) iterator to one-pass the end of patch.
     */
    ElementConstIterator cend() const;


    /**
     * This function returns the iterator to the last active element on the grid.
     */
    ElementIterator last();

    /**
     * This function returns the (const) iterator to the last active element on the grid.
     */
    ElementConstIterator last() const;

    /**
     * This function returns the (const) iterator to the last active element on the grid.
     */
    ElementConstIterator clast() const;
    ///@}


    /** @name Functions for the index transformations */
    ///@{
    /**
     * Transformation from a tensor-index to a flat-index.
     */
    Index tensor_to_flat(const TensorIndex<dim_> &tensor_index) const;

    /**
     * Transformation from a flat-index to a tensor-index.
     */
    TensorIndex<dim_> flat_to_tensor(const Index flat_index) const;
    ///@}


    ///@name Dealing with boundary information
    ///@{
    /**
     * Get the patch @p face boundary id.
     */
    boundary_id get_boundary_id(const int face) const;

    /**
     * Set the patch @p face to have the boundary @p id.
     */
    void set_boundary_id(const int face, const boundary_id id);


    template<int sub_dim>
    using BoundaryNormal = std::array<Points<dim_>, dim_-sub_dim>;
    /**
     * Returns the outward pointing
     * unit normal vector space to the element of sub dim_ k.
     */
    template<int sub_dim>
    BoundaryNormal<sub_dim> get_boundary_normals(const int s_id) const;

    template<int k>
    using InterGridMap = std::map<typename CartesianGrid<k>::ElementIterator, ElementConstIterator>;

    /**
     * Construct a sub grid of dimension k conforming to
     * the grid sub element sub_elem_id and a map from the elements of
     * the sub_element grid to the corresponding element of the current
     * grid.
     */
    template<int k>
    std::shared_ptr<CartesianGrid<k> >
    get_sub_grid(const int sub_elem_id, InterGridMap<k> &elem_map) const;

    ///@}

    /**
     * Given a vector of points, this function return a map with
     * entries indexed by the grid element each point belongs to
     * containing a list of indices of the points that belong to
     * this element.
     * For ex.
     * @code
     *
     * @endcode
     */
    std::map<ElementIterator, vector<int> >
    find_elements_of_points(const ValueVector<Points<dim_>> &points) const;

public:
    /**
     * Prints debug information of the CartesianGrid to a LogStream.
     */
    void print_info(LogStream &out) const;

    /**
     * Comparison operator. Returns true if the knot coordinates of two grid
     * are equal.
     */
    bool operator==(const CartesianGrid<dim_> &grid) const;

private:
    /** Type for the refinement signal. */
    using signal_refine_t =
        boost::signals2::signal<
        void (const std::array<bool,dim_> &,const CartesianGrid<dim_> &)>;

public:

    /** Slot type for the refinement signal. */
    using SignalRefineSlot = typename signal_refine_t::slot_type;

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
        const std::array<bool,dim_> &refinement_directions,
        const std::array<Size,dim_> &n_subdivisions);

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
     *  emitted whenever a refine() function is called by an object holding
     *  a CartesianGrid member.
     */
    boost::signals2::connection
    connect_refinement(const SignalRefineSlot &subscriber);

    /**
     * Returns the grid before the last refinement. If no refinement is
     * performed, this function returns a null pointer.
     */
    std::shared_ptr<const self_t > get_grid_pre_refinement() const;
    ///@}

private:
    /** Flag for optimization use */
    Kind kind_ = Kind::non_uniform;

    /** Boundary ids, one id per face */
    std::array<boundary_id, UnitElement<dim_>::template num_elem<dim_-1>()> boundary_id_;

    /**
     *  Knot coordinates along each coordinate direction.
     *  For each direction the knot coordinates are sorted in an increasing
     *  order and does not contain duplicate values.
     */
    KnotCoordinates knot_coordinates_;


private:
    /**
     * Container for the element ids having a certain property.
     *
     * The property is the key of the std::map.
     */
    std::map<ElementProperty,std::set<Index>> properties_elements_id_;

    /**
     * Returns true if the element identified with <p>elem_flat_id</p> has
     * the ElementProperty <p>property</p>.
     */
    bool test_if_element_has_property(const Index elem_flat_id, const ElementProperty &property) const;

public:

    /**
     * @name Functions related to the management/query of the element properties.
     */
    ///@{
    /**
     * Returns true if the element identified with <p>elem_flat_id</p> is active.
     */
    bool is_element_active(const Index elem_flat_id) const;

    /**
     * Returns the flat id of the elements having a certain @p property (non-const version).
     */
    std::set<Index> &get_elements_id_same_property(const ElementProperty &property);

    /**
     * Returns the flat id of the elements having a certain @p property (const version).
     */
    const std::set<Index> &get_elements_id_same_property(const ElementProperty &property) const;

    /**
     * Sets the @p status of the given @p property for the element with flat id @ elem_flat_id.
     */
    void set_element_property(const ElementProperty &property,
                              const Index elem_flat_id,
                              const bool status);
    ///@}


private:
    /**
     * Returns the flat ids of the sub-elements corresponding to the element with index @p elem_id,
     * referred to a CartesianGrid built as a refinement of the current one using
     * @p n_sub_elems for each element.
     */
    vector<Index> get_sub_elements_id(const TensorSize<dim_> &n_sub_elems, const Index elem_id) const;

    /**
     * Perform a uniform refinement of the knots along the @p direction_id
     * direction,
     * dividing each interval in the knot vector into @p n_subdivisions
     * intervals.
     * @param[in] direction_id Direction along which the refinement is
     * performed.
     * @param[in] n_subdivisions Number of subdivision in which each interval
     * in the knot vector is
     * divided. This value must be >= 2.
     */
    void refine_knots_direction(const int direction_id,
                                const Size n_subdivisions);

    /**
     * This class member is the grid before the last refinement. If no
     * refinement is performed, this is a null pointer.
     */
    std::shared_ptr<const self_t> grid_pre_refinement_ = nullptr;

    /**
     * Signals for the h-refinement. It can be viewed as a FIFO list of
     * function pointers.
     */
    signal_refine_t refine_signals_;

    friend class CartesianGridElement<dim_>;
};

//template<int dim_>
//std::shared_ptr<CartesianGrid::SubGrid<dim_-1>>
//get_face_grid(const int sub_elem_id, InterGridMap<k> &elem_map) const;

IGA_NAMESPACE_CLOSE

#endif /* CARTESIAN_GRID_H_ */
