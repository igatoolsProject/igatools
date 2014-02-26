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

#ifndef CARTESIAN_GRID_H_
#define CARTESIAN_GRID_H_

#include <igatools/base/config.h>
#include <igatools/base/tensor.h>
#include <igatools/utils/cartesian_product_array.h>
#include <igatools/base/logstream.h>
#include <igatools/geometry/grid_forward_iterator.h>
#include <igatools/geometry/unit_element.h>

#include <array>
#include <vector>
#include <memory>

#include <boost/signals2.hpp>

IGA_NAMESPACE_OPEN

template <int> class CartesianGridElementAccessor;


/**
 * This class represents a the knotspan or knot vector.
 * View with a tensor product structure it is also known as
 * a rectilinear grid in the dim-dimensional space.
 *
 *
 * \section sec1 Creating a patch
 * - Example 1: creating a 3-dimensional grid in the \f$ [0,1] \times [0,1] \times [0,1] \f$ domain
 * with 4 uniform intervals in each coordinate direction
 * \code
 *    CartesianGrid<3> grid( 5 ); // 4 intervals == 5 knot values
 * \endcode
 * - Example 2: creating a 2-dimensional grid in the \f$ [0,1] \times [0,1] \f$ domain
 * with 4 uniform intervals in the coordinate direction 0 and 2 uniform intervals in
 * the coordinate direction 1
 * \code
 *   array< int, 2 > num_knots;
 *   num_knots[0] = 5;
 *   num_knots[1] = 3;
 *   CartesianGrid<2> grid( num_knots );
 * \endcode
 *
 *
 * \section sec2 Accessing data through Iterators and Accessors
 * The data of the patch is access through iterators (STL like).
 * The main advantage of using iterators are:
 * - you can write code independent of the dimension
 * - you can work on the patch data without worrying of its internal data structure
 *
 * \subsection ex Examples
 * - Count the number of elements on a patch
 * \code
 *   CartesianGrid<dim>::ElementIterator element = patch.begin();
 *   CartesianGrid<dim>::ElementIterator endelement = patch.end();
 *   int n_elements=0;
 *   for (; element!=endelement; ++element)
 *       ++n_elements;
 * \endcode
 *
 *   \author Massimiliano Martinelli (massimiliano.martinelli@gmail.com), 2012, 2013
 *   \author pauletti 2012, 2013
 *
 * @ingroup refinement
 *
 */
template< int dim_ >
class CartesianGrid
{
public:
    /** Dimensionality of the grid. */
    static const int dim = dim_ ;

    /**
     * Types of grid for future optimization
     */
    enum class Kind
    {
        uniform, direction_uniform, non_uniform
    };

    /** Type for the face of the grid */
    using FaceType = Conditional<(dim>0),CartesianGrid<dim-1>,CartesianGrid<0> >;


    /** Type for iterator over the elements.*/
    using ElementIterator = GridForwardIterator<CartesianGridElementAccessor<dim> >;

    /** @name Constructors*/
    ///@{
    /**
     * Construct a uniform cartesian grid of the unit <b>d</b>-dimensional
     * cube [0,1]^d, with \b n knots (equally spaced) in each dimension.
     */
    explicit CartesianGrid(const Size n = 2);

    /**
     * Construct a uniform cartesian grid of the unit <b>d</b>-dimensional
     * cube [0,1]^d, with <b>n</b>[0],..,<b>n</b>[dim-1] knots in each dimension
     * respectively.
     */
    explicit CartesianGrid(const TensorSize<dim> &n_knots);

    /**
     * Construct a cartesian grid where the knot coordinate in each
     * direction is provided.
     *
     * The knot coordinate in each direction must be sorted and without
     * repetition.
     * @note In debug mode, a check for this precondition (up to machine precision)
     * is perform and if not satistified an exception is raised.
     */
    explicit
    CartesianGrid(const CartesianProductArray<Real,dim> &knot_coordinates);


    /**
     * @todo Document me
     */
    explicit CartesianGrid(const BBox<dim> &end_points,
                           const TensorSize<dim> &n_knots);

    /**
     * Copy constructor.
     * Perform a deep copy of the member variables except the
     * signal_refine_ variable, that is not copied at all.
     */
    CartesianGrid(const CartesianGrid<dim_> &grid);

    /**  Move constructor */
    CartesianGrid(CartesianGrid<dim_> &&grid) = default;

    /** Destructor */
    ~CartesianGrid() = default;
    ///@}


    /** @name Creators */
    ///@{
    /**
     * Create a uniform cartesian grid of the unit <b>d</b>-dimensional
     * cube [0,1]^d, with <b>n</b>[0],..,<b>n</b>[dim-1] knots in each dimension
     * respectively.
     */
    static std::shared_ptr< CartesianGrid<dim_> > create(const TensorSize<dim> &n);


    /**
     * Construct a cartesian grid where the knot coordinate in each
     * direction is provided.
     *
     * The knot coordinate in each direction must be sorted and without
     * repetition.
     * \note n debug mode, a check for this precondition (up to machine precision)
     * is perform and if not satistified an exception is raised.
     */
    static std::shared_ptr< CartesianGrid<dim_> >
    create(const CartesianProductArray<Real,dim> &knot_coordinates) ;



    /**
     * Create a uniform cartesian grid of the unit <b>d</b>-dimensional
     * cube [0,1]^d, with \b n knots in each dimension.
     */
    static std::shared_ptr<CartesianGrid<dim_> > create(const Index n = 2);



    static std::shared_ptr< CartesianGrid<dim_> >
    create(const BBox<dim> &end_points,
           const TensorSize<dim> &n);
    ///@}


    /**
     * @name Assignment operators
     */
    ///@{

    /**
     * Copy assignment operator.
     */
    CartesianGrid<dim_> &operator=(const CartesianGrid<dim_> &grid) = default;

    /**
     * Move assignment operator.
     */
    CartesianGrid<dim_> &operator=(CartesianGrid<dim_> &&grid) = default;

    ///@}


    ///@name Getting grid information
    ///@{
    /**
     * Query the total number of elements.
     */
    Size get_num_elements() const;

    /**
     * Query the total number of elements along each coordinate direction.
     */
    TensorSize<dim> get_num_elements_dim() const;

    /**
     * Query the number of knot values along each coordinate direction
     * represented by the CartesianGrid.
     */
    TensorSize<dim> get_num_knots_dim() const;

    /**
     * Returns the knot coordinates along the direction \b i.
     */
    std::vector<Real> const &get_knot_coordinates(const int i) const;

    /**
     * Returns the knot coordinates along all the directions (const version).
     */
    CartesianProductArray<Real, dim> const &get_knot_coordinates() const;



    /**
     * Returns the knot coordinates along the direction \b i.
     */
    CartesianProductArray<Real, dim> get_element_lengths() const ;
    ///@}


    ///@name Iterating of grid elements
    ///@{
    /**
     * This function returns a element iterator to the first element of the patch.
     */
    ElementIterator begin() const;

    /**
     * This function returns a element iterator to one-pass the end of patch.
     */
    ElementIterator end() const;
    ///@}

    ///@name Dealing with boundary information
    ///@{
    /**
     * Get the patch \b face boundary id.
     */
    boundary_id get_boundary_id(const int face) const;

    /**
     * Set the patch \b face to have the boundary \b id.
     */
    void set_boundary_id(const int face, const boundary_id id);

    /**
     * Returns the outward pointing
     * unit normal vector to the face number @p face_no.
     */
    Point<dim> get_face_normal(const int face_no) const;

    FaceType get_face_grid(const int face_id, std::map<int,int> &elem_map) const;

    ///@}


    /**
     * Prints the CartesianGrid on a LogStream.
     */
    void print_info(LogStream &out) const;



private:

    Kind kind_ = Kind::non_uniform;

    /**
     * Boundary ids, one id per face
     */
    std::array< boundary_id, UnitElement<dim>::faces_per_element > boundary_id_;

    /**
     *  Knot coordinates along each coordinate direction.
     *  For each direction the knot coordinates are sorted in an increasing
     *  order and does not contain duplicate values.
     */
    CartesianProductArray<Real,dim> knot_coordinates_;

    /** Type for the refinement signal. */
    using signal_refine_t = boost::signals2::signal<
                            void (const std::array<bool,dim> &,const CartesianGrid<dim> &)>;

public:

    /** Slot type for the refinement signal. */
    using SignalRefineSlot = typename signal_refine_t::slot_type;

    /** @name Functions for performing grid refinement */
    ///@{
    /**
     * Perform a uniform refinement of the grid along the @p direction_id direction,
     * dividing each interval into @p n_subdivisions intervals.
     * @param[in] direction_id Direction along which the refinement is performed.
     * @param[in] n_subdivisions Number of subdivision in which each interval in the grid
     * along the specified direction is divided. This value must be >= 2.
     */
    void refine_direction(const int direction_id, const Size n_subdivisions);

    /**
     * Refine the cartesian grid and the objects connected to it (if any),
     * e.g. maps, spaces, etc.
     * along the directions specified by the true values in the entries of the
     * array of bools @p refinement_directions,
     * and with a number of subdivisions for each interval (along each direction)
     * specified by @p n_subdivisions.
     *
     * @note If the i-th direction is not active for the refinement
     * (i.e. <tt>refinement_directions[i] == false</tt>),
     *  then the corresponding value <tt>n_subdivisions[i]</tt> will be ignored.
     */
    void refine_directions(
        const std::array<bool,dim> &refinement_directions,
        const std::array<Size,dim> &n_subdivisions) ;

    /**
     * Refine the cartesian grid and the objects connected to it (if any),
     * e.g. maps, spaces, etc.
     *
     * Each interval in the unrefined grid is uniformly divided in @p n_subdivisions
     * sub-intervals.
     */
    void refine(const Size n_subdivisions = 2);

    /**
     *  Connect a slot (i.e. a function pointer) to the refinement signals which will be
     *  emitted whenever a refine() function is called by an object holding a CartesianGrid member.
     */
    boost::signals2::connection connect_refinement(const SignalRefineSlot &subscriber) ;
    ///@}


private:
    /**
     * Perform a uniform refinement of the knots along the @p direction_id direction,
     * dividing each interval in the knot vector into @p n_subdivisions intervals.
     * @param[in] direction_id Direction along which the refinement is performed.
     * @param[in] n_subdivisions Number of subdivision in which each interval in the knot vector is
     * divided. This value must be >= 2.
     */
    void refine_knots_direction(const int direction_id,
                                const Size n_subdivisions) ;

    /**
     * Signals for the h-refinement. It can be viewed as a FIFO list of function pointers.
     */
    signal_refine_t refine_signals_;

    friend class CartesianGridElementAccessor< dim >;
};


IGA_NAMESPACE_CLOSE

#endif /* CARTESIAN_GRID_H_ */
