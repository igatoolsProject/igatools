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

#ifndef __NEW_BSPLINE_SPACE_H_
#define __NEW_BSPLINE_SPACE_H_

#include <igatools/base/config.h>
#include <igatools/base/logstream.h>

#include <igatools/basis_functions/spline_space.h>
#include <igatools/basis_functions/dof_distribution.h>
#include <igatools/basis_functions/bernstein_extraction.h>
#include<igatools/geometry/new_mapping.h>
#include<igatools/geometry/new_push_forward.h>
#include <igatools/basis_functions/new_physical_space.h>

IGA_NAMESPACE_OPEN

class SpaceManager;

template < int, int, int> class BSplineElement;
template < int, int, int> class BSplineElementHandler;
/**
 * Multivariate (tensor product) scalar, vector or k-tensor
 * valued B-spline space.
 * This object can be thought of as providing the
 * B-spline basis functions of a spline space.
 * The space is determined by:
 * - the knot vectors (implemented in the class CartesianGrid)
 * - the multiplicity vectors
 * - and the degree
 *
 * NewBSplineSpace allows the use of different multiplicity vectors
 * and degrees for each direction and for each component.
 *
 * \section const Constructing a NewBSplineSpace
 * Similar to the mechanism use in CartesianGrid we use
 * the create technique to create a smartpointer for each constructor
 * the class provides.
 * @todo enter a glossary for create idiom technique and refer from here
 * \code
 * auto space = NewBSplineSpace<dim>::create();
 * \endcode
 *
 * \section eval Evaluating basis function
 * Basis function are evaluated iterating on the elements.
 * Similarly we can obtain other information such as the local to global.
 *
 * \section dofs Degrees of freedom (dofs)
 * Each basis function is assigned a degree of freedom (an integer),
 * the order of this assignment is done following some prescribed
 * recipe or policy whose knowledge is not Really needed by most users.
 * They are internally stored in a grid-like multiarray container,
 * called the index space.
 * It works together with the index_space_offset which for each element
 * provides a view of the index
 * space to know which
 * are the non-zero basis function on each element, what we generally
 * refer to as the local to global mapping.
 *
 * \section bezier Storage of the basis functions (Bezier Extraction)
 * The basis functions on each element are stored in the Bspline space
 * as the 1D Bezier extraction operator.
 * When they need to be evaluated the operator applied to the
 * Berenstein polynomials,
 * allows to compute the values at the quadrature points.
 *
 * @todo write a module about cache optimization and handling.
 *
 * \section hom_range Optimizing homogeneous range type vector spaces
 *
 * \author martinelli, 2012, 2013, 2014
 * \author pauletti, 2012, 2013, 2014
 *
 *
 * @ingroup containers
 */
template<int dim_, int range_ = 1, int rank_ = 1>
class NewBSplineSpace :
    public std::enable_shared_from_this<NewBSplineSpace<dim_,range_,rank_> >,
    public SplineSpace<dim_, range_, rank_>
{
private:
    using BaseSpace = SplineSpace<dim_, range_, rank_>;

    /** Type for current class. */
    using self_t = NewBSplineSpace<dim_,range_,rank_>;

public:
    /** see documentation in \ref FunctionSpaceOnGrid */
    using PushForwardType = NewPushForward<Transformation::h_grad,dim_,0>;
    using PhysSpace = NewPhysicalSpace<self_t, 0, Transformation::h_grad>;

    /** Required type for space templated functions */
    using RefSpace = self_t;

    using GridType = CartesianGrid<dim_>;
    using ElementHandler = BSplineElementHandler<dim_, range_, rank_>;

    static const int dim       = PushForwardType::dim;
    static const int codim     = PushForwardType::codim;
    static const int space_dim = PushForwardType::space_dim;
    static const int range     = range_;
    static const int rank      = rank_;

    // TODO (pauletti, Oct 16, 2014): delete this
    static const iga::RefSpaceType ref_space_type = iga::RefSpaceType(0);

    using BaseSpace::n_components;
    using BaseSpace::components;
    using BaseSpace::dims;

    static const bool has_weights = false;

public:
    using typename BaseSpace::Func;
    using typename BaseSpace::Point;
    using typename BaseSpace::Value;
    using typename BaseSpace::Derivative;
    using typename BaseSpace::Div;

    using RefPoint = Point;

public:
#if 0
    /** Type for the reference face space.*/
    using RefFaceSpace = Conditional<(dim>0),
          NewBSplineSpace<dim-1,range,rank>,
          NewBSplineSpace<0,range,rank> >;

    using FaceSpace = NewPhysicalSpace<RefFaceSpace,
          typename PushForwardType::FacePushForward>;
#endif

    /** Type for the element accessor. */
    using ElementAccessor = BSplineElement<dim,range,rank>;

    /** Type for iterator over the elements.  */
    using ElementIterator = GridForwardIterator<ElementAccessor>;


    using typename BaseSpace::InteriorReg;

    using typename BaseSpace::DegreeTable;
    using typename BaseSpace::MultiplicityTable;
    using typename BaseSpace::KnotsTable;
    using typename BaseSpace::SpaceDimensionTable;
    using typename BaseSpace::EndBehaviour;
    using typename BaseSpace::EndBehaviourTable;


public:
    /**
     * @name Creators.
     */
    ///@{
    /**
     * Builds and returns a maximum regularity BSpline space
     * over CartesianGrid
     * @p knots for the given @p degree in all directions and homogeneous
     * in all components.
     */
    static std::shared_ptr<self_t>
    create(const int degree, std::shared_ptr<GridType> knots);

    /**
     * Builds and returns a maximum regularity BSpline space over CartesianGrid
     * @p knots for the given @p degree[i] in the i-th direction and homogeneous
     * in all components.
     */
    static std::shared_ptr<self_t>
    create(const TensorIndex<dim> &degree, std::shared_ptr<GridType> knots);

    /**
     * Builds and returns a maximum regularity BSpline space over CartesianGrid
     * @p knots for the given @p degree for each direction and for each
     * component.
     */
    static std::shared_ptr<self_t>
    create(const DegreeTable &degree,
           std::shared_ptr<GridType> knots,
           const bool homogeneous_range = false);

    /**
     * Builds and returns a BSpline space over the CartesianGrid
     * @p knots with the given multiplicity vector @p mult_vectors
     * for each component
     * and the given @p degree for each direction and for each
     * component.
     */
    static std::shared_ptr<self_t>
    create(const DegreeTable &deg,
           std::shared_ptr<GridType> knots,
           std::shared_ptr<const MultiplicityTable> interior_mult,
           const EndBehaviourTable &ends = EndBehaviourTable());
    ///@}

    /** Destructor. */
    ~NewBSplineSpace() = default;

protected:
    /** @name Constructors */
    ///@{
    /**
     * Constructs a maximum regularity BSpline space over CartesianGrid
     * @p knots for the given @p degree in all directions and homogeneous
     * in all components.
     */
    explicit NewBSplineSpace(const int degree, std::shared_ptr<GridType> knots);

    /**
     * Constructs a maximum regularity BSpline space over CartesianGrid
     * @p knots for the given @p degree[i] in the i-th direction and homogeneous
     * in all components.
     */
    explicit NewBSplineSpace(const TensorIndex<dim> &degree,
                             std::shared_ptr<GridType> knots);

    /**
     * Constructs a maximum regularity BSpline space over CartesianGrid
     * @p knots for the given @p degree for each direction and for each
     * component.
     */
    explicit NewBSplineSpace(const DegreeTable &degree,
                             std::shared_ptr<GridType> knots,
                             const bool homogeneous_range = false);

    /**
     * Constructs a BSpline space over the CartesianGrid
     * @p knots with the given multiplicity vector @p mult_vectors
     * for each component
     * and the given @p degree for each direction and for each
     * component.
     */
    explicit NewBSplineSpace(const DegreeTable &deg,
                             std::shared_ptr<GridType> knots,
                             std::shared_ptr<const MultiplicityTable> interior_mult,
                             const EndBehaviourTable &ends);

    /**
     * Copy constructor. Not allowed to be used.
     */
    NewBSplineSpace(const self_t &space) = delete;
    ///@}

    /** @name Assignment operators */
    ///@{
    /** Copy assignment. Not allowed to be used. */
    self_t &
    operator=(const self_t &space) = delete;
    ///@}

public:
    // TODO (pauletti, Oct 16, 2014): need to be documented or deleted, check!
    vector<Index> get_loc_to_global(const CartesianGridElement<dim> &element) const;

    vector<Index> get_loc_to_patch(const CartesianGridElement<dim> &element) const;

    /** @name Functions involving the element iterator */
    ///@{
    /**
     * Returns a element iterator to the first element of the patch
     */
    ElementIterator begin() const;

    /**
     * Returns a element iterator to the last element of the patch
     */
    ElementIterator last() const;


    /**
     * Returns a element iterator to one-pass the end of patch.
     */
    ElementIterator end() const;
    ///@}


#if 0
    /** Getting some underlying objects */
    ///@{
    std::shared_ptr<RefFaceSpace>
    get_ref_face_space(const Index face_id,
                       vector<Index> &face_to_element_dofs,
                       typename GridType::FaceGridMap &elem_map) const;


    std::shared_ptr<FaceSpace>
    get_face_space(const Index face_id,
                   vector<Index> &face_to_element_dofs) const;
#endif
    std::shared_ptr<const self_t >
    get_reference_space() const;

#if 0
    /** Return the push forward (non-const version). */
    std::shared_ptr<PushForwardType> get_push_forward();


    /** Return the push forward (const version). */
    std::shared_ptr<const PushForwardType> get_push_forward() const;

#endif
    std::shared_ptr<SpaceManager> get_space_manager();


    std::shared_ptr<const SpaceManager> get_space_manager() const;


    /** Returns the container with the global dof distribution (const version). */
    const DofDistribution<dim, range, rank> &
    get_dof_distribution_global() const;

    /** Returns the container with the global dof distribution (non const version). */
    DofDistribution<dim, range, rank> &
    get_dof_distribution_global();

    /** Returns the container with the patch dof distribution (const version). */
    const DofDistribution<dim, range, rank> &
    get_dof_distribution_patch() const;


    /** Returns the container with the patch dof distribution (non const version). */
    DofDistribution<dim, range, rank> &
    get_dof_distribution_patch();
    ///@}

    /**
     * Adds an @p offset to the values of the dof ids.
     */
    void add_dofs_offset(const Index offset);

    /**
     * This function returns the global dof id corresponding to the basis function
     * with tensor index <p>tensor_index</p> on the @p comp component of the space.
     */
    Index
    get_global_dof_id(const TensorIndex<dim> &tensor_index,
                      const Index comp) const;

    /**
     * Prints internal information about the space.
     * @note Mostly used for debugging and testing.
     */
    void print_info(LogStream &out) const;

private:

    /** Container with the local to global basis indices
     * @note The concept of global indices refers to a global numeration of the
     * dofs of all the spaces.
     */
    DofDistribution<dim, range, rank> dof_distribution_global_;

    /** Container with the local to patch basis indices
     * @note The concept of patch indices refers to the numeration at patch
     * level of the dofs.
     */
    DofDistribution<dim, range, rank> dof_distribution_patch_;

    /** @name Bezier extraction operator. */
    BernsteinExtraction<dim, range, rank> operators_;


    friend class BSplineElement<dim, range, rank>;
    friend class BSplineElementHandler<dim, range, rank>;

    /**
     * Refines the function space after a grid uniform refinement.
     *
     * @param[in] refinement_directions Directions along which the refinement is performed.
     * @param[in] grid_old Grid before the refinement.
     *
     * @pre Before invoking this function, must be invoked the function grid_->refine().
     * @note This function is connected to the CartesianGrid's signal for the refinement, and
     * it is necessary in order to avoid infinite loops in the refine() function calls.
     *
     * @ingroup h_refinement
     */
    void refine_h_after_grid_refinement(
        const std::array<bool,dim> &refinement_directions,
        const GridType &grid_old) ;

public:
    DeclException1(ExcScalarRange, int,
                   << "Range " << arg1 << "should be 0 for a scalar valued"
                   << " space.");
};

IGA_NAMESPACE_CLOSE

#endif /* __BSPLINE_SPACE_H_ */
