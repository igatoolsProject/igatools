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

#ifndef __BSPLINE_SPACE_H_
#define __BSPLINE_SPACE_H_

#include <igatools/base/config.h>
#include <igatools/base/logstream.h>

#include <igatools/basis_functions/spline_space.h>
#include <igatools/basis_functions/dof_distribution.h>
#include <igatools/basis_functions/bernstein_extraction.h>

#include<igatools/geometry/mapping.h>
#include<igatools/geometry/push_forward.h>
#include <igatools/basis_functions/physical_space.h>

IGA_NAMESPACE_OPEN

class SpaceManager;

template < int, int, int> class BSplineElementAccessor;

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
 * BSplineSpace allows the use of different multiplicity vectors
 * and degrees for each direction and for each component.
 *
 * \section const Constructing a BSplineSpace
 * Similar to the mechanism use in CartesianGrid we use
 * the create technique to create a smartpointer for each constructor
 * the class provides.
 * @todo enter a glossary for create idiom technique and refer from here
 * \code
 * auto space = BSplineSpace<dim>::create();
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
 * When they need to be evaluated the operator applied to the Berenstein polynomials,
 * allows to compute the values at the quadrature points.
 *
 * @todo write a module about cache optimization and handling.
 *
 * \section hom_range Optimizing homogeneous range type vector spaces
 *
 * \author martinelli, 2012, 2013, 2014
 * \author pauletti, 2012, 2013
 *
 * @tparam dim Dimensionality of the parametric space (must be equal to the dimensionality
 * of the grid used top build the space
 *
 * @ingroup containers
 */
template<int dim_, int range_ = 1, int rank_ = 1>
class BSplineSpace :
    public std::enable_shared_from_this<BSplineSpace<dim_,range_,rank_> >,
    public SplineSpace<dim_, range_, rank_>
{
private:
    using BaseSpace = SplineSpace<dim_, range_, rank_>;

    /** Type for current class. */
    using self_t = BSplineSpace<dim_,range_,rank_>;

public:
    /** see documentation in \ref FunctionSpaceOnGrid */
    using PushForwardType = PushForward<Transformation::h_grad,dim_,0>;
    using PhysSpace = PhysicalSpace<self_t, PushForwardType>;

    /** Required type for space templated functions */
    using RefSpace = self_t;

    using GridType = typename PushForwardType::GridType;

    static const int dim = PushForwardType::dim;

    static const int codim = PushForwardType::codim;

    static const int space_dim = PushForwardType::space_dim;

    static const int range = range_;

    static const int rank = rank_;

    using BaseSpace::n_components;
    using BaseSpace::components;
    using BaseSpace::dims;

    static const bool has_weights = false;

public:
    using typename BaseSpace::Func;
    using typename BaseSpace::Derivative;
    using typename BaseSpace::Point;
    using typename BaseSpace::Value;
    using typename BaseSpace::Div;

public:
    /** Type for the reference face space.*/
    using RefFaceSpace = Conditional<(dim>0),
          BSplineSpace<dim-1,range,rank>,
          BSplineSpace<0,range,rank> >;

    using FaceSpace = PhysicalSpace<RefFaceSpace, typename PushForwardType::FacePushForward>;


    /** Type for the element accessor. */
    using ElementAccessor = BSplineElementAccessor<dim,range,rank>;

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
    ~BSplineSpace() = default;

protected:
    /** @name Constructors */
    ///@{
    /**
     * Constructs a maximum regularity BSpline space over CartesianGrid
     * @p knots for the given @p degree in all directions and homogeneous
     * in all components.
     */
    explicit BSplineSpace(const int degree, std::shared_ptr<GridType> knots);


    /**
     * Constructs a maximum regularity BSpline space over CartesianGrid
     * @p knots for the given @p degree[i] in the i-th direction and homogeneous
     * in all components.
     */
    explicit BSplineSpace(const TensorIndex<dim> &degree,
                          std::shared_ptr<GridType> knots);


    /**
     * Constructs a maximum regularity BSpline space over CartesianGrid
     * @p knots for the given @p degree for each direction and for each
     * component.
     */
    explicit BSplineSpace(const DegreeTable &degree,
                          std::shared_ptr<GridType> knots,
                          const bool homogeneous_range = false);


    /**
     * Constructs a BSpline space over the CartesianGrid
     * @p knots with the given multiplicity vector @p mult_vectors
     * for each component
     * and the given @p degree for each direction and for each
     * component.
     */
    explicit BSplineSpace(const DegreeTable &deg,
                          std::shared_ptr<GridType> knots,
                          std::shared_ptr<const MultiplicityTable> interior_mult,
                          const EndBehaviourTable &ends);
    ///@}


    /** @name Assignment operators */
    ///@{

    /** Copy assignment. Not allowed to be used. */
    self_t &
    operator=(const self_t &space) = delete;
    ///@}

public:
    /** @name Getting information about the space */
    ///@{
    std::vector<Index> get_loc_to_global(const TensorIndex<dim> &j) const;

    std::shared_ptr<const self_t >
    get_reference_space() const;
    ///@}

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

    /**
     * Prints internal information about the space.
     * @note Mostly used for debugging and testing.
     */
    void print_info(LogStream &out) const;





    std::shared_ptr<RefFaceSpace>
    get_ref_face_space(const Index face_id,
                       std::vector<Index> &face_to_element_dofs,
                       typename GridType::FaceGridMap &elem_map) const;


    std::shared_ptr<FaceSpace>
    get_face_space(const Index face_id,
                   std::vector<Index> &face_to_element_dofs) const;


    /** Return the push forward (non-const version). */
    std::shared_ptr<PushForwardType> get_push_forward();


    /** Return the push forward (const version). */
    std::shared_ptr<const PushForwardType> get_push_forward() const;


    /**
     * Adds an @p offset to the values of the dof ids.
     */
    void add_dofs_offset(const Index offset);



    /**
     * @note try not to use as plans are to make it private
     */
    TensorIndex<dim>
    basis_flat_to_tensor(const Index index, const Index comp) const;


    /**
     * @note try not to use as plans are to make it private
     */
    Index
    basis_tensor_to_flat(const TensorIndex<dim> &tensor_index,
                         const Index comp) const;


    std::shared_ptr<SpaceManager> get_space_manager();

    std::shared_ptr<const SpaceManager> get_space_manager() const;

private:

    /** Container with the local to global basis indices */
    DofDistribution<dim, range, rank> basis_indices_;

    /** @name Bezier extraction operator. */
    BernsteinExtraction<dim, range, rank> operators_;


    friend class BSplineElementAccessor<dim, range, rank>;


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


    /** Returns the container with the local to global basis indices (const version). */
    const DofDistribution<dim, range, rank> &
    get_basis_indices() const;

    /** Returns the container with the local to global basis indices (non-const version). */
    DofDistribution<dim, range, rank> &
    get_basis_indices();

};


IGA_NAMESPACE_CLOSE

#endif /* __BSPLINE_SPACE_H_ */
