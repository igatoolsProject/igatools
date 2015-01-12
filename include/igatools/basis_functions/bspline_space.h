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

#ifndef __BSPLINE_SPACE_H_
#define __BSPLINE_SPACE_H_

#include <igatools/base/config.h>
#include <igatools/base/logstream.h>

#include <igatools/basis_functions/reference_space.h>
#include <igatools/basis_functions/dof_distribution.h>
#include <igatools/basis_functions/bernstein_extraction.h>
#include <igatools/geometry/mapping.h>
#include <igatools/geometry/push_forward.h>
#include <igatools/basis_functions/physical_space.h>

IGA_NAMESPACE_OPEN

class SpaceManager;

template <int, int, int> class BSplineElement;
template <int, int, int> class BSplineElementHandler;
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
class BSplineSpace :
    public std::enable_shared_from_this<BSplineSpace<dim_,range_,rank_> >,
    public ReferenceSpace<dim_, range_, rank_>
{
private:
    using BaseSpace = ReferenceSpace<dim_, range_, rank_>;

    /** Type for current class. */
    using self_t = BSplineSpace<dim_,range_,rank_>;

public:
    /** see documentation in \ref FunctionSpaceOnGrid */

    using GridType = CartesianGrid<dim_>;
    using ElementHandler = BSplineElementHandler<dim_, range_, rank_>;

    static const int dim       = dim_;
    static const int codim     = 0;
    static const int space_dim = dim_;
    static const int range     = range_;
    static const int rank      = rank_;

//    using BaseSpace::n_components;
//    using BaseSpace::components;
//    using BaseSpace::dims;

public:
    using typename BaseSpace::Func;
    using typename BaseSpace::Point;
    using typename BaseSpace::Value;
    using typename BaseSpace::Derivative;
    using typename BaseSpace::Div;

    /**
     * See documentation in \ref FunctionSpaceOnGrid
     *
     * @see FunctionSpaceOnGrid
     */
    using PushForwardType = typename BaseSpace::PushForwardType;

    using RefSpace = typename BaseSpace::RefSpace;

    using RefPoint = Point;

public:
    /** Type for the element accessor. */
    using ElementAccessor = BSplineElement<dim,range,rank>;

    /** Type for iterator over the elements.  */
    using ElementIterator = CartesianGridIterator<ReferenceElement<dim,range,rank>>;


    using SpaceData = SplineSpace<dim_,range_,rank_>;

    using Degrees = typename SpaceData::Degrees;
    using Multiplicity = typename SpaceData::Multiplicity;
    using EndBehaviour = typename SpaceData::EndBehaviour;
    using Periodicity = typename SpaceData::Periodicity;

    using KnotsTable = typename SpaceData::KnotsTable;
    using DegreeTable = typename SpaceData::DegreeTable;
    using MultiplicityTable = typename SpaceData::MultiplicityTable;
    using SpaceDimensionTable = typename SpaceData::SpaceDimensionTable;
    using PeriodicTable = typename SpaceData::PeriodicTable;
    using EndBehaviourTable = typename SpaceData::EndBehaviourTable;

    using BCTable = typename SpaceData::BCTable;


    using BaseSpace::ComponentContainer;

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
    create(const int degree,
           std::shared_ptr<GridType> knots,
           const InteriorReg interior_reg = InteriorReg::maximum,
           const bool periodic = false,
           const BasisEndBehaviour end_b = BasisEndBehaviour::interpolatory);

    /**
     * Builds and returns a maximum regularity BSpline space over CartesianGrid
     * @p knots for the given @p degree[i] in the i-th direction and homogeneous
     * in all components.
     */
    static std::shared_ptr<self_t>
    create(const Degrees &degree,
           std::shared_ptr<GridType> knots,
           const InteriorReg interior_reg = InteriorReg::maximum,
           const Periodicity &periodic = filled_array<bool, dim>(false),
           const EndBehaviour &end_b =
               filled_array<BasisEndBehaviour, dim>(BasisEndBehaviour::interpolatory));

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
           const PeriodicTable &periodic,
           const EndBehaviourTable &end_b);
    ///@}

    /**
     * Create an element (defined on this space) with a given flat_index.
     */
    virtual std::shared_ptr<ReferenceElement<dim_,range_,rank_> > create_element(const Index flat_index) const override final;


    /** Destructor. */
    virtual ~BSplineSpace() = default;

protected:
    /** @name Constructors */
    ///@{
    /**
     * Constructs a maximum regularity BSpline space over CartesianGrid
     * @p knots for the given @p degree in all directions and homogeneous
     * in all components.
     */
    explicit BSplineSpace(const int degree,
                          std::shared_ptr<GridType> knots,
                          const InteriorReg interior_reg,
                          const bool periodic,
                          const BasisEndBehaviour end_b);

    /**
     * Constructs a maximum regularity BSpline space over CartesianGrid
     * @p knots for the given @p degree[i] in the i-th direction and homogeneous
     * in all components.
     */
    explicit BSplineSpace(const Degrees &degree,
                          std::shared_ptr<GridType> knots,
                          const InteriorReg interior_reg,
                          const Periodicity &periodic,
                          const EndBehaviour &end_b);

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
                          const PeriodicTable &periodic,
                          const EndBehaviourTable &end_b);


    /**
     * Copy constructor. Not allowed to be used.
     */
    BSplineSpace(const self_t &space) = delete;
    ///@}

    /** @name Assignment operators */
    ///@{
    /** Copy assignment. Not allowed to be used. */
    self_t &
    operator=(const self_t &space) = delete;
    ///@}

public:
    // TODO (pauletti, Oct 16, 2014): need to be documented or deleted, check!
    virtual vector<Index> get_loc_to_global(const CartesianGridElement<dim> &element) const override final;

    virtual vector<Index> get_loc_to_patch(const CartesianGridElement<dim> &element) const override final;

    /** @name Functions involving the element iterator */
    ///@{
    /**
     * Returns a element iterator to the first element of the patch
     */
    virtual ElementIterator begin() const override final;

    /**
     * Returns a element iterator to the last element of the patch
     */
    virtual ElementIterator last() const override final;


    /**
     * Returns a element iterator to one-pass the end of patch.
     */
    virtual ElementIterator end() const override final;
    ///@}

#if 0
    ElementHandler get_element_handler() const;
#endif

    template <int k>
    using InterGridMap = typename GridType::template InterGridMap<k>;

    template <int k>
    using InterSpaceMap = vector<Index>;

    template <int k>
    using SubRefSpace = ReferenceSpace<k, range, rank>;

    template <int k>
    using SubSpace = PhysicalSpace<k,range,rank,dim-k, Transformation::h_grad>;

    /**
     * Construct a sub space of dimension k conforming to
     * the subspace sub element sub_elem_id and a map from the elements of
     * the sub_element grid to the corresponding element of the current
     * grid.
     */
    template<int k>
    std::shared_ptr<SubRefSpace<k> >
    get_ref_sub_space(const int sub_elem_id,
                      InterSpaceMap<k> &dof_map,
                      std::shared_ptr<CartesianGrid<k>> sub_grid = nullptr) const;

    template<int k>
    std::shared_ptr<SubSpace<k> >
    get_sub_space(const int s_id, InterSpaceMap<k> &dof_map,
                  std::shared_ptr<CartesianGrid<k>> sub_grid,
                  std::shared_ptr<typename GridType::template InterGridMap<k>> elem_map) const;

    std::shared_ptr<const self_t > get_reference_space() const;


    virtual std::shared_ptr<SpaceManager> get_space_manager() override final ;
    virtual std::shared_ptr<const SpaceManager> get_space_manager() const override final;


    /** Returns the container with the global dof distribution (const version). */
    virtual const DofDistribution<dim, range, rank> &
    get_dof_distribution_global() const override final;

    /** Returns the container with the global dof distribution (non const version). */
    virtual DofDistribution<dim, range, rank> &
    get_dof_distribution_global() override final;

    /** Returns the container with the patch dof distribution (const version). */
    virtual const DofDistribution<dim, range, rank> &
    get_dof_distribution_patch() const override final;


    /** Returns the container with the patch dof distribution (non const version). */
    virtual DofDistribution<dim, range, rank> &
    get_dof_distribution_patch() override final;
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
    virtual void print_info(LogStream &out) const override final;

private:


    /** Container with the local to global basis indices
     * @note The concept of global indices refers to a global numeration of the
     * dofs of all the spaces.
     */
    EndBehaviourTable end_b_;

    DofDistribution<dim, range, rank> dof_distribution_global_;

    /** Container with the local to patch basis indices
     * @note The concept of patch indices refers to the numeration at patch
     * level of the dofs.
     */
    DofDistribution<dim, range, rank> dof_distribution_patch_;

    /** @name Bezier extraction operator. */
    BernsteinExtraction<dim, range, rank> operators_;


    /** If end knots are not in the repeated knot vector */
    using EndIntervalTable = typename BaseSpace::template
                             ComponentContainer<std::array<std::pair<Real, Real>, dim>>;
    EndIntervalTable end_interval_;

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

    void create_connection_for_h_refinement(std::shared_ptr<self_t> space);

public:
    DeclException1(ExcScalarRange, int,
                   << "Range " << arg1 << "should be 0 for a scalar valued"
                   << " space.");


    virtual bool is_bspline() const override final
    {
        return true;
    }

    virtual std::shared_ptr<typename ReferenceSpace<dim_,range_,rank_>::ElementHandler> create_elem_handler() const override final
    {
        const auto this_space = std::enable_shared_from_this<self_t>::shared_from_this();
        return ElementHandler::create(this_space);
    }

};

IGA_NAMESPACE_CLOSE

#endif /* __BSPLINE_SPACE_H_ */
