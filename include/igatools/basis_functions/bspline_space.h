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
#include <igatools/geometry/cartesian_grid.h>
#include <igatools/basis_functions/multiplicity.h>
#include <igatools/utils/static_multi_array.h>
#include <igatools/utils/dynamic_multi_array.h>
#include <igatools/basis_functions/function_space.h>
#include <igatools/linear_algebra/dense_matrix.h>
#include<igatools/geometry/mapping.h>
#include<igatools/geometry/push_forward.h>

IGA_NAMESPACE_OPEN

//Forward declaration to avoid including the header
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
 */
template<int dim_, int range_ = 1, int rank_ = 1>
class BSplineSpace :
    public std::enable_shared_from_this<BSplineSpace<dim_,range_,rank_> >,
    public FunctionSpaceOnGrid<CartesianGrid<dim_> >
{
private:
    using BaseSpace = FunctionSpaceOnGrid<CartesianGrid<dim_>>;

    /** Type for current class. */
    using self_t = BSplineSpace<dim_,range_,rank_>;

public:
    /** see documentation in \ref FunctionSpaceOnGrid */
    using PushForwardType = PushForward<Transformation::h_grad,dim_,0>;

    using RefSpace = self_t;

    using GridType = typename PushForwardType::GridType;

    static const int dim = PushForwardType::dim;

    static const int codim = PushForwardType::codim;

    static const int space_dim = PushForwardType::space_dim;

    static const int range = range_;

    static const int rank = rank_;

    static constexpr int n_components = constexpr_pow(range, rank);

    static const bool has_weights = false;

public:
    /** Type for the reference face space.*/
    using RefFaceSpace = BSplineSpace<dim-1,range,rank>;

    /** Type for the element accessor. */
    using ElementAccessor = BSplineElementAccessor<dim,range,rank>;

    /** Type for iterator over the elements.  */
    using ElementIterator = GridForwardIterator<ElementAccessor>;

public:
    /** Container indexed by the components of the space */
    template< class T>
    using ComponentTable = StaticMultiArray<T,range,rank>;

    using DegreeTable = ComponentTable<TensorIndex<dim>>;

    using MultiplicityTable = ComponentTable<Multiplicity<dim> >;

    /** @name Constructor and destructor */
    ///@{
    /**
     * Constructs a maximum regularity BSpline space over CartesianGrid
     * @p knots for the given @p degree in all directions and homogeneous
     * in all components.
     */
    explicit BSplineSpace(std::shared_ptr<GridType> knots, const int degree);

    /**
     * Smart pointer create construction technique, see more detail
     * in the corresponding wrapped constructor before.
     */
    static std::shared_ptr<self_t>
    create(std::shared_ptr<GridType> knots, const int degree);

    /**
     * Constructs a maximum regularity BSpline space over CartesianGrid
     * @p knots for the given @p degree[i] in the i-th direction and homogeneous
     * in all components.
     */
    explicit BSplineSpace(std::shared_ptr<GridType> knots,
                          const TensorIndex<dim> &degree);

    /**
     * Smart pointer create construction technique, see more detail
     * in the corresponding wrapped constructor before.
     */
    static std::shared_ptr<self_t>
    create(std::shared_ptr<GridType> knots, const TensorIndex<dim> &degree);

    /**
     * Constructs a maximum regularity BSpline space over CartesianGrid
     * @p knots for the given @p degree for each direction and for each
     * component.
     */
    explicit BSplineSpace(std::shared_ptr<GridType> knots,
                          const DegreeTable &degree);

    /**
     * Smart pointer create construction technique, see more detail
     * in the corresponding wrapped constructor before.
     */
    static std::shared_ptr<self_t>
    create(std::shared_ptr<GridType> knots,
           const DegreeTable &degree);


    /**
     * Constructs a BSpline space over the CartesianGrid
     * @p knots with the given multiplicity vector @p mult_vectors
     * for each component
     * and the given @p degree for each direction and for each
     * component.
     */
    explicit BSplineSpace(std::shared_ptr<GridType> knots,
                          const MultiplicityTable &mult_vectors,
                          const DegreeTable &degree);

    /**
     * Smart pointer create construction technique, see more detail
     * in the corresponding wrapped constructor before.
     */
    static std::shared_ptr<self_t>
    create(std::shared_ptr<GridType> knots,
           const MultiplicityTable &mult_vectors,
           const DegreeTable &degree);

    /** Destructor */
    ~BSplineSpace() = default;
    ///@}


    /** @name Assignment operators */
    ///@{

    /** Copy assignment. Not allowed to be used. */
    self_t &
    operator=(const self_t &space) = delete;
    ///@}

    /** @name Getting information about the space */
    ///@{
    /**
     * Returns true if all component belong to the same scalar valued
     * space.
     */
    bool is_range_homogeneous() const;

    /**
     * Total number of basis functions. This is the dimensionality
     * of the space.
     */
    Size get_num_basis() const;

    /**
     * Total number of basis functions
     * for the comp space component.
     */
    Size get_num_basis(const int comp) const;

    /**
     *  Total number of basis functions for the comp space component
     *  and the dir direction.
     */
    Size get_num_basis(const int comp, const int dir) const;

    /**
     * Component-direction indexed table with the number of basis functions
     * in each direction and component
     */
    DegreeTable get_num_basis_table() const;

    /**
     * Returns the number of dofs per element.
     */
    Size get_num_basis_per_element() const;

    /**
     *  Return the number of dofs per element for the i-th space component.
     */
    Size get_num_basis_per_element(int i) const;



    /**
     * Returns the degree of the BSpline space for each component and for each coordinate direction.
     * \return The degree of the BSpline space for each component and for each coordinate direction.
     * The first index of the returned object is the component id, the second index is the direction id.
     */
    const ComponentTable<TensorIndex<dim>> &get_degree() const;

    ///@}

    /** @name Getting the space data */
    ///@{
    /**
     * Return the knots with repetitions, in each direction, for each component of the space.
     */
    const ComponentTable<CartesianProductArray<Real,dim> > &
    get_knots_with_repetitions() const;

    std::shared_ptr<const self_t >
    get_reference_space() const;

    /**
     * Returns a reference to the dense multi array storing the global dofs.
     * Each element has a statically defined zone to read their dofs from,
     * independent of the distribution policy in use.
     */
    const ComponentTable<DynamicMultiArray<Index,dim>> &get_index_space() const;


    /**
     * @todo Missing documentation
     */
    const std::vector<std::vector<Index>> &get_element_global_dofs() const;
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

    /**
     * Transforms basis flat index of the component comp to a basis
     * tensor index.
     */
    TensorIndex<dim>
    flat_to_tensor(const Index index, const Index comp = 0) const;

    /**
     * Transforms a basis tensor index of the component comp to the
     * corresponding basis flat index.
     */
    Index
    tensor_to_flat(const TensorIndex<dim> &tensor_index,
                   const Index comp = 0) const ;
    ///@}

    /**
     * Prints internal information about the space.
     * @note Mostly used for debugging and testing.
     */
    void print_info(LogStream &out) const;




    /**
     * Return the knot multiplicities for each component of the space.
     */
    const ComponentTable<Multiplicity<dim> > &
    get_multiplicities() const;





    /** Return the push forward (non-const version). */
    std::shared_ptr<PushForwardType> get_push_forward();


    /** Return the push forward (const version). */
    std::shared_ptr<const PushForwardType> get_push_forward() const;



private:

    /**
     * @name Common code to initialize a bspline space
     */
    ///@{
    void init();
    void init_dofs();
    void fill_num_dof_per_element();
    void fill_bezier_extraction_operator();
    void fill_index_space_standard_policy();
    void fill_element_dofs_from_index_space();
    ///@}

    /**
     * Degree of the BSpline Space, they can be different for each direction
     * and for each component (see documentation on deg_table_t for index order).
     */
    const ComponentTable<TensorIndex<dim>> degree_;

    /**
     * Number of dofs per element.
     */
    int num_dofs_per_element_ = 0;


    /**
     * Multiplicities of the knots.
     */
    ComponentTable<Multiplicity<dim>> mult_;


    /**
     * Dense multi array storing the global dofs.
     * Each element has a statically defined zone to read their dofs from,
     * independent of the distribution policy in use.
     */
    ComponentTable<DynamicMultiArray<Index,dim> > index_space_;


    /** Where to read the global dofs of a given element */
    ComponentTable<Multiplicity<dim> > index_space_offset_;


    /**
     * Table container with one row per element, each row
     * contains the global dof indices corresponding to each element.
     * This is element_global_dofs_[elem_index][local_basis_indes]
     * gives the global dof index corresponding
     * to the local_basis_index-th on the element with index elem_index.
     * @todo: has to be removed, superseeded by index_space
     */
    using elem_dof_table_t = std::vector<std::vector<Index>>;
    elem_dof_table_t element_global_dofs_;


    /**
     * Number of one dimensional basis function in each component and in
     * each direction.
     * num_dofs_[comp][dir]
     */
    DegreeTable num_dofs_;


    /** @name Bezier extraction operator. */
    ///@{
    template<class T>
    using comp_p_array_table_t = ComponentTable<CartesianProductArray<T,dim>>;

    comp_p_array_table_t<DenseMatrix>         bezier_op_data_;
    comp_p_array_table_t<const DenseMatrix *> bezier_op_;
    ///@}


    /** @name Range optimization data */
    ///@{
protected:
    /**
     * True if each component of the vector valued space belongs
     * to the same scalar valued space.
     */
    const bool homogeneous_range_;

    //TODO(pauletti, Apr 27, 2014): make this private w/getter
public:
    /**
     * Knots (with repetitions) along each direction of each space component.
     */
    ComponentTable<CartesianProductArray<Real,dim>> knots_with_repetitions_;

private:

    /**
     * Stores information to allow an optimize treatment
     * of homogeneous range spaces.
     *
     */
    ComponentTable<int> map_component_to_active_data_;

    /**
     * Number of active components of the space.
     * If the space is range-homogeneous the number of active components is 1, otherwise
     * is n_components = pow(range,rank).
     */
    Size num_active_components_;
    ///@}


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


    //TODO(pauletti, Apr 27, 2014): make this private and use a getter
public:
    /** Knots with repetitions before refinement */
    ComponentTable<CartesianProductArray<Real,dim>> knots_with_repetitions_pre_refinement_;

public:
    DeclException1(ExcScalarRange, int,
                   << "Range " << arg1 << "should be 0 for a scalar valued"
                   << " space.");
};

IGA_NAMESPACE_CLOSE

#endif /* __BSPLINE_SPACE_H_ */
