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

#ifndef NURBS_SPACE_H_
#define NURBS_SPACE_H_

#include <igatools/base/config.h>

#include <igatools/basis_functions/bspline_space.h>
#include <igatools/basis_functions/bspline_element.h>
#include <igatools/base/ig_function.h>

#ifdef NURBS

IGA_NAMESPACE_OPEN

template <int, int, int> class NURBSElement;
template <int, int, int> class NURBSElementHandler;

/**
 * Multivariate (tensor product) scalar, vector or k-tensor
 * valued NURBS space.
 *
 * @ingroup containers
 */
template <int dim_, int range_ = 1, int rank_ = 1>
class NURBSSpace :
    public std::enable_shared_from_this<NURBSSpace<dim_,range_,rank_> >,
    public ReferenceSpace<dim_,range_,rank_>
{
private:
    using BaseSpace = ReferenceSpace<dim_,range_,rank_>;
    using self_t = NURBSSpace<dim_, range_, rank_>;

public:
    using SpSpace = BSplineSpace<dim_, range_, rank_>;

    /** see documentation in \ref FunctionSpaceOnGrid */

    using GridType = CartesianGrid<dim_>;
    static const int dim       = dim_;
    static const int codim     = 0;
    static const int space_dim = dim_;
    static const int range     = range_;
    static const int rank      = rank_;
    static const bool is_physical_space = false;

    static const auto n_components = SplineSpace<dim_, range_, rank_>::n_components;


    /**
     * See documentation in \ref FunctionSpaceOnGrid
     *
     * @see FunctionSpaceOnGrid
     */
    using PushForwardType = typename BaseSpace::PushForwardType;

    using RefSpace = typename BaseSpace::RefSpace;



public:
    using Func = typename SpSpace::Func;
    template <int order>
    using Derivative = typename SpSpace::template Derivative<order>;
    using Point = typename SpSpace::Point;
    using Value = typename SpSpace::Value;
    using Div   = typename SpSpace::Div;

    using RefPoint = typename SpSpace::RefPoint;

public:

    /** Type for the element accessor. */
    using ElementAccessor = NURBSElement<dim, range, rank> ;

    /** Type for iterator over the elements.  */
    using ElementIterator = CartesianGridIterator<ReferenceElement<dim,range,rank> >;

    using ElementHandler = NURBSElementHandler<dim_, range_, rank_>;



    template <int k>
    using InterGridMap = typename GridType::template InterGridMap<k>;

    template <int k>
    using InterSpaceMap = vector<Index>;

    template <int k>
    using SubRefSpace = NURBSSpace<k, range, rank>;

    template <int k>
    using SubSpace = PhysicalSpace<k,range,rank, dim-k, Transformation::h_grad>;

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

public:
//    /** Container indexed by the components of the space */
    template< class T>
    using ComponentContainer = typename SpSpace::template ComponentContainer<T>;

    using Degrees = typename SpSpace::Degrees;
    using Multiplicity = typename SpSpace::Multiplicity;
    using EndBehaviour = typename SpSpace::EndBehaviour;
    using Periodicity = typename SpSpace::Periodicity;

    using KnotsTable = typename SpSpace::KnotsTable;
    using DegreeTable = typename SpSpace::DegreeTable;
    using MultiplicityTable = typename SpSpace::MultiplicityTable;
    using TensorSizeTable = typename SpSpace::TensorSizeTable;
    using PeriodicityTable = typename SpSpace::PeriodicityTable;
    using EndBehaviourTable = typename SpSpace::EndBehaviourTable;

    using BCTable = typename SpSpace::BCTable;

    using WeightSpace = BSplineSpace<dim_,1,1>;
    using WeightFunction = IgFunction<ReferenceSpace<dim_,1,1> >;
    using WeightFunctionPtr = std::shared_ptr<WeightFunction>;
    using WeightFunctionPtrTable = ComponentContainer<WeightFunctionPtr>;
    using Weights = DynamicMultiArray<Real,dim>;
    using WeightsTable = ComponentContainer<Weights>;

public:
    /**
     * Returns a shared_ptr wrapping a NURBSSpace from a BSplineSpace and a scalar weight function.
     */
    static std::shared_ptr<self_t>
    create(std::shared_ptr<SpSpace> bs_space,
           const WeightFunctionPtrTable &weight_func_table);

    ///@}

    /**
     * Create an element (defined on this space) with a given flat_index.
     */
    virtual std::shared_ptr<ReferenceElement<dim_,range_,rank_> > create_element(const Index flat_index) const override final;

    /** Destructor */
    virtual ~NURBSSpace() = default;


protected:
    /** @name Constructor */
    ///@{
    /**
     * Construct a NURBSSpace from a BSplineSpace and a table of weights.
     */
    explicit NURBSSpace(std::shared_ptr<SpSpace> bs_space,
                        const WeightFunctionPtrTable &weight_func_table);

    /**
     * Copy constructor. Not allowed to be used.
     */
    NURBSSpace(const self_t &space) = delete;

    ///@}

public:
    /** @name Getting information about the space */
    ///@{

    virtual bool is_bspline() const override final
    {
        return false;
    }

    virtual const DegreeTable &get_degree() const override final
    {
        return this->sp_space_->get_degree();
    }

    virtual vector<Index> get_element_dofs(
        const CartesianGridElement<dim> &element,
        const std::string &dofs_property = DofProperties::none) const override final
    {
        return this->sp_space_->get_element_dofs(element,dofs_property);
    }


#if 0
    /**
     * Returns the multiplicity of the internal knots that defines the BSpline
     * space for each component and for each coordinate direction.
     * \return The multiplicity of the internal knots that defines the BSpline
     * space for each component and for each coordinate direction.
     */
    std::shared_ptr<const MultiplicityTable> get_interior_mult() const;
#endif

    ///@}

    const std::shared_ptr<SpSpace> get_spline_space() const;


#if 0
    /** Return the push forward (non-const version). */
    std::shared_ptr<PushForwardType> get_push_forward();

    /** Return the push forward (const version). */
    std::shared_ptr<const PushForwardType> get_push_forward() const;

    std::shared_ptr<const self_t >
    get_reference_space() const;
#endif


    /**
     * Get the weights of the NURBSSpace.
     */
    const WeightFunctionPtrTable &get_weights() const;

    const PeriodicityTable &get_periodicity() const override final
    {
        return sp_space_->get_periodicity();
    }

#if 0
    /**
     * Returns a reference to the end behaviour table of the BSpline space.
     */
    virtual EndBehaviourTable &get_end_behaviour_table() override final
    {
        return sp_space_->get_end_behaviour_table();
    };
#endif
    /**
     * Returns a const reference to the end behaviour table of the BSpline space.
     */
    virtual const EndBehaviourTable &get_end_behaviour_table() const override final
    {
        return sp_space_->get_end_behaviour_table();
    };


    /**
     * Prints internal information about the space.
     * @note Mostly used for debugging and testing.
     */
    virtual void print_info(LogStream &out) const override final;

#if 0
    /**
     * Returns a const-reference to the table containing
     * the boundary conditions on each face of each scalar component of the space.
     *
     * For example, with the code
     * @code{.cpp}
       const auto &bc_table = space.get_boundary_conditions_table();

       BoundaryConditionType bc_id = bc_table[1][3]; // boundary condition on face 3 of space's component 1
       @endcode
     * we copy to the variable <tt>bc_id</tt> the value of the boundary condition
     * on the face 3 of the space component 1.
     *
     * @sa BoundaryConditionType
     */
    const BCTable &get_boundary_conditions_table() const
    {
        return sp_space_->get_boundary_conditions_table();
    }

    /**
     * Returns a reference to the table containing
     * the boundary conditions on each face of each scalar component of the space.
     *
     * For example, with the code
     * @code{.cpp}
       auto &bc_table = space.get_boundary_conditions_table();

       bc_table[1][3] = BoundaryConditionType::DirichletHomogeneous; // setting Dirichlet homogeneous boundary condition on face 3 of space's component 1
       @endcode
     * we assign the value <tt>BoundaryConditionType::DirichletHomogeneous</tt> to the
     * boundary condition on the face 3 of the space component 1.
     *
     * @sa BoundaryConditionType
     */
    BCTable &get_boundary_conditions_table()
    {
        return sp_space_->get_boundary_conditions_table();
    }
#endif
private:
    /**
     * B-spline space
     */
    std::shared_ptr<SpSpace> sp_space_;

    /**
     * Weight function component table. Each component of the table is a shared_ptr to an IgFunction
     * built using a scalar BSplineSpace.
     */
    WeightFunctionPtrTable weight_func_table_;


#if 0
    /**
     * Refines the NURBSSpace after the uniform refinement of the BSplineSpace.
     *
     * @param[in] refinement_directions Directions along which the refinement is performed.
     * @param[in] grid_old Grid before the refinement.
     *
     * @pre Before invoking this function, must be invoked the function BSplineSpace::refine().
     * @note This function is connected to the CartesianGrid's signal for the refinement, and
     * it is necessary in order to avoid infinite loops in the CartesianGrid::refine() function calls.
     * @note The implementation of this function is based on "The NURBS Book" Algorithm A5.4.
     *
     * @ingroup h_refinement
     */
    void refine_h_weights(const std::array<bool,dim> &refinement_directions,
                          const GridType &grid_old);

    /**
     * Create a signal and a connection for the refinement.
     */
    void create_refinement_connection();
    /**
     * Performs checks after the construction of the object.
     * In debug mode, if something is going wrong, an assertion will be raised.
     *
     * @warning This function should be used as last line in the implementation of each constructor.
     */
    void perform_post_construction_checks() const;


#endif
    friend ElementAccessor;
    friend ElementHandler;

    /**
     * Returns the weight coefficient associated with a given basis function.
     */
    Real get_weight_coef_from_basis_id(const Index basis_id) const;


public:
    virtual std::shared_ptr<typename BaseSpace::ElementHandler> create_elem_handler() const override final
    {
        const auto this_space = std::enable_shared_from_this<self_t>::shared_from_this();
        return ElementHandler::create(this_space);
    }

};

IGA_NAMESPACE_CLOSE
#endif /* #ifdef NURBS */

#endif /* NURBS_SPACE_H_ */


