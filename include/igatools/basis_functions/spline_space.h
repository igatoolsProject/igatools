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

#ifndef SPLINE_SPACE_H_
#define SPLINE_SPACE_H_

#include <igatools/base/config.h>
#include <igatools/basis_functions/reference_space.h>

IGA_NAMESPACE_OPEN







/**
 * @brief Tensor product spline space
 *
 * A one dimensional polynomial spline space is determined by:
 * - a domain, the interval [a,b]
 * - the polynomial order
 * - a partition of [a,b], the knots
 * - the interior knots smoothness (given by the interior multiplicity)
 *
 * This class provides the realization of a scalar, vector or k-tensor
 * value tensor product spline space.
 *
 * @note This is independent of the basis functions one may wish to use
 * for the given space.
 *
 * @author pauletti, 2014
 *
 */
template<int dim, int range = 1, int rank = 1>
class SplineSpace :
    public ReferenceSpace<dim,range,rank>
{

private:
    using RefSpace = ReferenceSpace<dim,range,rank>;
    using GridSpace = FunctionSpaceOnGrid<CartesianGrid<dim>>;
    using typename GridSpace::GridType;


public:
    using GridSpace::dims;




    using Func = Function<dim, 0, range, rank>;

public:
    template <int order>
    using Derivative = typename Func::template Derivative<order>;
    using Point = typename Func::Point;
    using Value = typename Func::Value;
    using Div   = typename Func::Div;

public:
    using RefSpace::n_components;
//    static constexpr int n_components = RefSpace::template ComponentContainer<Size>::n_entries;
//    static const std::array<Size, n_components> components;


public:
    using KnotCoordinates = typename GridType::KnotCoordinates;
    using BoundaryKnots = std::array<CartesianProductArray<Real,2>, dim>;
    using Multiplicity = CartesianProductArray<Size, dim>;

    using DegreeTable = typename RefSpace::DegreeTable;

    using MultiplicityTable = typename RefSpace::template ComponentContainer<Multiplicity>;
    using BoundaryKnotsTable = typename RefSpace::template ComponentContainer<BoundaryKnots>;
    using KnotsTable = typename RefSpace::template ComponentContainer<KnotCoordinates>;
    using PeriodicTable = typename RefSpace::template ComponentContainer<std::array<bool, dim> >;

    using IndexSpaceTable = typename RefSpace::template ComponentContainer<DynamicMultiArray<Index,dim>>;
    using IndexSpaceMarkTable = Multiplicity;


    using SpaceDimensionTable = typename RefSpace::SpaceDimensionTable;

    using BCTable = typename RefSpace::BCTable;

    enum class EndBehaviour
    {
        /**
         * Interpolatory basis functions at knots bounday (i.e. open knot vector).
         */
        interpolatory,

        periodic,

        end_knots
    };
    using EndBehaviourTable = typename RefSpace::template ComponentContainer<std::array<EndBehaviour, dim> >;

    // For the interior multiplicities
    // maximum regularity
    // minimul regularity discontinous
    enum class InteriorReg
    {
        maximum, minimun
    };

public:
    /**
     * Most general constructor
     */
    explicit SplineSpace(const DegreeTable &deg,
                         std::shared_ptr<GridType> knots,
                         std::shared_ptr<const MultiplicityTable> interior_mult,
                         const PeriodicTable &periodic = PeriodicTable(filled_array<bool,dim>(false)));

    explicit SplineSpace(const DegreeTable &deg,
                         std::shared_ptr<GridType> knots,
                         const InteriorReg &interior_mult,
                         const PeriodicTable &periodic = PeriodicTable(filled_array<bool,dim>(false)))
        :SplineSpace(deg, knots, fill_max_regularity(deg, knots), periodic)
    {}

    const DegreeTable &get_degree() const
    {
        return deg_;
    }


    virtual const std::array<Index,n_components> &get_components_map() const override final
    {
        return interior_mult_->get_comp_map();
    }

    /** @name Getting information about the space */
    ///@{
    /**
     * Total number of basis functions. This is the dimensionality
     * of the space.
     */
    virtual Size get_num_basis() const override final
    {
        return space_dim_.total_dimension;
    }

    /**
     * Total number of basis functions
     * for the comp space component.
     */
    virtual Size get_num_basis(const int comp) const override final
    {
        return space_dim_.comp_dimension[comp];
    }

    /**
     *  Total number of basis functions for the comp space component
     *  and the dir direction.
     */
    virtual Size get_num_basis(const int comp, const int dir) const override final
    {
        return  space_dim_[comp][dir];
    }

    /**
     * Component-direction indexed table with the number of basis functions
     * in each direction and component
     */
    virtual const SpaceDimensionTable &get_num_basis_table() const override final
    {
        return space_dim_;
    }

    virtual SpaceDimensionTable get_num_all_element_basis() const override
    {
        typename RefSpace::template ComponentContainer<TensorSize<dim>> n_basis(deg_.get_comp_map());
        for (auto comp : deg_.get_active_components_id())
            n_basis[comp] = TensorSize<dim>(deg_[comp]+1);

        return SpaceDimensionTable(n_basis);
    }

    /**
     * Component table with the offset of basis functions
     * in each component.
     */
    typename RefSpace::template ComponentContainer<Size> get_basis_offset() const
    {
        typename RefSpace::template ComponentContainer<Size> offset;
        offset[0] = 0;
        for (int comp = 1; comp < n_components; ++comp)
            offset[comp] = offset[comp-1] + space_dim_.comp_dimension[comp];

        return offset;
    }

    ///@}

    template<int k>
    using SubSpace = SplineSpace<k, range, rank>;

    template<int k>
    std::shared_ptr<typename SubSpace<k>::MultiplicityTable>
    get_sub_space_mult(const Index s_id) const;

    template<int k>
    typename SubSpace<k>::DegreeTable
    get_sub_space_degree(const Index s_id) const;


#if 0
    /**
     * Returns the multiplicity of the face space face_id
     */
    std::shared_ptr<typename FaceSpace::MultiplicityTable>
    get_face_mult(const Index face_id) const;

    /**
     * Returns the multiplicity of the face space face_id
     */
    typename FaceSpace::DegreeTable
    get_face_degree(const Index face_id) const;
#endif

    KnotsTable compute_knots_with_repetition(const BoundaryKnotsTable &boundary_knots) const;

    KnotsTable compute_knots_with_repetition(const EndBehaviourTable &ends) const;

    /**
     * For each element and for each component there is an initial
     * tensor index in the Index space from where all non-zero basis
     * function can be determined.
     */
    MultiplicityTable accumulated_interior_multiplicities() const;




    virtual void print_info(LogStream &out) const override;


private:
    /**
     * Fill the multiplicy for the maximum possible regularity
     *  of the given number of knots
     */
    std::shared_ptr<MultiplicityTable> fill_max_regularity(const DegreeTable &deg, std::shared_ptr<const GridType> grid);

#if 0
    BoundaryKnotsTable interpolatory_end_knots() const;
#endif

    CartesianProductArray<Real,2> interpolatory_end_knots(const int comp_id,const int dir) const;


private:
    std::shared_ptr<const MultiplicityTable> interior_mult_;

    DegreeTable deg_;

    /** Table with the dimensionality of the space in each component and direction */
    SpaceDimensionTable space_dim_;

    EndBehaviourTable end_behaviour_;

    PeriodicTable periodic_;



    /**
     * Boundary conditions on each face of each scalar component of the space.
     */
    BCTable boundary_conditions_table_;

public:

    /** Returns the multiplicity of the internal knots that defines the space. */
    std::shared_ptr<const MultiplicityTable> get_interior_mult() const
    {
        return interior_mult_;
    }

    const EndBehaviourTable &get_end_behaviour() const
    {
        return end_behaviour_;
    }

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
        return boundary_conditions_table_;
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
        return boundary_conditions_table_;
    }



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

    std::shared_ptr<const SplineSpace<dim,range,rank> > spline_space_previous_refinement_;

public:
    std::shared_ptr<const SplineSpace<dim,range,rank> >
    get_spline_space_previous_refinement() const
    {
        return spline_space_previous_refinement_;
    }

protected:

    /** This function initialize the member variables from the constructor
     * arguments or after an h-refinement. */
    void init();

public:
    virtual bool is_bspline() const override
    {
        Assert(false,ExcMessage("This class should not have this function."))
        return true;
    }

    virtual vector<Index> get_loc_to_global(const CartesianGridElement<dim> &element) const override
    {
        Assert(false,ExcMessage("This class should not have this function."))
        return vector<Index>();
    }

    virtual vector<Index> get_loc_to_patch(const CartesianGridElement<dim> &element) const override
    {
        Assert(false,ExcMessage("This class should not have this function."))
        return vector<Index>();
    }

    using ElementIterator = CartesianGridIterator<ReferenceElement<dim,range,rank>>;

    virtual ElementIterator begin() const override
    {
        Assert(false,ExcMessage("This class should not have this function."));
        return ElementIterator(nullptr, 0);
    }

    virtual ElementIterator end() const override
    {
        Assert(false,ExcMessage("This class should not have this function."));
        return ElementIterator(nullptr, 0);
    }

    virtual ElementIterator last() const override
    {
        Assert(false,ExcMessage("This class should not have this function."));
        return ElementIterator(nullptr, 0);
    }


    /** Returns the container with the global dof distribution (const version). */
    virtual const DofDistribution<dim, range, rank> &
    get_dof_distribution_global() const override
    {
        Assert(false,ExcMessage("This class should not have this function."));
        return *reinterpret_cast<const DofDistribution<dim,range,rank> *>(this);
    }

    /** Returns the container with the global dof distribution (non const version). */
    virtual DofDistribution<dim, range, rank> &
    get_dof_distribution_global() override
    {
        Assert(false,ExcMessage("This class should not have this function."));
        return *reinterpret_cast<DofDistribution<dim,range,rank> *>(this);
    }

    /** Returns the container with the patch dof distribution (const version). */
    virtual const DofDistribution<dim, range, rank> &
    get_dof_distribution_patch() const override
    {
        Assert(false,ExcMessage("This class should not have this function."));
        return *reinterpret_cast<const DofDistribution<dim,range,rank> *>(this);
    }


    /** Returns the container with the patch dof distribution (non const version). */
    virtual DofDistribution<dim, range, rank> &
    get_dof_distribution_patch() override
    {
        Assert(false,ExcMessage("This class should not have this function."));
        return *reinterpret_cast<DofDistribution<dim,range,rank> *>(this);
    }


public:
    virtual std::shared_ptr<typename RefSpace::ElementHandler> create_elem_handler() const override
    {
        Assert(false,ExcMessage("This class should not have this function."));
        return nullptr;
    }

};









IGA_NAMESPACE_CLOSE

#endif
