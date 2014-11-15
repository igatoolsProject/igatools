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

#ifndef __NEW_PHYSICAL_SPACE_H_
#define __NEW_PHYSICAL_SPACE_H_

#include <igatools/base/config.h>
#include <igatools/geometry/new_mapping.h>
#include <igatools/geometry/new_push_forward.h>
#include <igatools/geometry/grid_forward_iterator.h>
#include <igatools/basis_functions/function_space.h>
#include <igatools/basis_functions/dof_distribution.h>
#include <igatools/utils/static_multi_array.h>

IGA_NAMESPACE_OPEN

class SpaceManager;

template <class> class PhysicalSpaceElement;

template <class> class SpaceElementHandler;

/**
 *
 * @sa FunctionSpace
 *
 * @ingroup containers
 */
template <class RefSpace_, int codim_, Transformation type_= Transformation::h_grad>
class NewPhysicalSpace :
    public std::enable_shared_from_this<NewPhysicalSpace<RefSpace_, codim_, type_>>,
            public FunctionSpaceOnGrid<CartesianGrid<RefSpace_::dim> >
{
private:
    using BaseSpace = FunctionSpaceOnGrid<CartesianGrid<RefSpace_::dim> >;
    using self_t = NewPhysicalSpace<RefSpace_, codim_, type_>;

public:
    ///@{
    /** see documentation in \ref FunctionSpaceOnGrid */
    using PushForwardType = NewPushForward<type_, RefSpace_::dim, codim_>;

    using RefSpace = RefSpace_;

    //using GridType = typename PushForwardType::GridType;
    ///@}
    using ElementHandler = SpaceElementHandler<self_t>;

    static const int dim = PushForwardType::dim;

    static const int codim = PushForwardType::codim;

    static const int space_dim = PushForwardType::space_dim;

    static const int range = PushForwardType::template PhysRange<RefSpace::range>::value;

    static const int rank = RefSpace::rank;

    using MapFunc =  MapFunction<dim, space_dim>;

    static constexpr int n_components = constexpr_pow(range, rank);

    static const std::array<int, n_components> components;

    /**
     * Type alias for the boundary conditions on each face of each scalar component of the space.
     */
    using BCTable = typename RefSpace::BCTable;

public:
    using Func = NewFunction<dim, codim, range, rank>;
    template <int order>
    using Derivative = typename Func::template Derivative<order>;
    using Point = typename Func::Point;
    using Value = typename Func::Value;
    using Gradient = typename Func::Gradient;
    using Div   = typename Func::Div;

    using RefPoint = typename RefSpace::Point;


public:
    template< class T>
    using ComponentContainer = typename RefSpace::template ComponentContainer<T>;

    using SpaceDimensionTable = typename RefSpace::SpaceDimensionTable;

    using DegreeTable = typename RefSpace::DegreeTable;

public:


    using ElementAccessor = PhysicalSpaceElement<self_t>;
    using ElementIterator = GridForwardIterator<ElementAccessor>;


    NewPhysicalSpace(const self_t &phys_space) = delete;

    static std::shared_ptr<self_t>
    create(std::shared_ptr<RefSpace> ref_space,
           std::shared_ptr<MapFunc> map_func);

    /**
     * Total number of dofs of the space.
     */
    Index get_num_basis() const;

    /**
     * Returns a element iterator to the first element of the patch.
     */
    ElementIterator begin() const;

    /**
     * Returns a element iterator to the last element of the patch.
     */
    ElementIterator last() const;

    /**
     * Returns a element iterator to one-pass the end of patch.
     */
    ElementIterator end() const;

    /** Returns the container with the global dof distribution (const version). */
    const DofDistribution<dim, range, rank> &get_dof_distribution_global() const;

    /** Returns the container with the global dof distribution (non const version). */
    DofDistribution<dim, range, rank> &get_dof_distribution_global();

    /** Returns the container with the patch dof distribution (const version). */
    const DofDistribution<dim, range, rank> &get_dof_distribution_patch() const;

    /** Returns the container with the patch dof distribution (non const version). */
    DofDistribution<dim, range, rank> &get_dof_distribution_patch();

    auto get_num_all_element_basis() const
    {
        return ref_space_->get_num_all_element_basis();
    }


#if 0
    const DegreeTable &get_degree() const;
#endif

    vector<Index> get_loc_to_global(const CartesianGridElement<dim> &element) const;

    vector<Index> get_loc_to_patch(const CartesianGridElement<dim> &element) const;


#if 0
    /**
     * Returns the element accessor with its flat id corresponding to @p elem_flat_id.
     *
     * @warning This function creates a new ElementAccessor object,
     * so it could be computationally expensive.
     */
    ElementAccessor get_element(const Index elem_flat_id) const;


    std::shared_ptr<const PushForwardType> get_push_forward() const;
#endif
    std::shared_ptr<const RefSpace> get_reference_space() const;
    std::shared_ptr<RefSpace> get_reference_space();
    std::shared_ptr<MapFunc> get_map_func() const
    {
        return map_func_;
    }

#if 0
    std::shared_ptr<FaceSpace>
    get_face_space(const Index face_id,
                   vector<Index> &face_to_element_dofs) const;

    Index get_id() const;


#endif
    void print_info(LogStream &out) const;



    std::shared_ptr<SpaceManager> get_space_manager();

    std::shared_ptr<const SpaceManager> get_space_manager() const;
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
        return ref_space_->get_boundary_conditions_table();
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
        return ref_space_->get_boundary_conditions_table();
    }
#endif

private:
    NewPhysicalSpace(std::shared_ptr<RefSpace> ref_space,
                     std::shared_ptr<MapFunc>  map_func);


    std::shared_ptr<RefSpace> ref_space_;
    std::shared_ptr<MapFunc>  map_func_;


    friend ElementAccessor;
};

IGA_NAMESPACE_CLOSE

#endif
