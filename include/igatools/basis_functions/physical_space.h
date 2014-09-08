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

#ifndef __PHYSICAL_SPACE_H_
#define __PHYSICAL_SPACE_H_

#include <igatools/base/config.h>
#include <igatools/base/function.h>
#include <igatools/geometry/mapping.h>
#include <igatools/geometry/push_forward.h>
#include <igatools/geometry/grid_forward_iterator.h>
#include <igatools/basis_functions/function_space.h>
#include <igatools/basis_functions/dof_distribution.h>
#include <igatools/utils/static_multi_array.h>

IGA_NAMESPACE_OPEN


class SpaceManager;

//Forward declaration to avoid including the header
template < class > class PhysicalSpaceElementAccessor;

/**
 *
 * @sa FunctionSpace
 *
 * @ingroup containers
 */
template <class RefSpace_, class PushForward_>
class PhysicalSpace :
    public std::enable_shared_from_this<PhysicalSpace<RefSpace_, PushForward_>>,
            public FunctionSpaceOnGrid<CartesianGrid<RefSpace_::dim> >
{
private:
    using BaseSpace = FunctionSpaceOnGrid<CartesianGrid<RefSpace_::dim> >;
    using self_t = PhysicalSpace<RefSpace_, PushForward_>;

public:
    ///@{
    /** see documentation in \ref FunctionSpaceOnGrid */
    using PushForwardType = PushForward_;

    using RefSpace = RefSpace_;

    using GridType = typename PushForwardType::GridType;
    ///@}
    static const int dim = PushForwardType::dim;

    static const int codim = PushForwardType::codim;

    static const int space_dim = PushForwardType::space_dim;

    static const int range = PushForwardType::template PhysRange<RefSpace::range>::value;

    static const int rank = RefSpace::rank;

    static constexpr int n_components = constexpr_pow(range, rank);

public:
    using Func = Function<space_dim, range, rank>;
    template <int order>
    using Derivative = typename Func::template Derivative<order>;
    using Point = typename Func::Point;
    using Value = typename Func::Value;
    using Div   = typename Func::Div;

    using RefPoint = typename RefSpace::Point;


public:
    template< class T>
    using ComponentContainer = typename RefSpace::template ComponentContainer<T>;

    using SpaceDimensionTable = typename RefSpace::SpaceDimensionTable;

    using DegreeTable = typename RefSpace::DegreeTable;

public:
    /** Type for the reference space on the face. */
    using RefFaceSpace = typename RefSpace_::RefFaceSpace;
    using FaceSpace = PhysicalSpace<RefFaceSpace, typename PushForwardType::FacePushForward>;
    /**
     * Type for the element accessor.
     */
    using ElementAccessor = PhysicalSpaceElementAccessor<self_t>;

    /**
     * Typedef for the element iterator
     */
    typedef GridForwardIterator<ElementAccessor> ElementIterator;


    PhysicalSpace(std::shared_ptr<RefSpace> ref_space,
                  std::shared_ptr<PushForwardType> push_forward);

    PhysicalSpace(const self_t &phys_space) = delete;

    static std::shared_ptr<self_t> create(
        std::shared_ptr<RefSpace> ref_space,
        std::shared_ptr<PushForwardType> push_forward);

    /**
     * Total number of dofs of the space.
     */
    Index get_num_basis() const;

    /** Returns the container with the global dof distribution (const version). */
    const DofDistribution<dim, range, rank> &get_dof_distribution_global() const;

    /** Returns the container with the global dof distribution (non const version). */
    DofDistribution<dim, range, rank> &get_dof_distribution_global();

    /** Returns the container with the patch dof distribution (const version). */
    const DofDistribution<dim, range, rank> &get_dof_distribution_patch() const;

    /** Returns the container with the patch dof distribution (non const version). */
    DofDistribution<dim, range, rank> &get_dof_distribution_patch();

    SpaceDimensionTable get_num_all_element_basis() const
    {
    	return ref_space_->get_num_all_element_basis();
    }

    const DegreeTable &get_degree() const;

    vector<Index> get_loc_to_global(const CartesianGridElement<dim> &element) const;

    vector<Index> get_loc_to_patch(const CartesianGridElement<dim> &element) const;

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

    /**
     * Returns the element accessor with its flat id corresponding to @p elem_flat_id.
     *
     * @warning This function creates a new ElementAccessor object,
     * so it could be computationally expensive.
     */
    ElementAccessor get_element(const Index elem_flat_id) const;


    std::shared_ptr<const PushForwardType> get_push_forward() const;

    std::shared_ptr<const RefSpace> get_reference_space() const;


    std::shared_ptr<FaceSpace>
    get_face_space(const Index face_id,
                   vector<Index> &face_to_element_dofs) const;


    void print_info(LogStream &out) const;



    Index get_id() const;

    // TODO (pauletti, Jun 12, 2014): if we are using this it should be
    // implemented in all library classes
    void print_memory_info(LogStream &out) const;


    std::shared_ptr<SpaceManager> get_space_manager();

    std::shared_ptr<const SpaceManager> get_space_manager() const;




private:
    std::shared_ptr<RefSpace> ref_space_;

    std::shared_ptr<PushForwardType> push_forward_;

    friend ElementAccessor;
};

IGA_NAMESPACE_CLOSE

#endif
