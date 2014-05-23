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
#include <igatools/geometry/mapping.h>
#include <igatools/geometry/push_forward.h>
#include <igatools/geometry/grid_forward_iterator.h>
#include <igatools/basis_functions/function_space.h>
#include <igatools/utils/static_multi_array.h>

IGA_NAMESPACE_OPEN

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


    /** Container indexed by the components of the space */
    template< class T>
    using ComponentTable = StaticMultiArray<T,range,rank>;

public:
    /** Type for the reference space on the face. */
    using RefFaceSpace = typename RefSpace_::RefFaceSpace;

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

    std::shared_ptr<self_t> clone() const;

    static std::shared_ptr<self_t> create(
        std::shared_ptr<RefSpace> ref_space,
        std::shared_ptr<PushForwardType> push_forward);

    /**
     * Total number of dofs of the space.
     */
    Index get_num_basis() const;

    /**
     * Returns the number of dofs per element.
     */
    int get_num_basis_per_element() const;

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

    std::shared_ptr<const PushForwardType> get_push_forward() const;

    std::shared_ptr<const RefSpace> get_reference_space() const;

    void print_info(LogStream &out) const;

    void print_memory_info(LogStream &out) const;


    /**
     * Returns the degree of the BSpline space for each component and for each coordinate direction.
     * The first index of the returned object is the component id, the second index is the direction id.
     */
    const ComponentTable<TensorIndex<dim>> &get_degree() const;


    /**
     * @todo Missing documentation
     */
    const std::vector<std::vector<Index>> &get_element_global_dofs() const;

private:
    std::shared_ptr<RefSpace> ref_space_;

    std::shared_ptr<PushForwardType> push_forward_;

    /**
     * Map between the element accessors defined by the reference space and the element accessors defined by the push-forward.
     */
    std::vector<int> map_elements_;


    friend ElementAccessor;
};

IGA_NAMESPACE_CLOSE

#endif
