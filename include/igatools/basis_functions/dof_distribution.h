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

#ifndef DOF_DISTRIBUTION_H_
#define DOF_DISTRIBUTION_H_

#include <igatools/base/config.h>
#include <igatools/utils/tensor_sized_container.h>
#include <igatools/basis_functions/spline_space.h>
#include <igatools/geometry/cartesian_grid_element_accessor.h>
#include <igatools/utils/concatenated_iterator.h>

IGA_NAMESPACE_OPEN

/**
 *
 * Class to handle the distribution of the basis function
 * indices storing what is known as the local to global map.
 *
 * It is a helper class for the BSplineSpace.
 *
 * @author pauletti, 2014
 *
 */
template<int dim, int range = 1, int rank = 1>
class DofDistribution : public TensorSizedContainer<dim>
{
public:
    using Space = SplineSpace<dim, range, rank>;
    using MultiplicityTable = typename Space::MultiplicityTable;
    using DegreeTable = typename Space::DegreeTable;
    using SpaceDimensionTable = typename Space::SpaceDimensionTable;
    using DofsPerElementTable = typename Space::template ComponentContainer<Index>;
    using IndexDistributionTable =
        StaticMultiArray<DynamicMultiArray<Index,dim>,range,rank>;


    /** Type alias for the dofs container used in each scalar component of a single-patch space. */
    using DofsComponentContainer = std::vector<Index>;

    /** Type alias for the View on the dofs in each scalar component of a single-patch space. */
    using DofsComponentView = ContainerView<DofsComponentContainer>;

    /** Type alias for the ConstView on the dofs in each scalar component of a single-patch space. */
    using DofsComponentConstView = ConstContainerView<DofsComponentContainer>;

    /** Type alias for a concatenated iterator defined on several component views. */
    using DofsIterator = ConcatenatedIterator<DofsComponentView>;

    /** Type alias for a concatenated const-iterator defined on several component views. */
    using DofsConstIterator = ConcatenatedConstIterator<DofsComponentView,DofsComponentConstView>;

    /** Type alias for the View on the dofs held by the single-patch space. */
    using DofsView = View<DofsIterator,DofsConstIterator>;

    /** Type alias for the ConstView on the dofs held by the single-patch space. */
    using DofsConstView = ConstView<DofsIterator,DofsConstIterator>;




    enum class DistributionPolicy
    {
        standard, component, other
    };

    DofDistribution() = delete;

    DofDistribution(std::shared_ptr<CartesianGrid<dim> > grid,
                    const MultiplicityTable &accum_mult,
                    const SpaceDimensionTable &n_basis,
                    const DegreeTable &degree_table,
                    DistributionPolicy pol = DistributionPolicy::standard);

    void reassign_dofs(const IndexDistributionTable &index_distribution, const DistributionPolicy pol);

    std::vector<Index> get_loc_to_global_indices(const TensorIndex<dim> &elem_tensor_id) const;

    std::vector<Index> get_loc_to_global_indices(const Index &elem_flat_id) const;


    TensorIndex<dim>
    basis_flat_to_tensor(const Index index, const Index comp) const;


    Index
    basis_tensor_to_flat(const TensorIndex<dim> &tensor_index,
                         const Index comp) const;

    /**
     * Print the class content
     */
    void print_info(LogStream &out) const;


    /** Add an @p offset to the dofs. */
    void add_dofs_offset(const Index offset);

    /** Returns the minimum dof id. */
    Index get_min_dof_id() const;

    /** Returns the maximum dof id. */
    Index get_max_dof_id() const;



    /** @name Getting information of a specific element */
    ///@{
    /**
     * Returns the number of dofs in each component, for each direction, relative to the element @p elem_flat_id.
     *
     * @note This make sense only for elements in which the dofs have (component-wise) a tensor-product structure.
     */
    DofsPerElementTable get_num_dofs_per_element_table(const Index elem_flat_id) const;

    /**
     * Returns the number of dofs of the element @p elem_flat_id.
     */
    Size get_num_dofs_per_element(const Index elem_flat_id) const;

    /**
     *  Return the number of dofs of the element @p elem_flat_id, for the i-th space component.
     */
    Size get_num_dofs_per_element(const Index elem_flat_id,const int i) const;
    ///@}

private:

    /**
     * Container used to store the dofs ids of each component of a single patch space.
     *
     * @warning This object can have a BIG memory footprint, therefore its copy is discouraged: please
     * use the associated View instead!
     */
    IndexDistributionTable index_distribution_;

    /**
     * View of the active dofs ids on a given single-patch space.
     */
    DofsView dofs_view_;



    void create_element_loc_to_global_view(
        std::shared_ptr<const CartesianGrid<dim> > grid,
        const MultiplicityTable &accum_mult,
        const SpaceDimensionTable &n_elem_basis);

    /**
     * Pointer to a vector in which each entry is a const view to the global dofs on a single element
     * of a space.
     * The size of the vector is equal to the number of active elements in the space.
     *
     * @note We use the pointer because this object can be used by other classes (@see SpaceManager),
     * and we want to keep the syncronization of the element views without the expense of successive copies.
     */
    std::shared_ptr<std::vector<DofsConstView>> elements_loc_to_global_flat_view_;

    DistributionPolicy policy_;

public:

    /**
     * Returns the container used to store the dofs ids of each component of a single patch space.
     *
     * @warning This object can have a BIG memory footprint, therefore its copy is discouraged: please
     * use the associated View instead!
     */
    const IndexDistributionTable &get_index_distribution() const;

    /**
     * Returns a view of the active dofs ids on a given single-patch space (non-const version).
     */
    DofsView &get_dofs_view();


    /**
     * Returns a view of the active dofs ids on a given single-patch space (const version).
     */
    const DofsView &get_dofs_view() const;



    /**
     * Returns a pointer to a vector in which each entry is a const view to the global dofs
     * on a single element of a space.
     * The size of the vector is equal to the number of active elements in the space.
     */
    std::shared_ptr<const std::vector<DofsConstView>> get_elements_view() const;
};

IGA_NAMESPACE_CLOSE

#endif
