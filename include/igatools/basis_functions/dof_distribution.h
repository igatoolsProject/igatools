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
#include <igatools/geometry/cartesian_grid_element.h>
#include <igatools/utils/concatenated_iterator.h>

IGA_NAMESPACE_OPEN

/**
 *
 * Class to handle the distribution of the basis function
 * indices that are defined on a single patch space,
 * storing what is known as the local to global map.
 *
 * It is a helper class for the BSplineSpace.
 *
 * This class basically has two main (private) member:
 * - index_distribution_ that is the container for the basis function indices of a single-patch space
 * - elements_loc_to_global_flat_view_ that represent the views of the dofs that are active on each element of the space.
 *
 *
 *
 * @author pauletti, 2014
 * @author M.Martinelli, 2014
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
    using DofsComponentContainer = vector<Index>;

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


    /** @name Constructors */
    ///@{
    /** Default constructor. Not allowed to be used. */
    DofDistribution() = delete;

    //TODO: document this constructor
    DofDistribution(std::shared_ptr<CartesianGrid<dim> > grid,
                    const MultiplicityTable &accum_mult,
                    const SpaceDimensionTable &n_basis,
                    const DegreeTable &degree_table,
                    DistributionPolicy pol = DistributionPolicy::standard);

    /** Copy constructor.*/
    DofDistribution(const DofDistribution &dof_ditribution) = delete;

    /** Move constructor.*/
    DofDistribution(DofDistribution &&dof_ditribution) = default;

    /** Destructor. */
    ~DofDistribution() = default;
    //@}


    /** Assignment operators */
    ///@{
    /** Copy assignment operator. Not allowed to be used. */
    DofDistribution &operator=(const DofDistribution &dof_distribution) = delete;

    /** Move assignment operator.*/
    DofDistribution &operator=(DofDistribution &&dof_distribution) = default;
    ///@}



    /**
     * This function looks for a @p dof_id and (if found) gives back its component id @p comp
     * and its TensorIndex<dim> @p tensor_index within the component.
     *
     * @returns TRUE if the @p dof_id is found in the DofDistribution.
     *
     * @warning If the @p dof_id is NOT found in the DofDistribution the values @p comp_id and
     * @p tensor_id are UNDETERMINED.
     */
    bool find_dof_id(const Index dof_id, int &comp_id, TensorIndex<dim> &tensor_index) const;

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
     * Returns the number of active dofs of the @p element.
     */
    Size get_num_dofs_element(const CartesianGridElement<dim> &element) const;


    /** Returns the active dofs of the @p element.*/
    vector<Index> get_loc_to_global_indices(const CartesianGridElement<dim> &element) const;
    ///@}


    /**
     * Returns the container used to store the dofs ids of each component of a single patch space.
     *
     * @warning This object can have a BIG memory footprint, therefore its copy is discouraged: please
     * use the associated View instead!
     */
    const IndexDistributionTable &get_index_table() const;

    /**
     * Returns a view of the active dofs ids on a given single-patch space (non-const version).
     */
    DofsView &get_dofs_view();


    /**
     * Returns a view of the active dofs ids on a given single-patch space (const version).
     */
    const DofsView &get_dofs_view() const;



    /**
     * Returns a pointer to a std::map in which the key is the element flat_id and the value
     * is a const view to the global dofs on a single element of a space.
     * The size of the map is equal to the number of active elements in the space.
     */
    std::shared_ptr<const std::map<Index,DofsConstView>> get_elements_view() const;


    /**
     * Converts a @p global_dof_id into the correspondent local (patch) representation.
     *
     * @note The @p global_dof_id must be in the DofDistribution, otherwise in DEBUG mode
     * an assertion will be raised.
     */
    Index global_to_patch_local(const Index global_dof_id) const;

private:

    /**
     * Container used to store the dofs ids of each component of a single patch space.
     *
     * @warning This object can have a BIG memory footprint, therefore its copy is discouraged: please
     * use the associated View instead!
     */
    IndexDistributionTable index_table_;

    /**
     * View of the active dofs ids on a given single-patch space.
     */
    DofsView dofs_view_;


    /**
     * This functions uses the indices stored in the index_distribution_ member variable and
     * creates the views relative to the elements in the space.
     */
    void create_element_loc_to_global_view(
        std::shared_ptr<const CartesianGrid<dim> > grid,
        const MultiplicityTable &accum_mult,
        const SpaceDimensionTable &n_elem_basis);

    /**
     * Pointer to a std::map in which the key is the element flat_id and the value
     * is a const view to the global dofs on a single element of a space.
     * The size of the map is equal to the number of active elements in the space.
     *
     * @note We use the pointer because this object can be used by other classes (@see SpaceManager),
     * and we want to keep the synchronization of the element views without the expense of successive copies.
     */
    std::shared_ptr<std::map<Index,DofsConstView>> elements_loc_to_global_flat_view_;

    DistributionPolicy policy_;
};

IGA_NAMESPACE_CLOSE

#endif
