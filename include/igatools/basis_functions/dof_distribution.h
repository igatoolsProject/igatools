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

#ifndef DOF_DISTRIBUTION_H_
#define DOF_DISTRIBUTION_H_

#include <igatools/base/config.h>
#include <igatools/basis_functions/spline_space.h>
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
 * This class basically has two one (private) member:
 * - index_distribution_ that is the container for the basis function indices
 *   of a single-patch space
 *
 * @author pauletti, 2014
 * @author M.Martinelli, 2014, 2015
 *
 */
template<int dim, int range = 1, int rank = 1>
class DofDistribution
{
public:
    using Space = SplineSpace<dim, range, rank>;
    using DegreeTable = typename Space::DegreeTable;
    using PeriodicTable = typename Space::PeriodicTable;
    using TensorSizeTable = typename Space::TensorSizeTable;
    using OffsetTable = typename Space::template ComponentContainer<Size>;
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
    DofDistribution(const TensorSizeTable &n_basis,
                    const DegreeTable &degree_table,
                    const PeriodicTable &periodic,
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
     * This function looks for a @p dof_id and (if found) gives back its
     * component id @p comp
     * and its TensorIndex<dim> @p tensor_index within the component.
     *
     * @returns TRUE if the @p dof_id is found in the DofDistribution.
     *
     * @warning If the @p dof_id is NOT found in the DofDistribution the values
     *  @p comp_id and
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



    /**
     * Returns the container used to store the dofs ids of each component
     * of a single patch space.
     *
     * @warning This object can have a BIG memory footprint, therefore
     * its copy is discouraged: please use the associated View instead!
     */
    const IndexDistributionTable &get_index_table() const;

    /**
     * Component table with the offset of unique dofs
     * in each component.
     */
    OffsetTable get_dofs_offset() const;


    /**
     * Returns a view of the active dofs ids on a given single-patch space (non-const version).
     */
    DofsView &get_dofs_view();


    /**
     * Returns a view of the active dofs ids on a given single-patch space (const version).
     */
    const DofsView &get_dofs_view() const;


    /**
     * Return the container holding the number of unique dofs, component-by-component and direction-by-direction.
     */
    const TensorSizeTable &get_num_dofs_table() const;

    /**
     * Converts a @p global_dof_id into the correspondent local (patch) representation.
     *
     * @note The @p global_dof_id must be in the DofDistribution, otherwise in DEBUG mode
     * an assertion will be raised.
     */
    Index global_to_patch_local(const Index global_dof_id) const;


    /**
     * @name Functions related to the management/query of the dof properties.
     */
    ///@{

    /**
     * Returns true if the dof with id @p dof_id has the asked @p property.
     */
    bool test_if_dof_has_property(const Index dof_id, const std::string &property) const;


    /**
     * Adds a new <tt>property</tt> definition for the dofs in the DofDistribution.
     *
     * @note If the <tt>property</tt> is already present, n assertion will be raised (in Debug mode).
     */
    void add_dofs_property(const std::string &property);


    /**
     * Returns the id of the dofs having a certain @p property (non-const version).
     */
    std::set<Index> &get_dofs_id_same_property(const std::string &property);

    /**
     * Returns the id of the dofs having a certain @p property (const version).
     */
    const std::set<Index> &get_dofs_id_same_property(const std::string &property) const;


    /**
     * Sets the @p status of the given @p property for the dof with the specified @p dof_id.
     */
    void set_dof_property_status(const std::string &property, const Index dof_id, const bool status);

    /**
     * Sets the @p status of the given @p property for all the dofs in the DofDistribution.
     */
    void set_all_dofs_property_status(const std::string &property, const bool status);
    ///@}

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
     * Number of unique dofs, component-by-component and direction-by-direction.
     */
    TensorSizeTable num_dofs_table_;



    /**
     * Size of the index_table_, component-by-component and direction-by-direction.
     */
    TensorSizeTable index_table_size_;

    DistributionPolicy policy_;


    /**
     * Container for the dofs having a certain property.
     *
     * The property name is the key of the std::map.
     */
    PropertiesIdContainer properties_dofs_;
};

IGA_NAMESPACE_CLOSE

#endif
