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

#ifndef DOFS_MANAGER_H_
#define DOFS_MANAGER_H_

#include <igatools/base/config.h>
#include <igatools/base/logstream.h>
#include <igatools/base/linear_constraint.h>
#include <igatools/utils/concatenated_iterator.h>
#include <igatools/basis_functions/equality_constraint.h>
#include <igatools/linear_algebra/sparsity_pattern.h>
//#include <boost/graph/adjacency_list.hpp>

#include <memory>

IGA_NAMESPACE_OPEN




/**
 * @brief The purpose of this class is to provide an unified way to access the dofs information provided
 * by different type of function spaces.
 *
 * For example, let consider this
 * @code{.cpp}
   auto space = BSplineSpace<3,3,1>:create();
   @endcode
 * The space represented by the object <tt>space</tt> has 3 scalar components and its dofs are
 * organized as 3 object of type DynamicMultiArray<Index,3>.
 *
 * @todo: complete the documentation
 * @author M. Martinelli
 * @date 10 Jun 2014
 */
class DofsManager
{
public:
    /** Type alias for the dofs container used in each scalar component of a single-patch space. */
    using DofsComponentContainer = std::vector<Index>;

    /** Type alias for the View on the dofs in each scalar component of a single-patch space. */
    using DofsComponentView = ContainerView<DofsComponentContainer>;

    /** Type alias for the ConstView on the dofs in each scalar component of a single-patch space. */
    using DofsComponentConstView = ConstContainerView<DofsComponentContainer>;

    /** Type alias for a concatenated iterator defined on several compoenent views. */
    using DofsIterator = ConcatenatedIterator<DofsComponentView>;

    /** Type alias for a concatenated const-iterator defined on several compoenent views. */
    using DofsConstIterator = ConcatenatedConstIterator<DofsComponentConstView>;

    /** Type alias for the View on the dofs held by each space in the DofsManager object. */
    using SpaceDofsView = View<DofsIterator,DofsConstIterator>;

    /** Type alias for the View on the dofs held by the DofsManager object. */
    using DofsView = View<DofsIterator,DofsConstIterator>;

    using DofsConstView = ConstView<DofsIterator,DofsConstIterator>;


    /** Default constructor. */
    DofsManager();

    /**
     * Prints internal information about the DofsManager.
     * @note Mostly used for debugging and testing.
     */
    void print_info(LogStream &out) const;


    /** @name Functions for changing the internal state of the DofsManager */
    ///@{
    /**
     * Sets the DofsManager in a state that can receive the views of the dofs by some spaces.
     */
    void dofs_arrangement_open();

    /**
     * Sets the DofsManager in a state that cannot receive anymore the views of the dofs by some spaces.
     *
     * If the input argument @p automatic_dofs_renumbering is set to TRUE (the default value)
     * then the dofs in each space are renumbered by the DofsManager.
     * The renumbering is made in ascending order processing the dofs space views as inserted
     * using the function add_dofs_space_view.
     *
     * If the input argument @p automatic_dofs_renumbering is set to FALSE, no renumbering is performed.
     */
    void dofs_arrangement_close(const bool automatic_dofs_renumbering = true);


    /**
     * Sets the DofsManager in a state that can receive new equality constraints.
     */
    void equality_constraints_open();

    /**
     * Coomunicate the DofsManager that the insertion of the equality constraints is
     * completed.
     */
    void equality_constraints_close();

    /**
     * Sets the DofsManager in a state that can receive new linear constraints.
     */
    void linear_constraints_open();

    /**
     * Coomunicate the DofsManager that the insertion of the linear constraints is
     * completed.
     */
    void linear_constraints_close();

    /**
     * Add an equality constraint between @dof_id_master and @p dof_id_slave.
     */
    void add_equality_constraint(const Index dof_id_master,const Index dof_id_slave);

    /**
     * Adds the view to the dofs of a space to the vector of views held by the DofsManager.
     * @param[in] space_id The identifier of the space.
     * @param[in] num_dofs_space Number of dofs that are represented by the view that is added by this function.
     * @param[in] dofs_space_view View of the dofs that must be added to the DofsManager.
     * @pre In order to call this function, the DofsManager must be be in the state that permits to receive
     * the dofs view. In other words, the user should call dofs_arrangement_open().
     */
    void add_dofs_space_view(const int space_id,
                             const Index num_dofs_space,
                             const SpaceDofsView &dofs_space_view);
    ///@}

    /** @name Functions for querying dofs information */
    ///@{
    DofsView &get_dofs_view();

    DofsConstView get_dofs_view() const;


    /**
     * Returns the global dof corresponding to the @p local_dof
     * in the space with id equal to @p space_id.
     */
    Index get_global_dof(const int space_id, const Index local_dof) const;

    /**
     * Returns the global dofs corresponding to the @p local_dofs
     * in the space with id equal to @p space_id.
     */
    std::vector<Index> get_global_dofs(const int space_id, const std::vector<Index> &local_dof) const;
    ///@}


    /** Return the number of unique dofs in the MultiPatchSpace. */
    Index get_num_dofs() const;


    /** Returns the number of linear constraints. */
    Index get_num_linear_constraints() const;


    /** Returns the number of equality constraints. */
    Index get_num_equality_constraints() const;



    /**
     * Removes the equality constraints redundancies.
     */
    void remove_equality_constraints_redundancies();


    /** Returns the sparsity pattern of the dofs handled by the DofsManager. */
    SparsityPattern get_sparsity_pattern() const;

private:
    bool is_dofs_arrangement_open_ = false;

    bool are_equality_constraints_open_ = false;
    bool are_linear_constraints_open_ = false;

    std::vector<DofsComponentView> dofs_components_view_;

    struct SpaceInfo
    {
        SpaceInfo() = delete;
        SpaceInfo(const Index n_dofs, const SpaceDofsView &dofs_view);

        Index n_dofs_;
        Index offset_;
        SpaceDofsView dofs_view_;
    };

    std::map<int,SpaceInfo> spaces_info_;



    std::unique_ptr<DofsView> dofs_view_;


    std::vector<DofsConstView> elements_dofs_view_;



    std::vector<LinearConstraint> linear_constraints_;




    std::vector<EqualityConstraint> equality_constraints_;


    /** Counts and return the number of unique dofs in the MultiPatchSpace. */
    Index count_unique_dofs() const;

    /** Number of unique dofs in the MultiPatchSpace. */
    Index num_unique_dofs_;
};


IGA_NAMESPACE_CLOSE


#endif // #ifndef DOFS_MANAGER_H_
