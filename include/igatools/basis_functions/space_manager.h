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

#ifndef SPACE_MANAGER_H_
#define SPACE_MANAGER_H_

#include <igatools/base/config.h>
#include <igatools/base/logstream.h>
#include <igatools/base/equality_constraint.h>
#include <igatools/base/linear_constraint.h>
#include <igatools/utils/concatenated_iterator.h>
#include <igatools/basis_functions/nurbs_space.h>
#include <igatools/basis_functions/physical_space.h>


#include <igatools/contrib/variant.h>
#include <boost/variant.hpp>

#include <memory>
#include <set>

IGA_NAMESPACE_OPEN


template<class RefSpace,class PushFwd>
using PhysSpacePtr = std::shared_ptr<PhysicalSpace<RefSpace,PushFwd>>;

static const int rank = 1;

using SpacePtrVariant =
    Variant<
    std::shared_ptr<BSplineSpace<0,1,rank>>,
    std::shared_ptr<BSplineSpace<1,1,rank>>,
    std::shared_ptr<BSplineSpace<2,1,rank>>,
    std::shared_ptr<BSplineSpace<3,1,rank>>,
    std::shared_ptr<BSplineSpace<0,2,rank>>,
    std::shared_ptr<BSplineSpace<1,2,rank>>,
    std::shared_ptr<BSplineSpace<2,2,rank>>,
    std::shared_ptr<BSplineSpace<3,2,rank>>,
    std::shared_ptr<BSplineSpace<0,3,rank>>,
    std::shared_ptr<BSplineSpace<1,3,rank>>,
    std::shared_ptr<BSplineSpace<2,3,rank>>,
    std::shared_ptr<BSplineSpace<3,3,rank>>,
    std::shared_ptr<NURBSSpace<0,1,rank>>,
    std::shared_ptr<NURBSSpace<1,1,rank>>,
    std::shared_ptr<NURBSSpace<2,1,rank>>,
    std::shared_ptr<NURBSSpace<3,1,rank>>,
    std::shared_ptr<NURBSSpace<0,2,rank>>,
    std::shared_ptr<NURBSSpace<1,2,rank>>,
    std::shared_ptr<NURBSSpace<2,2,rank>>,
    std::shared_ptr<NURBSSpace<3,2,rank>>,
    std::shared_ptr<NURBSSpace<0,3,rank>>,
    std::shared_ptr<NURBSSpace<1,3,rank>>,
    std::shared_ptr<NURBSSpace<2,3,rank>>,
    std::shared_ptr<NURBSSpace<3,3,rank>>,
    PhysSpacePtr<BSplineSpace<0,1,rank>,PushForward<Transformation::h_grad,0,0>>,
    PhysSpacePtr<BSplineSpace<0,2,rank>,PushForward<Transformation::h_grad,0,0>>,
    PhysSpacePtr<BSplineSpace<0,3,rank>,PushForward<Transformation::h_grad,0,0>>,
    PhysSpacePtr<NURBSSpace<0,1,rank>,PushForward<Transformation::h_grad,0,0>>,
    PhysSpacePtr<NURBSSpace<0,2,rank>,PushForward<Transformation::h_grad,0,0>>,
    PhysSpacePtr<NURBSSpace<0,3,rank>,PushForward<Transformation::h_grad,0,0>>,
    PhysSpacePtr<BSplineSpace<1,1,rank>,PushForward<Transformation::h_grad,1,0>>,
    PhysSpacePtr<BSplineSpace<0,1,rank>,PushForward<Transformation::h_grad,0,1>>,
    PhysSpacePtr<BSplineSpace<1,2,rank>,PushForward<Transformation::h_grad,1,0>>,
    PhysSpacePtr<BSplineSpace<1,3,rank>,PushForward<Transformation::h_grad,1,0>>,
    PhysSpacePtr<BSplineSpace<0,2,rank>,PushForward<Transformation::h_grad,0,1>>,
    PhysSpacePtr<BSplineSpace<0,3,rank>,PushForward<Transformation::h_grad,0,1>>,
    PhysSpacePtr<NURBSSpace<1,1,rank>,PushForward<Transformation::h_grad,1,0>>,
    PhysSpacePtr<NURBSSpace<0,1,rank>,PushForward<Transformation::h_grad,0,1>>,
    PhysSpacePtr<NURBSSpace<1,2,rank>,PushForward<Transformation::h_grad,1,0>>,
    PhysSpacePtr<NURBSSpace<1,3,rank>,PushForward<Transformation::h_grad,1,0>>,
    PhysSpacePtr<NURBSSpace<0,2,rank>,PushForward<Transformation::h_grad,0,1>>,
    PhysSpacePtr<NURBSSpace<0,3,rank>,PushForward<Transformation::h_grad,0,1>>,
    PhysSpacePtr<BSplineSpace<1,1,rank>,PushForward<Transformation::h_grad,1,1>>,
    PhysSpacePtr<BSplineSpace<1,2,rank>,PushForward<Transformation::h_grad,1,1>>,
    PhysSpacePtr<BSplineSpace<2,1,rank>,PushForward<Transformation::h_grad,2,0>>,
    PhysSpacePtr<BSplineSpace<2,2,rank>,PushForward<Transformation::h_grad,2,0>>,
    PhysSpacePtr<BSplineSpace<0,1,rank>,PushForward<Transformation::h_grad,0,2>>,
    PhysSpacePtr<BSplineSpace<0,2,rank>,PushForward<Transformation::h_grad,0,2>>,
    PhysSpacePtr<BSplineSpace<2,3,rank>,PushForward<Transformation::h_grad,2,0>>,
    PhysSpacePtr<BSplineSpace<1,3,rank>,PushForward<Transformation::h_grad,1,1>>,
    PhysSpacePtr<NURBSSpace<1,1,rank>,PushForward<Transformation::h_grad,1,1>>,
    PhysSpacePtr<NURBSSpace<1,2,rank>,PushForward<Transformation::h_grad,1,1>>,
    PhysSpacePtr<NURBSSpace<2,1,rank>,PushForward<Transformation::h_grad,2,0>>,
    PhysSpacePtr<NURBSSpace<2,2,rank>,PushForward<Transformation::h_grad,2,0>>,
    PhysSpacePtr<NURBSSpace<0,1,rank>,PushForward<Transformation::h_grad,0,2>>,
    PhysSpacePtr<NURBSSpace<0,2,rank>,PushForward<Transformation::h_grad,0,2>>,
    PhysSpacePtr<NURBSSpace<2,3,rank>,PushForward<Transformation::h_grad,2,0>>,
    PhysSpacePtr<NURBSSpace<1,3,rank>,PushForward<Transformation::h_grad,1,1>>,
    PhysSpacePtr<BSplineSpace<1,1,rank>,PushForward<Transformation::h_grad,1,2>>,
    PhysSpacePtr<BSplineSpace<1,3,rank>,PushForward<Transformation::h_grad,1,2>>,
    PhysSpacePtr<BSplineSpace<2,1,rank>,PushForward<Transformation::h_grad,2,1>>,
    PhysSpacePtr<BSplineSpace<2,3,rank>,PushForward<Transformation::h_grad,2,1>>,
    PhysSpacePtr<BSplineSpace<3,1,rank>,PushForward<Transformation::h_grad,3,0>>,
    PhysSpacePtr<BSplineSpace<3,3,rank>,PushForward<Transformation::h_grad,3,0>>,
    PhysSpacePtr<BSplineSpace<0,1,rank>,PushForward<Transformation::h_grad,0,3>>,
    PhysSpacePtr<BSplineSpace<0,3,rank>,PushForward<Transformation::h_grad,0,3>>,
    PhysSpacePtr<NURBSSpace<1,1,rank>,PushForward<Transformation::h_grad,1,2>>,
    PhysSpacePtr<NURBSSpace<1,3,rank>,PushForward<Transformation::h_grad,1,2>>,
    PhysSpacePtr<NURBSSpace<2,1,rank>,PushForward<Transformation::h_grad,2,1>>,
    PhysSpacePtr<NURBSSpace<2,3,rank>,PushForward<Transformation::h_grad,2,1>>,
    PhysSpacePtr<NURBSSpace<3,1,rank>,PushForward<Transformation::h_grad,3,0>>,
    PhysSpacePtr<NURBSSpace<3,3,rank>,PushForward<Transformation::h_grad,3,0>>,
    PhysSpacePtr<NURBSSpace<0,1,rank>,PushForward<Transformation::h_grad,0,3>>,
    PhysSpacePtr<NURBSSpace<0,3,rank>,PushForward<Transformation::h_grad,0,3>>>;



/**
 * @brief The purpose of this class is to provide an unified way to access the dofs information provided
 * by different type of function spaces.
 *
 * In order to correctly use the Space manager, the user is FORCED to perform (in order) the
 * following MANDATORY steps:
 *   1. <b>Insertion of the spaces.</b>
 *      The purpose of this phase is to populate the SpaceManager with the spaces that must handle
 *      and coordinate. This phase starts with the execution function spaces_insertion_open() and it is
 *      concluded after the execution of the function spaces_insertion_close().
 *      Between the execution of these two functions, the spaces are added to the SpaceManager with
 *      the function add_space().
 *      The function spaces_insertion_close() has an (optional) input argument that is used to
 *      choose if the global dofs of the spaces must keep the original numbering or must be renumbered
 *      in a way that no (global) dof id is used in different spaces (this is the default behaviour).
 *      The "no-renumbering" option can be useful if the global dof numbering needs to be kept
 *      (for any reason).
 *      @code{.cpp}
        SpaceManager sp_manager;

        // Phase 1
        sp_manager.spaces_insertion_open(); // starts the spaces insertion phase (phase 1)

        sp_manager.add_space(space_0); // adds space_0
        sp_manager.add_space(space_1); // adds space_1
        sp_manager.add_space(space_2); // adds space_2
        ...
        sp_manager.add_space(space_last);  // adds space_last

        sp_manager.spaces_insertion_close(); // ends the spaces insertion phase (phase 1) and renumbers the global dof id

        @endcode
 *
 *      <em>After this phase is not possible to insert any new space to the SpaceManager.</em>
 *
 *   2. <b>Definition of the connectivity between the spaces.</b>
 *      The purpose of this phase is to define relation between the dofs in the spaces and therefore
 *      the active blocks in the <em>sparsity pattern</em> of the system matrix.
 *      This phase starts with the execution function spaces_connectivity_open() and it is
 *      concluded after the execution of the function spaces_connectivity_close().
 *      Between the execution of these two functions, the spaces connectivities
 *      (i.e. the blocks in the sparsity pattern) are added to the SpaceManager with the function
 *      add_spaces_connection(). Please note that this last function takes two arguments: one argument
 *      the for the space of test functions and the other for the space of trial functions.
 *      This means that if a user wants to add symmetric blocks, he needs to call this function twice with
 *      the spaces reverted. The blocks on the diagonal are added calling add_spaces_connection()
 *      with the same spaces on both argument.
 *
 *      For example, if we have 3 spaces <tt>space_0</tt>, <tt>space_1</tt> and <tt>space_2</tt> and the
 *      connectivity is defined as follows:
 *        - <tt>space_0</tt> with <tt>space_2</tt> and it symmetric counterpart;
 *        - <tt>space_1</tt> with <tt>space_2</tt>;
 *        - <tt>space_0</tt> with itself (diagonal block);
 *        - <tt>space_1</tt> with itself (diagonal block)
 *
 *      we should use:
 *      @code{.cpp}
        // Phase 2
        sp_manager.spaces_connectivity_open();  // starts the spaces connectivity phase (phase 2)

        sp_manager.add_spaces_connection(space_0,space_2); // block 0-2
        sp_manager.add_spaces_connection(space_2,space_0); // block 2-0

        sp_manager.add_spaces_connection(space_1,space_2); // block 1-2

        sp_manager.add_spaces_connection(space_0,space_0); // block 0-0
        sp_manager.add_spaces_connection(space_1,space_1); // block 1-1

        sp_manager.spaces_connectivity_close();  // ends the spaces connectivity phase (phase 2)

        @endcode
 *
 *   3. <b>Definition of the <em>global dofs connectivity</em> in each spaces pair</b>.
 *   Once the relations between the spaces are defined (Phase 2), it remains to define the dofs
 *   relations within a spaces pair.
 *     - For the case in which the spaces connection is not defined with a single space
 *   (i.e. we are not in the block-diagonal), the user has the responsability to manage
 *   (and set) the dofs connectivity. The function for setting the dofs connectivity with a spaces pair
 *   is spaces_connection_set_dofs_connectivity(). The first two arguments of this function are the
 *   spaces defining the block for which we want to set the dofs connectivity, the third arguments is a
 *   <tt>std::map<Index,std::set<Index>> </tt>
 *   in which a generic key represents an active global dofs in the block (belonging from the first space)
 *   and the value is a <tt>std::set<Index></tt> representing the global dofs in the second space that
 *   are in relation with the global dof in the first space represented by the key.
 *   @code{.cpp}
     sp_manager.spaces_connection_set_dofs_connectivity(space_0,space_2,dofs_connectivity_0_2);
     @endcode
 *     - In the particular case in which a space defines a connection
 *   with itself (i.e. a the connection is represented by a diagonal block),
 *   the dofs relations are automatically inferred (without taking any action) from the
 *   internal structure of the space, i.e. iterating over the active elements of the space and retrieving
 *   the dofs connectivity in each element.
 *   If this default behaviour is not appropriate to describe the dofs connectivity,
 *   the user can set up manually the dofs connectivity by using the function
 *   spaces_connection_set_dofs_connectivity().
 *
 * @todo: complete the documentation
 * @author M. Martinelli
 * @date 10 Jun 2014
 */
class SpaceManager
{
public:
    /** Type alias for the dofs container used in each scalar component of a single-patch space. */
    using DofsComponentContainer = vector<Index>;

    /** Type alias for the View on the dofs in each scalar component of a single-patch space. */
    using DofsComponentView = ContainerView<DofsComponentContainer>;

    /** Type alias for the ConstView on the dofs in each scalar component of a single-patch space. */
    using DofsComponentConstView = ConstContainerView<DofsComponentContainer>;

    /** Type alias for a concatenated iterator defined on several compoenent views. */
    using DofsIterator = ConcatenatedIterator<DofsComponentView>;

    /** Type alias for a concatenated const-iterator defined on several compoenent views. */
    using DofsConstIterator = ConcatenatedConstIterator<DofsComponentView,DofsComponentConstView>;

    /** Type alias for the View on the dofs held by the SpaceManager object. */
    using DofsView = View<DofsIterator,DofsConstIterator>;

    /** Type alias for the ConstView on the dofs held by the SpaceManager object. */
    using DofsConstView = ConstView<DofsIterator,DofsConstIterator>;

    /** Type alias for the LinearConstraint. */
    using LC = LinearConstraint;

    /** Type alias for the LinearConstraintType. */
    using LCType = LinearConstraintType;


    /** @name Constructors */
    ///@{
    /** Default constructor. */
    SpaceManager();

    /** Copy constructor. */
    SpaceManager(const SpaceManager &space_manager) = delete;

    /** Move constructor. */
    SpaceManager(SpaceManager &&space_manager) = default;
    ///@}

    /** Destructor. */
    ~SpaceManager() = default;

    /**
     * Prints internal information about the SpaceManager.
     * @note Mostly used for debugging and testing.
     */
    void print_info(LogStream &out) const;


    /** @name Functions for managing the spaces in the SpaceManager */
    ///@{
    /**
     * Sets the SpaceManager in a state that can receive the insertion of new spaces.
     */
    void spaces_insertion_open();


    /**
     * Adds a space to the SpaceManager.
     *
     * @note An assertion will be raised (in DEBUG mode)
     * if the passed <p>space</p> is already present in the SpaceManager.
     */
    template<class Space>
    void add_space(std::shared_ptr<Space> space);


    /**
     * Sets the SpaceManager in a state that cannot receive any new space.
     *
     * If the input argument @p automatic_dofs_renumbering is set to TRUE (the default value)
     * then the dofs in each space are renumbered by the SpaceManager.
     * The renumbering is made in ascending order processing the dofs space views as inserted
     * using the function add_dofs_space_view.
     *
     * If the input argument @p automatic_dofs_renumbering is set to FALSE, no renumbering is performed.
     */
    void spaces_insertion_close(const bool automatic_dofs_renumbering = true);


    /**
     * Returns true if the spaces insertion phase is open.
     */
    bool is_spaces_insertion_open() const;
    ///@}


    /** @name Functions for managing the spaces connecitivty in the SpaceManager */
    ///@{
    /**
     * Sets the SpaceManager in a state that can receive the connectivity between existing spaces.
     */
    void spaces_connectivity_open();

    /**
     * Defines connection between the global dofs of the @p space_test (i.e. the row dofs id)
     * and the global dofs of the @p space_trial (i.e. the column dofs id)
     */
    template<class SpaceTest,class SpaceTrial>
    void add_spaces_connection(std::shared_ptr<SpaceTest> space_test,std::shared_ptr<SpaceTrial> space_trial);

    /**
     * Sets the SpaceManager in a state that cannnot receive any new connectivity between existing spaces.
     */
    void spaces_connectivity_close();


    /**
     * Returns true if the spaces connectivity phase is open.
     */
    bool is_spaces_connectivity_open() const;
    ///@}




    /** @name Functions for the insertion of the equality constraints */
    ///@{
    /**
     * Sets the SpaceManager in a state that can receive new equality constraints.
     */
    void equality_constraints_open();

    /**
     * Add an equality constraint between @dof_id_master and @p dof_id_slave.
     *
     * @note An assertion will be raised (in DEBUG mode)
     * if the space manager is not set in the proper state by the function
     * equality_constraints_open().
     */
    void add_equality_constraint(const Index dof_id_master,const Index dof_id_slave);

    /**
     * Communicate the SpaceManager that the insertion of the equality constraints is
     * completed.
     */
    void equality_constraints_close();
    ///@}

    /** @name Functions for the management of the linear constraints */
    ///@{
    /**
     * Sets the SpaceManager in a state that can receive new linear constraints.
     */
    void linear_constraints_open();


    /**
     * Add a LinearConstraint of a given @p type to the SpaceManager,
     * where @p global_dof_id is the global dof associated with the linear constraint,
     * @p dofs are the global dofs id involved by the constraint,
     * @p coeffs their coefficients and
     * @p rhs is the right hand side that defines the linear constraint equation.
     *
     * @note An assertion will be raised (in DEBUG mode)
     * if the space manager is not set in the proper state by the function
     * linear_constraints_open().
     */
    void add_linear_constraint(const Index global_dof_id,
                               const LinearConstraintType &type,
                               const vector<Index> &dofs,
                               const vector<Real> &coeffs,
                               const Real rhs);


    /**
     * Add a LinearConstraint to the SpaceManager.
     *
     * @note An assertion will be raised (in DEBUG mode)
     * if the space manager is not set in the proper state by the function
     * linear_constraints_open().
     */
    void add_linear_constraint(std::shared_ptr<LinearConstraint> linear_constraint);

    /**
     * Communicate the SpaceManager that the insertion of the linear constraints is
     * completed.
     */
    void linear_constraints_close();

    /**
     * Returns the vector of linear constraints of a given @p type, defined in the SpaceManager.
     * If no @type is passed at the input argument, the function returns all the linear constraints.
     */
    vector<std::shared_ptr<LinearConstraint> > get_linear_constraints(
        const LinearConstraintType &type = LinearConstraintType::any) const;


    /**
     * This function tests if a vector of global dof values @p dof_values, satisfies the linear constraints
     * (up to the tolerance @p tol).
     * If all coefficients satisfies the linear constraints the function returns an empty vector,
     * otherwise it returns a vector containing the linear constraints that are not satisfied.
     */
    vector<std::shared_ptr<LinearConstraint> > verify_linear_constraints(
        const vector<Real> &dof_values,
        const Real tol = 1.0e-13) const;
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
    vector<Index> get_global_dofs(const int space_id, const vector<Index> &local_dof) const;
    ///@}


    /** Return the number of unique dofs in the SpaceManager. */
    Index get_num_dofs() const;


    /** Returns the number of linear constraints. */
    Index get_num_linear_constraints() const;


    /** Returns the number of equality constraints. */
    Index get_num_equality_constraints() const;



    /**
     * Removes the equality constraints redundancies.
     */
    void remove_equality_constraints_redundancies();



    void add_element_dofs_view(const DofsConstView &element_dofs_view);


private:
    bool is_spaces_insertion_open_ = false;
    bool is_spaces_connectivity_open_ = false;
    bool are_equality_constraints_open_ = false;
    bool are_linear_constraints_open_ = false;


    /**
     * This is an helper class in order to store a pointer to a certain space variant,
     * and store some useful space information, without the necessity to know the template parameters
     * defining the space.
     *
     * This is useful for instance, if the type of space is not known at compile time
     * (e.g. loading the space from a file) and/or use the space in a class that is
     * non templatized w.r.t. the space type.
     *
     * All data in this class (except the pointer to the space itself, that is managed by the class Variant),
     * does not depends on the template parameters defining the space.
     */
    class SpaceInfo
    {
    public:
        SpaceInfo() = delete;
        SpaceInfo(const SpacePtrVariant &space,
                  const Index id,
                  const int dim,
                  const int codim,
                  const int space_dim,
                  const int range,
                  const int rank,
                  const Index num_dofs,
                  const Index min_dofs_id,
                  const Index max_dofs_id,
                  const DofsView &dofs_view,
                  const std::shared_ptr<const std::map<Index,DofsConstView>> elements_dofs_view);

        SpaceInfo(const SpaceInfo &sp) = delete;
        SpaceInfo(SpaceInfo &&sp) = delete;

        SpaceInfo &operator=(const SpaceInfo &sp) = delete;
        SpaceInfo &operator=(SpaceInfo &&sp) = delete;

        /**
         * Returns the space global id.
         */
        Index get_id() const;

        /** Returns the number of dofs of the space. */
        Index get_num_dofs() const ;


        /*
         * Returns a View to the dofs active on the space (non-const version).
         */
        DofsView &get_dofs_view();

        /*
         * Returns a View to the dofs active on the space (const version).
         */
        const DofsView &get_dofs_view() const;

        /**
         * Returns a vector of size equal to the number of elements in the single-patch space,
         * for which each entry is a view of the global dofs ids active on the element.
         */
        const std::map<Index,DofsConstView> &get_elements_dofs_view() const;




        /** Returns the minimum dof id present in the space.*/
        Index get_min_dofs_id() const;

        /** Returns the maximum dof id present in the space.*/
        Index get_max_dofs_id() const;

        /**
         * Addas an @p offset to all dofs present in the space.
         */
        void add_dofs_offset(const Index offset);

        /**
         * Return an object containing a variant of a shared_pointer pointing to a certain single-patch space.
         * The allowed space type can be any valid BSplineSpace, NURBSSpace or PhysicalSpace.
         */
        SpacePtrVariant &get_space_variant();


        /**
         * Return true if the id of the lhs is equal to the id of the rhs.
         */
        bool operator==(const SpaceInfo &sp) const;

    private:
        /**
         * Pointer to a generic single-patch space (it can be any of the type allowed for BSplineSpace,
         * NURBSSpace and PhysicalSpace).
         */
        SpacePtrVariant space_;

        /**
         * Global space id.
         */
        const Index id_;


        int dim_;

        int codim_;

        int space_dim_;

        int range_;

        int rank_;

        /**
         * Nuber of active dofs in the space.
         */
        Index num_dofs_;

        /** Minumum dof id in the space. */
        Index min_dofs_id_;

        /** Maximum dof id in the space. */
        Index max_dofs_id_;

        /**
         * View of the dofs ids active on the space.
         */
        DofsView dofs_view_;

        /**
         * Map of size equal to the number of elements in the single-patch space,
         * for which each entry is a view of the global dofs ids active on the element.
         *
         * The std::map key represent the element flat-id for which we store the dofs view.
         *
         * @note We use a std:shared_ptr because this container can be very big and
         * it is already present
         * the the DofDistribution class instantiated in the space itself.
         */
        std::shared_ptr<const std::map<Index,DofsConstView>> elements_dofs_view_;
    };

    using SpaceInfoPtr = std::shared_ptr<SpaceInfo>;

    /**
     * This functor defines the relation order between the spaces.
     */
    struct space_comparator_less
    {
        bool operator()(const SpaceInfoPtr &lhs,const SpaceInfoPtr &rhs) const
        {
            return lhs->get_id() < rhs->get_id();
        }
    };



    /**
     * Map containing the pointers to the spaces handled by the SpaceManager,
     * and some useful informations that does not depends on the template
     * parameters needed to instantiate the spaces.
     *
     * The std::map key is the space id.
     */
    std::map<int,SpaceInfoPtr> spaces_info_;

#if 0
    /**
     * Connectivity between the different spaces/blocks.
     * The key in the std::map represents a block of row dofs
     * and the values of a give key represents the connected
     * blocks of dofs along the columns of the system matrix.
     */
    std::map<SpaceInfoPtr,std::set<SpaceInfoPtr>> spaces_connectivity_;
#endif

    /**
     * View to the active dofs ids of all single-patch spaces handled by the SpaceManager.
     */
    DofsView dofs_view_;


    std::multimap<LCType,std::shared_ptr<LC>> linear_constraints_;


    vector<EqualityConstraint> equality_constraints_;


    /** Counts and return the number of unique dofs in the SpaceManager. */
    Index count_unique_dofs() const;


    /** Number of unique dofs in the SpaceManager. */
    Index num_unique_dofs_;


    class SpacesConnection
    {
    public:
        SpacesConnection() = delete;

        SpacesConnection(const SpaceInfoPtr &space_row,const SpaceInfoPtr &space_col);

        ~SpacesConnection() = default;

        SpacesConnection(const SpacesConnection &in) = default;
        SpacesConnection(SpacesConnection &&in) = default;

        SpacesConnection &operator=(const SpacesConnection &in) = delete;
        SpacesConnection &operator=(SpacesConnection &&in) = delete;

        bool operator==(const SpacesConnection &conn) const;

        /** Returns true if the row space and the column space are equal. */
        bool is_unique_space() const;

        using DofsConnectivity = std::map<Index,std::set<Index>>;

        void add_dofs_connectivity(const DofsConnectivity &dofs_connectivity)
        {
            Assert(false,ExcNotImplemented());
            AssertThrow(false,ExcNotImplemented());
        }

        const SpaceInfo &get_space_row() const
        {
            return *space_row_;
        }

        const SpaceInfo &get_space_col() const
        {
            return *space_col_;
        }

        std::set<Index> get_row_dofs() const
        {
            const auto &sp_dofs_view = space_row_->get_dofs_view();

            std::set<Index> dofs(sp_dofs_view.begin(),sp_dofs_view.end());

            for (const auto &pair : extra_dofs_connectivity_)
                dofs.insert(pair.first);

            return dofs;
        }

        std::set<Index> get_col_dofs() const
        {
            const auto &sp_dofs_view = space_col_->get_dofs_view();

            std::set<Index> dofs(sp_dofs_view.begin(),sp_dofs_view.end());

            for (const auto &pair : extra_dofs_connectivity_)
                dofs.insert(pair.second.begin(),pair.second.end());

            return dofs;
        }

    private:
        SpaceInfoPtr space_row_;
        SpaceInfoPtr space_col_;


        DofsConnectivity extra_dofs_connectivity_;

    };


    vector<SpacesConnection> spaces_connections_;

public:

    const vector<SpacesConnection> &get_spaces_connections() const
    {
        return spaces_connections_;
    }


    /**
     * Returns a map containing the pointers to the spaces handled by the SpaceManager,
     * and some useful informations that does not depends on the template
     * parameters needed to instantiate the spaces.
     */
    const std::map<Index,SpaceInfoPtr> &get_spaces_info() const;



    template<class SpaceTest,class SpaceTrial>
    SpacesConnection &get_spaces_connection(
        std::shared_ptr<SpaceTest> space_test,
        std::shared_ptr<SpaceTrial> space_trial);


    std::set<Index> get_row_dofs() const
    {
        std::set<Index> row_dofs;

        for (const auto &sp_conn : spaces_connections_)
        {
            const auto row_dofs_current_space = sp_conn.get_row_dofs();
            row_dofs.insert(row_dofs_current_space.begin(),row_dofs_current_space.end());
        }
        return row_dofs;
    }

    std::set<Index> get_col_dofs() const
    {
        std::set<Index> col_dofs;

        for (const auto &sp_conn : spaces_connections_)
        {
            const auto col_dofs_current_space = sp_conn.get_col_dofs();
            col_dofs.insert(col_dofs_current_space.begin(),col_dofs_current_space.end());
        }
        return col_dofs;
    }

};



template<class Space>
inline
void
SpaceManager::
add_space(std::shared_ptr<Space> space)
{
    Assert(space != nullptr,ExcNullPtr());

#ifndef NDEBUG
    // check that the input space is not already added
    for (auto &space_info : spaces_info_)
    {
        Assert(space != get<std::shared_ptr<Space>>(space_info.second->get_space_variant()),
               ExcMessage("Space already added in the SpaceManager."));
    }
#endif

    //------------------------------------------------------------------------
    auto &dof_distribution = space->get_dof_distribution_global();

    spaces_info_[space->get_id()] = std::shared_ptr<SpaceInfo>(
                                        new SpaceInfo(space,
                                                      space->get_id(),
                                                      Space::dim,
                                                      Space::codim,
                                                      Space::space_dim,
                                                      Space::range,
                                                      Space::rank,
                                                      space->get_num_basis(),
                                                      dof_distribution.get_min_dof_id(),
                                                      dof_distribution.get_max_dof_id(),
                                                      dof_distribution.get_dofs_view(),
                                                      dof_distribution.get_elements_view()));
    //---------------------------------------------------------------------------------------------
}

template<class SpaceTest,class SpaceTrial>
inline
void
SpaceManager::
add_spaces_connection(std::shared_ptr<SpaceTest> space_test,std::shared_ptr<SpaceTrial> space_trial)
{
    using std::cout;
    using std::endl;
    cout << "SpaceManager::add_spaces_connection()" << endl;
    cout << "adding connection between space " << space_test->get_id()
         << " and space " << space_trial->get_id() << endl;

    Assert(is_spaces_insertion_open_ == false,ExcInvalidState());
    Assert(is_spaces_connectivity_open_ == true,ExcInvalidState());

    Assert(space_test !=nullptr,ExcNullPtr());
    Assert(space_trial !=nullptr,ExcNullPtr());

    auto sp_test  = spaces_info_.at(space_test ->get_id());
    auto sp_trial = spaces_info_.at(space_trial->get_id());

//    spaces_connectivity_[sp_test].emplace(sp_trial);

    SpacesConnection conn(sp_test,sp_trial);

    Assert(std::count(spaces_connections_.begin(),spaces_connections_.end(),conn) == 0,
           ExcMessage("Spaces connection already added."));
    spaces_connections_.push_back(conn);
}

template<class SpaceTest,class SpaceTrial>
inline
auto
SpaceManager::
get_spaces_connection(
    std::shared_ptr<SpaceTest> space_test,
    std::shared_ptr<SpaceTrial> space_trial)
-> SpacesConnection &
{
    Assert(is_spaces_insertion_open_ == false,ExcInvalidState());
//    Assert(is_spaces_connectivity_open_ == true,ExcInvalidState());

    Assert(space_test !=nullptr,ExcNullPtr());
    Assert(space_trial !=nullptr,ExcNullPtr());

    auto sp_test  = spaces_info_.at(space_test ->get_id());
    auto sp_trial = spaces_info_.at(space_trial->get_id());

    auto it = std::find(
        spaces_connections_.begin(),
        spaces_connections_.end(),
        SpacesConnection(space_test,sp_trial));

    Assert(it != spaces_connections_.end(),
    ExcMessage("Spaces connection not found."));
    return *it;
}



IGA_NAMESPACE_CLOSE


#endif // #ifndef SPACE_MANAGER_H_
