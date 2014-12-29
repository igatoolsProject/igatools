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
#include <igatools/basis_functions/reference_space.h>
#include <igatools/basis_functions/physical_space.h>


#include <igatools/contrib/variant.h>
#include <boost/variant.hpp>

#include <memory>
#include <set>


// TODO (pauletti, Oct 3, 2014): the next line is very developer stage
// should not be in this format for main branch
// Update the naming and the pyhton script as well as the install option
/**
 * Include file used for declare the Variant object that can take
 */
#include <igatools/basis_functions/space_manager.inst>


IGA_NAMESPACE_OPEN

template <class RefSpace, int codim, Transformation type>
using PhysSpacePtr = std::shared_ptr<PhysicalSpace<RefSpace,codim,type>>;





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
 *      add_spaces_connection().
 *      Please note that this last function can be used with two different signatures:
 *        - the first version takes two arguments for the spaces: one argument is used to define the
 *        space of test functions (i.e. the space defining the block of dofs along the rows)
 *        and the other is used to define the space of trial functions
 *        (i.e. the space defining the block of dofs along the columns).
 *        This means that if a user wants to add symmetric blocks, he needs to call this function twice with
 *        the spaces reverted.
 *        - the second version takes only one space: this can be used to define blocks along the diagonal.
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

        sp_manager.add_spaces_connection(space_0); // block 0-0
        sp_manager.add_spaces_connection(space_1); // block 1-1

        sp_manager.spaces_connectivity_close();  // ends the spaces connectivity phase (phase 2)

        @endcode
 *
 *   3. <b>Definition of the <em>global dofs connectivity</em> in each spaces pair</b>.
 *   Once the relations between the spaces are defined (Phase 2), it remains to define the dofs
 *   relations within a spaces pair.
 *     - For the case in which the SpacesConnection is not defined with a single space
 *   (i.e. we are not in the block-diagonal), the user has the responsability to manage
 *   (and set) the dofs connectivity. In order to do so, the SpacesConnection can be retrieved with
 *   the function get_spaces_connection(), and then the dofs connectivity can be added using
 *   SpacesConnection::add_dofs_connectivity().
 *   The input argument arguments of this last function is a
 *   <tt>std::map<Index,std::set<Index>> </tt>
 *   in which the <tt>key</tt> represents an active global dofs in the block
 *   (belonging from the first space, i.e. a row dof id)
 *   and the <tt>value</tt> is a <tt>std::set<Index></tt> representing the global dofs in the second space
 *   (i.e. a set colum dofs id) that are in relation with the global dof in the first space represented by the key.
 *   @code{.cpp}

     map<Index,set<Index>> dofs_connectivity_0_2;
     ...
     ... // here we fill in some way (this is problem specific) the dofs connectivity between space_0 and space_2
     ...

     auto & conn_0_2 = sp_manager.get_spaces_connection(space_0,space_2); //get the connection between space_0 and space_2
     conn_0_2.add_dofs_connectivity(dofs_connectivity_0_2);

     @endcode
 *
 *     - In the particular case in which a space defines a connection
 *   with itself (i.e. a the connection is represented by a diagonal block),
 *   the dofs relations are automatically (by default) inferred (without taking any action) from the
 *   internal structure of the space, i.e. iterating over the active elements of the space and retrieving
 *   the dofs connectivity in each element.
 *   If this default behaviour is not appropriate to describe the dofs connectivity,
 *   the user must define the space connectivity as follows:
 *   @code{.cpp}

     sp_manager.add_spaces_connection(space_0,false); // the internal dofs strucure of space_0 will not be used
     @endcode
 *    instead of
 *   @code{.cpp}

     sp_manager.add_spaces_connection(space_0); // the internal dofs strucure of space_0 will be used
     @endcode
 *   In any case extra dofs connectivity can be added as explained above, retrieving the SpacesConnection
 *   and then using SpacesConnection::add_dofs_connectivity().
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

    class DofsConnectivity
    {
    public:
        decltype(auto) begin()
        {
            return connectivity_.begin();
        }

        decltype(auto) end()
        {
            return connectivity_.end();
        }

        decltype(auto) begin() const
        {
            return connectivity_.begin();
        }

        decltype(auto) end() const
        {
            return connectivity_.end();
        }

        decltype(auto) operator[](const Index row_global_id)
        {
            return connectivity_[row_global_id];
        }

        void clear()
        {
            connectivity_.clear();
        }

        void merge(const DofsConnectivity &dofs_connectivity)
        {
            for (const auto &dofs_connectivity_map_entry : dofs_connectivity)
            {
                const auto row_dof = dofs_connectivity_map_entry.first;
                const auto &cols_dof = dofs_connectivity_map_entry.second;

                connectivity_[row_dof].insert(cols_dof.begin(),cols_dof.end());
            }
        }

        Index get_num_rows() const
        {
            return connectivity_.size();
        }

        void print_info(LogStream &out) const
        {
            Index counter = 0;
            for (const auto &dofs_connectivity_map_entry : connectivity_)
            {
                const auto row_dof = dofs_connectivity_map_entry.first;
                const auto &cols_dof = dofs_connectivity_map_entry.second;

                out << "row[" << counter++ << ", global id = " << row_dof << "] : [ ";
                for (const auto &col_dof : cols_dof)
                    out << col_dof << " ";
                out << "]" << std::endl;
            }
        }

    private:
        std::map<Index,std::set<Index>> connectivity_;
    };

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


    /** @name Functions for managing the spaces connectivity in the SpaceManager */
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
     * Defines connection between the global dofs of the @p space with itself.
     * If the (optional) argument @p use_dofs_connectivity_from_space is set to FALSE,
     * then the space's internal dofs will not
     * be taken into account to define the dofs connectivity.
     */
    template<class Space>
    void add_spaces_connection(std::shared_ptr<Space> space, const bool use_dofs_connectivity_from_space = true);


    /**
     * Sets the SpaceManager in a state that cannot receive any new connectivity between existing spaces.
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
     * Returns the map of linear constraints of a given @p type, defined in the SpaceManager.
     * If no @type is passed at the input argument, the function returns all the linear constraints.
     */
    std::map<Index,std::shared_ptr<const LinearConstraint> > get_linear_constraints(
        const LinearConstraintType &type = LinearConstraintType::any) const;


    /**
     * This function tests if a vector of global dof values @p dof_values, satisfies the linear constraints
     * (up to the tolerance @p tol).
     * If all coefficients satisfies the linear constraints the function returns an empty vector,
     * otherwise it returns a vector containing the linear constraints that are not satisfied.
     */
    vector<std::shared_ptr<const LinearConstraint> > verify_linear_constraints(
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
     * This class represent a <em>space</em> (ReferenceSpace or PhysicalSpace) without
     * the difficulties (and limits) of the template technology.
     *
     * This class is not templatized, and it stores a certain kind of space through the use of
     * a modified version of the boost::variant class.
     *
     *
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
                  const Transformation transf_type,
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

        /** Returns the dimension of the space.*/
        int get_dim() const;

        /** Returns the codimension of the space.*/
        int get_codim() const;

        /** Returns the space dimension of the space.*/
        int get_space_dim() const;

        /** Returns the range of the space.*/
        int get_range() const;

        /** Returns the rank of the space.*/
        int get_rank() const;

        /** Returns true if all the parameters passed as arguments are the
         *  same as the ones contained as private member.*/
        bool check_parameters(const int dim,
                              const int codim,
                              const int space_dim,
                              const int range,
                              const int rank,
                              const Transformation transf_type) const;

        /** Returns the transformation type of the push forward of the space.*/
        Transformation get_transformation_type() const;


        /**
         * Adds an @p offset to all dofs present in the space.
         */
        void add_dofs_offset(const Index offset);

        /**
         * Return an object containing a variant of a shared_pointer pointing
         * to a certain single-patch space.
         * The allowed space type can be any valid ReferenceSpace or PhysicalSpace.
         */
        SpacePtrVariant &get_space_variant();


        /**
         * Return true if the id of the lhs is equal to the id of the rhs.
         */
        bool operator==(const SpaceInfo &sp) const;


    private:
        /**
         * Pointer to a generic single-patch space (it can be any of the type allowed for
         * ReferenceSpace and PhysicalSpace).
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

        Transformation transf_type_;


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




        std::set<Index> inactive_dofs_;

    public:

        void set_inactive_dofs(const std::set<Index> &inactive_dofs)
        {
            inactive_dofs_ = inactive_dofs;
        }

        const std::set<Index> &get_inactive_dofs() const
        {
            return inactive_dofs_;
        }

    };

public:
    using SpaceInfoPtr = std::shared_ptr<SpaceInfo>;

private:

#if 0
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
#endif


    /**
     * Map containing the pointers to the spaces handled by the SpaceManager,
     * and some useful informations that does not depends on the template
     * parameters needed to instantiate the spaces.
     *
     * The <tt>key</tt> of the std::map is the space id.
     */
    std::map<int,SpaceInfoPtr> spaces_info_;


    /**
     * View to the active dofs ids of all single-patch spaces handled by the SpaceManager.
     */
    DofsView dofs_view_;


    std::map<LCType,std::map<Index,std::shared_ptr<const LC> > > linear_constraints_;


    vector<EqualityConstraint> equality_constraints_;

#if 0
    /** Counts and return the number of unique dofs in the SpaceManager. */
    Index count_unique_dofs() const;


    /** Number of unique dofs in the SpaceManager. */
    Index num_unique_dofs_;
#endif

    /**
     * This function update the dofs view from the underlying spaces.
     *
     * It performs a loop over the spaces and aggregate the partial views in one global view.
     */
    void update_dofs_view();



    /** @name Memebr variables used to keep track of dofs renumbering */
    ///@{
    /** List of spaces that have the original global dofs (after the space creation). */
    std::list<SpaceInfoPtr> spaces_with_original_dofs_;

    /** List of spaces that have the global dofs renumbered. */
    std::list<SpaceInfoPtr> spaces_with_renumbered_dofs_;
    ///@}

    /**
     * This function renumber the dofs of the spaces that are in the list spaces_with_original_dofs_
     * and then put the renumbered spaces in the list spaces_with_renumbered_dofs_
     */
    void perform_space_dofs_renumbering();

    Index space_dofs_offset_;


    class SpacesConnection
    {
    public:
        SpacesConnection() = delete;

        SpacesConnection(const SpaceInfoPtr &space,const bool use_dofs_connectivity_from_space = true);
        SpacesConnection(const SpaceInfoPtr &space_row,const SpaceInfoPtr &space_col);

        ~SpacesConnection() = default;

        SpacesConnection(const SpacesConnection &in) = default;
        SpacesConnection(SpacesConnection &&in) = default;

        SpacesConnection &operator=(const SpacesConnection &in) = delete;
        SpacesConnection &operator=(SpacesConnection &&in) = delete;

        bool operator==(const SpacesConnection &conn) const;

        /** Returns true if the row space and the column space are equal. */
        bool is_unique_space() const;


        void add_dofs_connectivity(const DofsConnectivity &dofs_connectivity);

    private:
        SpaceInfoPtr space_row_;
        SpaceInfoPtr space_col_;

        bool use_dofs_connectivity_from_space_;

        DofsConnectivity extra_dofs_connectivity_;

    public:
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
            std::set<Index> dofs;
            if (use_dofs_connectivity_from_space_)
            {
                const auto &sp_dofs_view = space_row_->get_dofs_view();
                dofs.insert(sp_dofs_view.begin(),sp_dofs_view.end());
            }

            for (const auto &pair : extra_dofs_connectivity_)
                dofs.insert(pair.first);

#if 0
            LogStream out;
            out.begin_item("dofs row: ");
            for (const auto dof : dofs)
                out <<dof << " ";
            out<< std::endl;
            out.end_item();
#endif

#if 0
            const auto &inactive_dofs = space_row_->get_inactive_dofs();

            out.begin_item("inactive dofs row: ");
            for (const auto dof : inactive_dofs)
                out <<dof << " ";
            out<< std::endl;
            out.end_item();

            if (!inactive_dofs.empty())
                dofs.erase(inactive_dofs.begin(),inactive_dofs.end());
#endif

            return dofs;
        }

        std::set<Index> get_col_dofs() const
        {
            std::set<Index> dofs;
            if (use_dofs_connectivity_from_space_)
            {
                const auto &sp_dofs_view = space_col_->get_dofs_view();

                dofs.insert(sp_dofs_view.begin(),sp_dofs_view.end());
            }

            for (const auto &pair : extra_dofs_connectivity_)
                dofs.insert(pair.second.begin(),pair.second.end());

#if 0
            LogStream out;
            out.begin_item("dofs col: ");
            for (const auto dof : dofs)
                out <<dof << " ";
            out<< std::endl;
            out.end_item();
#endif
#if 0
            const auto &inactive_dofs = space_col_->get_inactive_dofs();

            out.begin_item("inactive dofs col: ");
            for (const auto dof : inactive_dofs)
                out <<dof << " ";
            out<< std::endl;
            out.end_item();

            if (!inactive_dofs.empty())
                dofs.erase(inactive_dofs.begin(),inactive_dofs.end());
#endif

            return dofs;
        }

        /**
         * Returns the extra dofs connectivity added to the SpacesConnection
         */
        const DofsConnectivity &get_extra_dofs_connectivity() const
        {
            return extra_dofs_connectivity_;
        }



        bool use_dofs_connectivity_from_space() const
        {
            return use_dofs_connectivity_from_space_;
        }

    };


    vector<SpacesConnection> spaces_connections_;

    /**
     * Extra dofs connectivity that is not related to any SpacesConnection object.
     */
    DofsConnectivity extra_dofs_connectivity_;

public:

    const vector<SpacesConnection> &get_spaces_connections() const
    {
        return spaces_connections_;
    }

    /**
     * Returns the extra dofs connectivity added to the SpaceManager.
     * @note These connection are not related to any SpacesConnection object.
     */
    const DofsConnectivity &get_extra_dofs_connectivity() const
    {
        return extra_dofs_connectivity_;
    }

    /**
     * Add extra dofs connectivity added to the SpaceManager.
     * @note The added dofs connectivity will not be related to any SpacesConnection object.
     */
    void add_dofs_connectivity(const DofsConnectivity &dofs_connectivity);


    /**
     * Returns a map containing the pointers to the spaces handled by the SpaceManager,
     * and some useful informations that does not depends on the template
     * parameters needed to instantiate the spaces.
     */
    const std::map<Index,SpaceInfoPtr> &get_spaces_info() const;


    /**
     * Returns a map containing the pointers to the spaces handled by the SpaceManager,
     * and some useful informations that does not depends on the template
     * parameters needed to instantiate the spaces.
     */
    std::map<Index,SpaceInfoPtr> &get_spaces_info();

    /**
     * Returns the SpacesConnection corresponding to the row space @p space_test and
     * column space @p space_trial.
     */
    template<class SpaceTest,class SpaceTrial>
    SpacesConnection &get_spaces_connection(
        std::shared_ptr<SpaceTest> space_test,
        std::shared_ptr<SpaceTrial> space_trial);


    /**
     * Returns the row dofs id.
     */
    std::set<Index> get_row_dofs() const;

    /**
     * Returns the column dofs id.
     */
    std::set<Index> get_col_dofs() const;

    /**
     * Returns the number of row dofs id.
     */
    Index get_num_row_dofs() const;

    /**
     * Returns the number of column dofs id.
     */
    Index get_num_col_dofs() const;


    /**
     * Returns the sparsity pattern associated to the information stored into the SpaceManager.
     *
     * The sparsity pattern is a <tt>std::map<Index,std::set<Index>></tt> in which the
     * <tt>key-value</tt> pair represents a single row in the system matrix
     * (the <tt>key</tt> is the row dof id and the <tt>value</tt> are the column dofs id connected
     * with the row dof id.
     */
    std::shared_ptr<const DofsConnectivity> get_sparsity_pattern() const;
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
        Assert(space->get_id() != space_info.second->get_id(),
               ExcMessage("The space with id=" + std::to_string(space->get_id()) +
                          " is already added in the SpaceManager."));
    }
#endif

    //------------------------------------------------------------------------
    auto &dof_distribution = space->get_dof_distribution_global();

    auto space_info = std::shared_ptr<SpaceInfo>(
                          new SpaceInfo(space,
                                        space->get_id(),
                                        Space::dim,
                                        Space::codim,
                                        Space::space_dim,
                                        Space::range,
                                        Space::rank,
                                        Space::PushForwardType::type,
                                        space->get_num_basis(),
                                        dof_distribution.get_min_dof_id(),
                                        dof_distribution.get_max_dof_id(),
                                        dof_distribution.get_dofs_view(),
                                        dof_distribution.get_elements_view()));

    spaces_info_[space_info->get_id()] = space_info;

    spaces_with_original_dofs_.push_back(space_info);
    //---------------------------------------------------------------------------------------------
}

template<class SpaceTest,class SpaceTrial>
inline
void
SpaceManager::
add_spaces_connection(std::shared_ptr<SpaceTest> space_test,std::shared_ptr<SpaceTrial> space_trial)
{
    /*
    using std::cout;
    using std::endl;
    cout << "SpaceManager::add_spaces_connection()" << endl;
    cout << "adding connection between space " << space_test->get_id()
         << " and space " << space_trial->get_id() << endl;
    //*/
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

template<class Space>
inline
void
SpaceManager::
add_spaces_connection(std::shared_ptr<Space> space, const bool use_dofs_connectivity_from_space)
{
    Assert(is_spaces_insertion_open_ == false,ExcInvalidState());
    Assert(is_spaces_connectivity_open_ == true,ExcInvalidState());

    Assert(space !=nullptr,ExcNullPtr());

    auto sp = spaces_info_.at(space->get_id());

    SpacesConnection conn(sp,use_dofs_connectivity_from_space);

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
        SpacesConnection(sp_test,sp_trial));

    Assert(it != spaces_connections_.end(),
    ExcMessage("Spaces connection not found."));
    return *it;
}



IGA_NAMESPACE_CLOSE


#endif // #ifndef SPACE_MANAGER_H_
