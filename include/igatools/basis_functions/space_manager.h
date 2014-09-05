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
#include <igatools/base/linear_constraint.h>
#include <igatools/utils/concatenated_iterator.h>
#include <igatools/basis_functions/equality_constraint.h>
#include <igatools/basis_functions/nurbs_space.h>
#include <igatools/basis_functions/physical_space.h>


#include <igatools/contrib/variant.h>
#include <boost/variant.hpp>

#include <memory>

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
    void space_insertion_open();

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
    void space_insertion_close(const bool automatic_dofs_renumbering = true);


    /**
     * Adds a space to the SpaceManager.
     *
     * @note An assertion will be raised (in DEBUG mode)
     * if the passed <p>space</p> is already present in the SpaceManager.
     */
    template<class Space>
    void add_space(std::shared_ptr<Space> space);
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

    /** @name Functions for the insertion of the linear constraints */
    ///@{
    /**
     * Sets the SpaceManager in a state that can receive new linear constraints.
     */
    void linear_constraints_open();


    /**
     * Add a LinearConstraint to the SpaceManager,
     * where @p dofs are the dofs id involved by the constraint,
     * @p coeffs their coefficients and
     * @p rhs is the right hand side that defines the linear constraint equation.
     *
     * @note An assertion will be raised (in DEBUG mode)
     * if the space manager is not set in the proper state by the function
     * linear_constraints_open().
     */
    void add_linear_constraint(const vector<Index> &dofs, const vector<Real> &coeffs, const Real rhs);


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

    /**
     * Returns true if the space insertion is open.
     */
    bool is_space_insertion_open() const;


    void add_element_dofs_view(const DofsConstView &element_dofs_view);


private:
    bool is_space_insertion_open_ = false;
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
        SpaceInfo();
        SpaceInfo(const SpacePtrVariant &space,
                  const Index num_dofs,
                  const Index min_dofs_id,
                  const Index max_dofs_id,
                  const DofsView &dofs_view,
                  const std::shared_ptr<const std::map<Index,DofsConstView>> elements_dofs_view);

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

    private:
        /**
         * Pointer to a generic single-patch space (it can be any of the type allowed for BSplineSpace,
         * NURBSSpace and PhysicalSpace).
         */
        SpacePtrVariant space_;

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
         * @note We use a std:shared_ptr because this container can be very big and
         * it is already present
         * the the DofDistribution class instantiated in the space itself.
         */
        std::shared_ptr<const std::map<Index,DofsConstView>> elements_dofs_view_;
    };

    /**
     * Map containing the pointers to the spaces handled by the SpaceManager,
     * and some useful informations that does not depends on the template
     * parameters needed to instantiate the spaces.
     */
    std::map<int,SpaceInfo> spaces_info_;


    /**
     * View to the active dofs ids of all single-patch spaces handled by the SpaceManager.
     */
    DofsView dofs_view_;


    vector<std::shared_ptr<LinearConstraint>> linear_constraints_;


    vector<EqualityConstraint> equality_constraints_;


    /** Counts and return the number of unique dofs in the SpaceManager. */
    Index count_unique_dofs() const;

    /** Number of unique dofs in the SpaceManager. */
    Index num_unique_dofs_;


public:
    /**
     * Returns a map containing the pointers to the spaces handled by the SpaceManager,
     * and some useful informations that does not depends on the template
     * parameters needed to instantiate the spaces.
     */
    const std::map<Index,SpaceInfo> &get_spaces_info() const;
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
        Assert(space != get<std::shared_ptr<Space>>(space_info.second.get_space_variant()),
               ExcMessage("Space already added in the SpaceManager."));
    }
#endif

    //------------------------------------------------------------------------
    using RefSpace = typename Space::RefSpace;
    auto ref_space = std::const_pointer_cast<RefSpace>(space->get_reference_space());

    auto &dof_distribution = ref_space->get_dof_distribution_global();

    spaces_info_[ref_space->get_id()] =
        SpaceInfo(space,
                  ref_space->get_num_basis(),
                  dof_distribution.get_min_dof_id(),
                  dof_distribution.get_max_dof_id(),
                  dof_distribution.get_dofs_view(),
                  dof_distribution.get_elements_view());
    //---------------------------------------------------------------------------------------------
}


IGA_NAMESPACE_CLOSE


#endif // #ifndef SPACE_MANAGER_H_
