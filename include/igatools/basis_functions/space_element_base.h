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


#ifndef SPACE_ELEMENT_BASE_H_
#define SPACE_ELEMENT_BASE_H_

#include <igatools/base/config.h>
#include <igatools/geometry/cartesian_grid_element.h>
#include <igatools/basis_functions/function_space.h>


IGA_NAMESPACE_OPEN


template <class Grid> class FunctionSpaceOnGrid;

//template <class Accessor> class CartesianGridIterator;

template <int dim>
class SpaceElementBase : private CartesianGridElement<dim>
{
protected:
    using base_t = CartesianGridElement<dim>;
private:
    using self_t = SpaceElementBase<dim>;

    using Space = FunctionSpaceOnGrid<CartesianGrid<dim>>;

protected:
    std::shared_ptr<const Space> space_;

public:
    using base_t::get_flat_index;
    using base_t::get_tensor_index;
    using base_t::get_grid;
    using base_t::is_boundary;

    using Topology = typename base_t::Topology;

    /** @name Constructors */
    ///@{
    /**
     * Default constructor. Not allowed to be used.
     */
    SpaceElementBase() = delete;

    /**
     * Constructs an accessor to element number index of a
     * function space.
     */
    SpaceElementBase(const std::shared_ptr<const Space> space,
                     const Index elem_index);



    /**
     * Copy constructor.
     * It can be used with different copy policies (i.e. deep copy or shallow copy).
     * The default behaviour (i.e. using the proper interface of a classic copy constructor)
     * uses the deep copy.
     */
    SpaceElementBase(const self_t &elem,
                     const CopyPolicy &copy_policy = CopyPolicy::deep);

    /**
     * Move constructor.
     */
    SpaceElementBase(self_t &&elem) = default;

    /**
     * Destructor.
     */
    ~SpaceElementBase() = default;
    ///@}


    /** Return a reference to "*this" as being an object of type CartesianGridElementAccessor.*/
    CartesianGridElement<dim> &as_cartesian_grid_element_accessor();


    /** Return a const-reference to "*this" as being an object of type CartesianGridElementAccessor.*/
    const CartesianGridElement<dim> &as_cartesian_grid_element_accessor() const;

    void print_info(LogStream &out) const;

    void print_cache_info(LogStream &out) const;

    void
    copy_from(const SpaceElementBase<dim> &elem,const CopyPolicy &copy_policy);


    /**
     * @name Functions for getting information about the element connectivity.
     */
    ///@{
    /**
     * Returns the global dofs of the local (non zero) basis functions
     * on the element.
     *
     * @note The dofs can be filtered invoking the function with the argument @p dof_property.
     * If @p dof_property is equal to DofProperties::active, then no filter is applied.
     *
     * For example:
     * \code
       auto loc_to_glob_all = elem->get_local_to_global(DofProperties::active);
       // loc_to_glob_all[0] is the global id of the first basis function on the element
       // loc_to_glob_all[1] is the global id of the second basis function on the element
       // ...
       auto loc_to_glob_active = elem->get_local_to_global(DofProperties::active);
       // loc_to_glob_active[0] is the global id of the first active basis function on the element
       // loc_to_glob_active[1] is the global id of the second active basis function on the element
       // ...
      \endcode
     *
     */
    vector<Index>
    get_local_to_global(const std::string &dofs_property) const;

    /**
     * Returns the patch dofs of the local (non zero) basis functions
     * on the element.
     *
     * @note The dofs can be filtered invoking the function with the argument @p dof_property.
     * If @p dof_property is equal to DofProperties::active, then no filter is applied.
     *
     */
    vector<Index>
    get_local_to_patch(const std::string &dofs_property) const;


    vector<Index>
    get_local_dofs(const std::string &dofs_property) const;


    /**
     *  Number of non zero basis functions with the given @p dofs_property,
     *  over the current element.
     */
    Size get_num_basis(const std::string &dofs_property) const;


    /**
     *  Number of non zero basis functions (without any property),
     *  over the current element.
     */
    virtual Size get_num_basis() const = 0;
    ///@}



public:

    /**
     * @name Comparison operators.
     *
     * @brief The comparison operators compares the <em>position</em> of the element in the grid.
     *
     * @warning To be comparable, two SpaceElementBase objects must be defined on the same space
     * (and therefore on the same grid),
     * otherwise an assertion will be raised (in Debug mode).
     */
    ///@{
    /** Returns TRUE if the two elements have the same index on the grid. */
    bool operator==(const self_t &a) const;


    /** Returns TRUE if the two elements have different indices on the grid. */
    bool operator!=(const self_t &a) const;

    /**
     * Returns TRUE if the the index of the element on the left of the operator <tt> < </tt>
     * is smaller than the the index of the element on the right.
     * */
    bool operator<(const self_t &a) const;

    /**
     * Returns TRUE if the the index of the element on the left of the operator <tt> < </tt>
     * is bigger than the the index of the element on the right.
     * */
    bool operator>(const self_t &a) const;
    ///@}
};


IGA_NAMESPACE_CLOSE



#endif // #ifndef SPACE_ELEMENT_BASE_H_

