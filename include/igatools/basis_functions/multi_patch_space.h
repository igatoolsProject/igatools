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

#ifndef MULTI_PATCH_SPACE_H_
#define MULTI_PATCH_SPACE_H_

#include <igatools/base/config.h>
#include <igatools/base/logstream.h>
#include <igatools/utils/concatenated_forward_iterator.h>
#include <igatools/basis_functions/dofs_manager.h>

#include <boost/graph/adjacency_list.hpp>

#include <memory>

IGA_NAMESPACE_OPEN




/**
 * @brief This class represents a space built upon several patches, where each patch is a different
 * instance of a given PhysicalSpace type.
 *
 * In order to build/represent a multi-patch space, two ingredients are needed:
 * - a vector of patches;
 * - a vector of interfaces, i.e. a data structure that specifies if and how the patches are
 * ''glued'' together.
 *
 * By now, the building process of a MultiPatchSpace must follow this order:
 * - 1) insertion of the patches (mandatory);
 * - 2) insertion of the interfaces (optional).
 * If no interfaces are inserted, the patches in MultiPatchSpace will be considered as totally
 * independent one from the other.
 *
 * @note Some restrictions must be ensured for building a MultiPatchSpace:
 * - each patch is defined by one and only one PhysicalSpace object (i.e. a given PhysicalSpace object
 * cannot be used to define more than one patch).
 * - each PhysicalSpace object must be build by one and only one RefSpace object
 * (i.e. a given RefSpace object cannot be used to define more than one PhysicalSpace).
 * - each PhysicalSpace object must be build by one and only one PushFwd object
 * (i.e. a given PushFwd object cannot be used to define more than one PhysicalSpace).
 * - each PushFwd object must be build by one and only one Mapping object
 * (i.e. a given Mapping object cannot be used to define more than one PushFwd).
 *
 * @ingroup containers
 */
template <class PhysicalSpace>
class MultiPatchSpace
{
public:
    /** @names Type aliases used within this class */
    ///@{
    /** Type alias for the reference space. */
    using RefSpace = typename PhysicalSpace::RefSpace;

    /** Type alias for the push-forward. */
    using PushForward = typename PhysicalSpace::PushForwardType;

    /** Type alias for the mapping. */
    using Map = typename PushForward::Map;

    /** Type alias for a patch . */
    using Patch = PhysicalSpace;

    /** Type alias for the pointer to a patch . */
    using PatchPtr = std::shared_ptr<const Patch>;

    /** Dimensionality of the reference domain. */
    static const int dim = PhysicalSpace::dim;
    ///@}


    /** @name Constructors */
    ///@{
    /** Default constructor. */
    MultiPatchSpace() = default;

    /** Copy constructor. */
    MultiPatchSpace(const MultiPatchSpace<PhysicalSpace> &multi_patch_space) = delete;

    /** Move constructor. */
    MultiPatchSpace(MultiPatchSpace<PhysicalSpace> &&multi_patch_space) = delete;

    /** Destructor. */
    ~MultiPatchSpace() = default;
    ///@}


    /** @name Assignment operators */
    ///@{
    /** Copy assignment operator.*/
    MultiPatchSpace<PhysicalSpace> &operator=(const MultiPatchSpace<PhysicalSpace> &multi_patch_space) = delete;

    /** Move assignment operator.*/
    MultiPatchSpace<PhysicalSpace> &operator=(MultiPatchSpace<PhysicalSpace> &&multi_patch_space) = delete;
    ///@}

    /** @name Functions for the management of the patches and/or interfaces addition. */
    ///@{

    /**
     * Sets the object in a state that permits the user to add new patches/interfaces.
     */
    void arrangement_open();


    /**
     * Communicates that the insertion of patches/interfaces is completed.
     *
     * Moreover, performs the data analysis in order to set equality and linear constraints
     * for the degrees of freedom.
     *
     * @warning After calling this function, it will be not possible to add new patches/interfaces
     * to the MultiPatchSpace.
     */
    void arrangement_close();

    /**
     * Adds a patch to the space.
     * @note In Debug mode, an assertion will be raised if the patch is already added.
     * @pre Before calling this function, the member variable is_arrangement_open_ must be set to TRUE
     * (for example using the function arrangement_open()).
     */
    void add_patch(PatchPtr patch);

    /**
     * Adds an interface between two different patches.
     * @note In Debug mode, an assertion will be raised if the interface is already added.
     * @pre Before calling this function, the member variable is_arrangement_open_ must be set to TRUE
     * (for example using the function arrangement_open()).
     */
    void add_interface(const InterfaceType &type,
                       PatchPtr patch_0,const int side_id_patch_0,
                       PatchPtr patch_1,const int side_id_patch_1);
    ///@}


    /** Returns the number of patches used to define this space. */
    int get_num_patches() const;


    /** Returns the number of interfaces used to define this space. */
    int get_num_interfaces() const;



    /** Returns a const reference to the DofsManager used in the MultiPatchSpace. */
    const DofsManager &get_dofs_manager() const;


    /**
     * Prints internal information about the space.
     * @note Mostly used for debugging and testing.
     */
    void print_info(LogStream &out) const;

private:


    bool is_arrangement_open_ = false;

    /** Vector of patches defining the MultiPatchSpace. */
    std::vector< PatchPtr > patches_;

    /**
     * Add the proper offset to the dofs id in the reference spaces in order to avoid same
     * dof ids between different spaces.
     *
     * @warning This function modifies some internal variables of the reference spaces used to
     * define the physical spaces (i.e. the patches).
     */
    void perform_ref_spaces_add_dofs_offset();


    /**
     * @brief This class represent an interface between two patches.
     *
     * An interface is defined by the two patches and one side/face indicator for each patch, plus
     * a flag of type InterfaceType used to specify the type of interface.
     *
     */
    class Interface
    {
    public:
        /** @name Constructors and destructor */
        ///@{
        /** Default constructor. */
        Interface() = delete;

        Interface(const InterfaceType &type,
                  PatchPtr patch_0,const int side_id_patch_0,
                  PatchPtr patch_1,const int side_id_patch_1);


        /** Copy constructor. */
        Interface(const Interface &interface) = delete;


        /** Move constructor. */
        Interface(Interface &&interface) = delete;


        /** Destructor. */
        ~Interface() = default;
        ///@}


        /** @name Assignment operators */
        ///@{
        /** Copy assignment operator. */
        Interface &operator=(const Interface &interface) = delete;


        /** Move assignment operator. */
        Interface &operator=(Interface &&interface) = delete;
        ///@}


        /** @name Comparison operators */
        ///@{
        /** Compares for equality. */
        bool operator==(const Interface &interface_to_compare) const;

        /** Compares for inequality. */
        bool operator!=(const Interface &interface_to_compare) const;
        ///@}

        /**
         * Prints internal information about the Interface.
         * @note Mostly used for debugging and testing.
         */
        void print_info(LogStream &out) const;

        /** Interface type. */
        InterfaceType type_;

        /**
         * The two pairs patch/side_id defining the interface.
         */
        std::array<std::pair<PatchPtr,int>,2> patch_and_side_;
    };

    /** Type alias for the pointer to an interface. */
    using InterfacePtr = std::shared_ptr<const Interface>;

    std::vector<InterfacePtr> interfaces_;






    /** @name Stuff related to the multipatch graph(implemented with the boost::graph library) */
    ///@{
    /** Type of container used to hold the edges of the graph. */
    using out_edge_list_t = boost::vecS;

    /** Type of container used to hold the vertices of the graph. */
    using vertex_list_t = boost::vecS;

    /** A vertex of the graph represents a Patch. */
    using vertex_property_t = PatchPtr;

    /** An edge of the graph represents an interface. */
    using edge_property_t = InterfacePtr;

    /** Type of the graph. */
    using Graph = boost::adjacency_list<
                  out_edge_list_t,
                  vertex_list_t,
                  boost::undirectedS,
                  vertex_property_t,
                  edge_property_t> ;

    /** Graph container used to represent the tree of the elements. */
    Graph multipatch_graph_;
    ///*}


    DofsManager dofs_manager_;
};


IGA_NAMESPACE_CLOSE


#endif // #ifndef MULTI_PATCH_SPACE_H_
