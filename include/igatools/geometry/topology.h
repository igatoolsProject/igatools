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


#ifndef TOPOLOGY_H_
#define TOPOLOGY_H_

#include <igatools/base/config.h>

#include <vector>

IGA_NAMESPACE_OPEN


/**
 * @brief Base class for topology identification.
 *
 * Many quantities in the library can be computed at quadrature points defined on the interior of
 * an element or at quadrature points on the interior of an element-face.
 * In order to unify the interfaces (and then remove a lot of common code),
 * we provide a way to specify (as input argument) from which topology (element or face-id)
 * the quantities must be computed/retrieved.
 *
 * The purpose of this class is to establish an unified way to refer to the different topologies.
 * The mechanism we chosen is based on the value of the member variable TopologyId::id_ :
 * - if TopologyId::id_ == -1 then the topology is the <b>element</b>;
 * - if TopologyId::id_ >= 0 then the topology is the <b>face</b> identified by the id given by
 * the value of TopologyId::id_;
 *
 * For example: @code TopologyId(-1) @endcode refers to the <em>element</em> topology, while
 * @code TopologyId(4) @endcode refers to the <em>face</em> with id number 4.
 *
 * For user convenience, we provide two specializations of TopologyId class:
 * - ElemTopology that represent the <b>element</b> topology (and essentially is the
 * same as TopologyId(-1) )
 * - FaceTopology, that represent the <b>face</b> topology, for which the face-id must be specified
 * in the FaceTopology(const Index face_id) constructor.
 *
 *
 * For example, if we have an PhysicalSpaceElementAccessor object <tt>elem</tt>, we can ask for the
 * basis functions values on the element with the following code:
 * @code{.cpp}
   const auto values_on_element = elem.get_basis_values(ElemTopology());
   @endcode
 * or also without any argument because the ElemTopology() value is the default value
 * argument for this function, i.e.
 * @code{.cpp}
   const auto values_on_element = elem.get_basis_values();
   @endcode
 * To the other hand, if we want the basis functions values on the <b>face</b> with <tt>id == 4</tt>
 * we should write:
 * @code{.cpp}
   const auto values_on_face = elem.get_basis_values(FaceTopology(4));
   @endcode
 * @author M. Martinelli
 * @date 19 March 2014
 */
template< int dim >
class TopologyId
{
public:
    /** @name Constructors */
    ///@{
    /** Default constructor. Not allowed to be used. */
    TopologyId() = delete;

    /**
     * Constructor. The kind of topology depends on the value of @p id:
     * - if <tt>id == -1</tt>, then the built object represent an <b>element</b>;
     * - if <tt>id >= 0</tt>, then the built object represent a <b>face</b> with the given @p id.
     */
    TopologyId(const Index id);

    /** Copy constructor. */
    TopologyId(const TopologyId<dim> &id) = default;

    /** Move constructor. */
    TopologyId(TopologyId<dim> &&id) = default;

    /** Destructor. */
    ~TopologyId() = default;
    ///@}

    /** @name Assignment operators */
    ///@{
    /** Copy assignment operator. */
    TopologyId<dim> &operator=(const TopologyId<dim> &id) = default;

    /** Move assignment operator. */
    TopologyId<dim> &operator=(TopologyId<dim> &&id) = default;
    ///@}

    /** Returns the id. */
    Index get_id() const;

    /** Returns true if the object refers to the <b>element</b> topology. */
    bool is_element() const;

    /** Returns true if the object refers to the <b>face</b> topology. */
    bool is_face() const;

    /**
     * Returns a vector with the indices of the active directions of the topology.
     * - If the topology is the element, then the vector length is equal to @p dim and contains
     * all integers from 0 to @p dim - 1 (included).
     * - If the topology is the face, then the vector length is equal to @p dim -1 and its values
     * depends on the face index (i.e. the input parameter @p face_id in the constructor
     * FaceTopology<dim>::FaceTopology(const Index face_id) ).
     */
    vector<Index> get_active_directions() const;

private:
    Index id_;
};


/**
 * @brief Class for <b>element</b> topology representation.
 *
 * @see TopologyId
 * @author M. Martinelli
 * @date 19 March 2014
 */
template <int dim>
class ElemTopology : public TopologyId<dim>
{
public:
    /** @name Constructors */
    ///@{
    /** Default onstructor. Build an object representing the <b>element</b> topology. */
    ElemTopology();

    /** Copy constructor. */
    ElemTopology(const ElemTopology<dim> &elem) = default;

    /** Move constructor. */
    ElemTopology(ElemTopology<dim> &&elem) = default;

    /** Destructor. */
    ~ElemTopology() = default;
    ///@}


    /** @name Assignment operators */
    ///@{
    /** Copy assignment operator. */
    ElemTopology<dim> &operator=(const ElemTopology<dim> &id) = default;

    /** Move assignment operator. */
    ElemTopology<dim> &operator=(ElemTopology<dim> &&id) = default;
    ///@}

};


/**
 * @brief Class for <b>face</b> topology representation.
 *
 * @see TopologyId
 * @author M. Martinelli
 * @date 19 March 2014
 */
template <int dim>
class FaceTopology : public TopologyId<dim>
{
public:
    /** @name Constructors */
    ///@{
    /** Default constructor. Not allowed to be used. */
    FaceTopology() = delete;

    /** Constructor. Builds an object representing the face with id given by @p face_id. */
    FaceTopology(const Index face_id);

    /** Copy constructor. */
    FaceTopology(const FaceTopology<dim> &face) = default;

    /** Move constructor. */
    FaceTopology(FaceTopology<dim> &&face) = default;

    /** Destructor. */
    ~FaceTopology() = default;
    ///@}

    /** @name Assignment operators */
    ///@{
    /** Copy assignment operator. */
    FaceTopology<dim> &operator=(const FaceTopology<dim> &id) = default;

    /** Move assignment operator. */
    FaceTopology<dim> &operator=(FaceTopology<dim> &&id) = default;
    ///@}

};



IGA_NAMESPACE_CLOSE

#endif // TOPOLOGY_H_

#include <igatools/geometry/topology-inline.h>

