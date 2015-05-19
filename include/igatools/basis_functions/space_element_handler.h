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


#ifndef SPACE_ELEMENT_HANDLER_H_
#define SPACE_ELEMENT_HANDLER_H_

#include <igatools/base/config.h>
#include <igatools/basis_functions/space.h>
#include <igatools/base/tuple_utils.h>

IGA_NAMESPACE_OPEN



/**
 *
 * @ingroup serializable
 */
template <int dim,int codim,int range,int rank>
class SpaceElementHandler
{
public:
    using ElementAccessor = typename Space<dim,codim,range,rank>::ElementAccessor;
    using ElementIterator = typename Space<dim,codim,range,rank>::ElementIterator;

private:
    using eval_pts_variant = SubElemVariants<Quadrature,dim>;
    using topology_variant = TopologyVariants<dim>;


protected:
    /** @name Constructors */
    ///@{
    /**
     * Default constructor. It does nothing but it is needed for the
     * <a href="http://www.boost.org/doc/libs/release/libs/serialization/">boost::serialization</a>
     * mechanism.
     */
    SpaceElementHandler() = default;


    SpaceElementHandler(std::shared_ptr<const Space<dim,codim,range,rank>> space);

    /**
     * Copy constructor. Not allowed to be used.
     */
    SpaceElementHandler(const SpaceElementHandler<dim,codim,range,rank> &elem_handler) = delete;

    /**
     * Move constructor. Not allowed to be used.
     */
    SpaceElementHandler(SpaceElementHandler<dim,codim,range,rank> &&elem_handler) = delete;

public:

    /**
     * Destructor.
     */
    virtual ~SpaceElementHandler() = default;

    ///@}

public:

    /**
     * Resets all the internal data in order to use the
     * same quadrature scheme for each active element of the space.
     */
    void reset(const ValueFlags &flag, const eval_pts_variant &quad);


    /**
     * Resets all the internal data in order to use the
     * same quadrature scheme for the elements of the space with ID specified by
     * the input parameter <tt>elements_flat_id</tt>.
     *
     * @note This function is pure virtual and must be implemented in the class that are derived
     * from SpaceElementHandler.
     */
    virtual void reset_selected_elements(
        const ValueFlags &flag,
        const eval_pts_variant &eval_points,
        const SafeSTLVector<int> &elements_flat_id) = 0;


    virtual void init_cache(SpaceElement<dim,codim,range,rank> &elem,
                            const topology_variant &topology) = 0;

    template <int sub_elem_dim>
    void init_cache(SpaceElement<dim,codim,range,rank> &elem);

    /**
     * Allocates the space in the cache of ElementAccessor <tt>element</tt>
     * necessary for the given quadrature and flag combination.
     * It also fills the invariant (not changing) members of
     * the cache.
     */
    void init_element_cache(ElementAccessor &elem);

    void init_element_cache(ElementIterator &elem);

    void init_face_cache(ElementAccessor &elem);

    void init_face_cache(ElementIterator &elem);


    virtual void fill_cache(
        SpaceElement<dim,codim,range,rank> &elem,
        const topology_variant &topology,
        const int sub_elem_id) = 0;

    template<int sub_elem_dim>
    void fill_cache(ElementAccessor &elem, const int sub_elem_id);

    void fill_element_cache(ElementAccessor &elem);

    void fill_element_cache(ElementIterator &elem);

    void fill_face_cache(ElementAccessor &elem, const int face_id);

    void fill_face_cache(ElementIterator &elem, const int face_id);


    virtual void print_info(LogStream &out) const = 0;


    std::shared_ptr<const Space<dim,codim,range,rank>> get_space() const;

private:
    std::shared_ptr<const Space<dim,codim,range,rank>> space_;


#ifdef SERIALIZATION
    /**
     * @name Functions needed for boost::serialization
     * @see <a href="http://www.boost.org/doc/libs/release/libs/serialization/">boost::serialization</a>
     */
    ///@{
    friend class boost::serialization::access;

    template<class Archive>
    void
    serialize(Archive &ar, const unsigned int version);
    ///@}
#endif // SERIALIZATION
};

IGA_NAMESPACE_CLOSE

#endif // SPACE_ELEMENT_HANDLER_H_
