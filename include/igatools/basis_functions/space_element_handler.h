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

#if 0
template<class SpaceType>
class ElementHandler
{
public:
    using ElemIterator = typename SpaceType::ElementIterator;
    using ElemAccessor = typename SpaceType::ElementAccessor;
    using DerivedElemHandler = typename SpaceType::ElementHandler;

    static const int dim = SpaceType::dim;
    static const int range = SpaceType::range;
    static const int rank = SpaceType::rank;


    DerivedElemHandler &as_derived_class()
    {
        return static_cast<DerivedElemHandler &>(*this);
    }

    /**
     * @name Reset functions
     */
    ///@{
    ///@}

    /**
     * @name Init functions
     */
    ///@{
    template <int k>
    void init_cache(ElemIterator &elem)
    {
        this->as_derived_class().template init_cache<k>(*elem);
    }

    /**
     * Allocates the space in the cache of ElementAccessor <tt>element</tt>
     * necessary for the given quadrature and flag combination.
     * It also fills the invariant (not changing) members of
     * the cache.
     */
    void init_element_cache(ElemAccessor &elem)
    {
        this->as_derived_class().template init_cache<dim>(elem);
    }

    /**
     * Same as init_element_cache() but using the ElementIterator as input/output argument.
     *
     * @sa init_element_cache(ElemAccessor &elem)
     */
    void init_element_cache(ElemIterator &elem)
    {
        this->init_element_cache(*elem);
    }

    void init_face_cache(ElemAccessor &elem)
    {
        this->as_derived_class().template init_cache<(dim > 0)?dim-1:0>(elem);
    }

    void init_face_cache(ElemIterator &elem)
    {
        this->init_face_cache(*elem);
    }
    ///@}

    /**
     * @name Fill functions
     */
    ///@{
    template<int k>
    void fill_cache(ElemIterator &elem, const int j)
    {
        this->as_derived_class().template fill_cache<k>(*elem,j);
    }

    void fill_element_cache(ElemAccessor &elem)
    {
        this->as_derived_class().template fill_cache<dim>(elem,0);
    }

    void fill_element_cache(ElemIterator &elem)
    {
        this->fill_element_cache(*elem);
    }

    void fill_face_cache(ElemAccessor &elem, const int j)
    {
        this->as_derived_class().template fill_cache<(dim > 0)?dim-1:0>(elem,j);
    }

    void fill_face_cache(ElemIterator &elem, const int j)
    {
        this->fill_face_cache(*elem,j);
    }
    ///@}
};
#endif


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
    		const topology_variant & topology) = 0;

    template <int sub_elem_dim>
    void init_cache(SpaceElement<dim,codim,range,rank> &elem)
    {
    	this->init_cache(elem,Topology<sub_elem_dim>());
    }

    /**
     * Allocates the space in the cache of ElementAccessor <tt>element</tt>
     * necessary for the given quadrature and flag combination.
     * It also fills the invariant (not changing) members of
     * the cache.
     */
    void init_element_cache(SpaceElement<dim,codim,range,rank> &elem)
    {
        this->template init_cache<dim>(elem);
    }

    void init_element_cache(ElementIterator &elem)
    {
        this->template init_cache<dim>(*elem);
    }


    virtual void fill_cache(
    		SpaceElement<dim,codim,range,rank> &elem,
    		const topology_variant & topology,
			const int sub_elem_id) = 0;

    template<int sub_elem_dim>
    void fill_cache(SpaceElement<dim,codim,range,rank> &elem, const int sub_elem_id)
    {
    	this->fill_cache(elem,Topology<sub_elem_dim>(),sub_elem_id);
    }

    void fill_element_cache(SpaceElement<dim,codim,range,rank> &elem)
    {
        this->template fill_cache<dim>(elem,0);
    }

    void fill_element_cache(ElementIterator &elem)
    {
        this->template fill_cache<dim>(*elem,0);
    }


    virtual void print_info(LogStream &out) const
    {
        Assert(false,ExcNotImplemented());
    }

    template <int sub_elem_dim = dim>
    Size get_num_points() const
    {
        Assert(false,ExcNotImplemented());
//        return grid_handler_.template get_num_points<sub_elem_dim>();
        return 0;
    }

    std::shared_ptr<const Space<dim,codim,range,rank>> get_space() const;

private:
    std::shared_ptr<const Space<dim,codim,range,rank>> space_;

//#if 0
#ifdef SERIALIZATION
    /**
     * @name Functions needed for boost::serialization
     * @see <a href="http://www.boost.org/doc/libs/release/libs/serialization/">boost::serialization</a>
     */
    ///@{
    friend class boost::serialization::access;

    template<class Archive>
    void
    serialize(Archive &ar, const unsigned int version)
    {
        auto non_const_space = std::const_pointer_cast<Space<dim,codim,range,rank>>(space_);
        ar &boost::serialization::make_nvp("space_", non_const_space);
        space_ = non_const_space;
        Assert(space_ != nullptr,ExcNullPtr());
    }
    ///@}
#endif // SERIALIZATION
//#endif
};

IGA_NAMESPACE_CLOSE

#endif // SPACE_ELEMENT_HANDLER_H_
