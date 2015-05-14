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


#ifndef ELEMENT_HANDLER_H_
#define ELEMENT_HANDLER_H_

#include <igatools/base/config.h>
#include <igatools/basis_functions/space.h>

IGA_NAMESPACE_OPEN


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



template <int dim,int codim,int range,int rank>
class SpaceElementHandler
{
public:
    using ElementAccessor = typename Space<dim,codim,range,rank>::ElementAccessor;
    using ElementIterator = typename Space<dim,codim,range,rank>::ElementIterator;

private:
    using eval_pts_variant = SubElemVariants<Quadrature,dim>;

public:
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


    template <int sub_elem_dim>
    void init_cache(ElementAccessor &elem)
    {
        Assert(false,ExcNotImplemented());
//        this->as_derived_class().template init_cache<sub_elem_dim>(elem);
    }

    template<int sub_elem_dim>
    void fill_cache(ElementAccessor &elem, const int sub_elem_id)
    {
        Assert(false,ExcNotImplemented());
//        this->as_derived_class().template fill_cache<sub_elem_dim>(elem,sub_elem_id);
    }

    virtual void print_info(LogStream &out) const
    {
        Assert(false,ExcNotImplemented());
    }


};

IGA_NAMESPACE_CLOSE

#endif // ELEMENT_HANDLER_H_
