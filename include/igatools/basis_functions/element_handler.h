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
#include <igatools/base/function.h>

template <int,int,int>
class ReferenceElementHandler;

IGA_NAMESPACE_OPEN


template<class ElementContainer>
class ElementHandler
{
public:
    using ElemIterator = typename ElementContainer::ElementIterator;
    using ElemAccessor = typename ElementContainer::ElementAccessor;
    using DerivedElemHandler = typename ElementContainer::ElementHandler;

    static const int dim = ElementContainer::dim;
    static const int range = ElementContainer::range;
    static const int rank = ElementContainer::rank;


    DerivedElemHandler &as_derived_class()
    {
        return static_cast<DerivedElemHandler &>(*this);
    }
    /*
        static const int l = iga::max(0, dim-num_sub_elem);

        using v2 = typename seq<Int, l, dim>::type;
        using topology_variant = typename boost::make_variant_over<v2>::type;
    //*/


    /**
     * @name Init functions
     */
    ///@{
    template <int k>
    void init_cache(ElemIterator &elem)
    {
        this->as_derived_class().template init_cache<k>(*elem);
    }

    void init_element_cache(ElemAccessor &elem)
    {
        this->as_derived_class().template init_cache<dim>(elem);
    }

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
        this->as_derived_class().template fill_cache<dim-1>(elem,j);
    }

    void fill_face_cache(ElemIterator &elem, const int j)
    {
        this->fill_face_cache(*elem,j);
    }

    ///@}
};


IGA_NAMESPACE_CLOSE

#endif // ELEMENT_HANDLER_H_
