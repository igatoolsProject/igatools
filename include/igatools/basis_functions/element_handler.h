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



IGA_NAMESPACE_OPEN


template<class ElemContainer>
class ElementHandler
{
public:
    using ElemIterator = typename ElemContainer::ElementIterator;
    using ElemAccessor = typename ElemContainer::ElementAccessor;

    static const int dim = ElemContainer::dim;

    static const int l = iga::max(0, dim-num_sub_elem);

    using v2 = typename seq<Int, l, dim>::type;
    using topology_variant = typename boost::make_variant_over<v2>::type;

    /**
     * @name Fill functions
     */
    ///@{
    virtual void fill_cache(ElemAccessor &elem, const topology_variant &topology, const int j) = 0;

    template<int k>
    void fill_cache(ElemAccessor &elem, const int j)
    {
        this->fill_cache(elem,Int<k>(),j);
    }

    template<int k>
    void fill_cache(ElemIterator &elem, const int j)
    {
        this->template fill_cache<k>(*elem,j);
    }

    void fill_element_cache(ElemAccessor &elem)
    {
        this->template fill_cache<dim>(elem,0);
    }

    void fill_element_cache(ElemIterator &elem)
    {
        this->template fill_cache<dim>(*elem,0);
    }
    ///@}

};


IGA_NAMESPACE_CLOSE

#endif // ELEMENT_HANDLER_H_
