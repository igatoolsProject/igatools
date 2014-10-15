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

#ifndef NEW_PUSH_FORWARD_ELEMENT_ACCESSOR_H_
#define NEW_PUSH_FORWARD_ELEMENT_ACCESSOR_H_

#include <igatools/geometry/new_push_forward.h>
#include <igatools/geometry/new_mapping_element_accessor.h>

IGA_NAMESPACE_OPEN


template<Transformation type, int dim, int codim = 0>
class PushForwardElement
    : public MappingElement<dim, codim>
{
private:
    using self_t  = PushForwardElement<type, dim, codim>;
    using paren_t = MappingElement<dim, codim>;
    using PF      = NewPushForward<type, dim, codim>;
public:
   // using ContainerType = const PF;
    using paren_t::MappingElement;

#if 0
    template<int order>
    using InvDerivative = typename Map::template InvDerivative<order>;

    ValueVector<Real> const &get_measures() const;

    ValueVector<Real> const &get_w_measures() const;

    ValueVector<InvDerivative<1>> const &get_inverse_gradients() const;

    ValueVector<InvDerivative<2>> const &get_inverse_hessians() const;
#endif
private:


private:
    template <typename Accessor> friend class GridForwardIterator;
    friend class NewPushForward<type, dim, codim>;

};

IGA_NAMESPACE_CLOSE

#endif
