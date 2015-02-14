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


#include <igatools/basis_functions/reference_element.h>
#include <igatools/basis_functions/reference_space.h>
#include <igatools/basis_functions/space_element.h>
#include <igatools/basis_functions/reference_element_handler.h>


IGA_NAMESPACE_OPEN




template <int dim, int range, int rank>
void
ReferenceElement<dim, range, rank>::
move_to(const Index flat_index)
{
    parent_t::move_to(flat_index);
}

template <int dim, int range, int rank>
template <int deriv_order>
auto
ReferenceElement<dim, range, rank>::
evaluate_basis_derivatives_at_points(const EvaluationPoints<dim> &points) ->
ValueTable<
Conditional< deriv_order==0,
             Value,
             Derivative<deriv_order> > >
{
    auto elem_handler = ReferenceElementHandler<dim,range,rank>::create(this->get_space());

    ValueFlags flags;
    if (deriv_order == 0)
        flags = ValueFlags::value;
    else if (deriv_order == 1)
        flags = ValueFlags::gradient;
    else if (deriv_order == 2)
        flags = ValueFlags::hessian;
    else
    {
        Assert(false,ExcNotImplemented());
    }

    elem_handler->reset_one_element(flags,points,this->get_flat_index());
    elem_handler->template init_cache<dim>(*this);
    elem_handler->template fill_cache<dim>(*this,0);

//    Assert(false,ExcNotImplemented());

    return this->template get_values<deriv_order,dim>(0);
}

IGA_NAMESPACE_CLOSE

#include <igatools/basis_functions/reference_element.inst>


