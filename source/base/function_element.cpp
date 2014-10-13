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


#include <igatools/base/function_element.h>


IGA_NAMESPACE_OPEN

template<int dim, int range, int rank>
auto
FunctionElement<dim, range, rank>::
get_points() const -> ValueVector<Point>
{
    return CartesianGridElement<dim>::get_points();
}



template<int dim, int range, int rank>
auto
FunctionElement<dim, range, rank>::
get_values() const -> ValueVector<Value> const &
{
    return elem_cache_->values_;
}



template<int dim, int range, int rank>
template<int order>
auto const &
FunctionElement<dim, range, rank>::
get_derivative() const
{
    return std::get<order>(elem_cache_->derivatives_);
}



template<int dim, int range, int rank>
auto
FunctionElement<dim, range, rank>::
get_gradients() const -> ValueVector<Gradient> const &
{
    return get_derivative<1>();
}



template<int dim, int range, int rank>
auto
FunctionElement<dim, range, rank>::
get_hessians() const -> ValueVector<Hessian> const &
{
    return get_derivative<2>();
}

IGA_NAMESPACE_CLOSE

