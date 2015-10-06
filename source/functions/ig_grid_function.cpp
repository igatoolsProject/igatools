//-+--------------------------------------------------------------------
// Igatools a general purpose Isogeometric analysis library.
// Copyright (C) 2012-2015  by the igatools authors (see authors.txt).
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

#include <igatools/functions/ig_grid_function.h>
#include <igatools/functions/ig_grid_function_handler.h>


IGA_NAMESPACE_OPEN


template<int dim,int space_dim>
auto
IgGridFunction<dim,space_dim>::
create_cache_handler() const -> std::shared_ptr<typename parent_t::ElementHandler>
{
  return std::make_shared<IgGridFunctionHandler<dim,space_dim>>(
    std::dynamic_pointer_cast<const self_t>(this->shared_from_this()));
}



template<int dim,int space_dim>
auto
IgGridFunction<dim,space_dim>::
get_ig_space() const -> std::shared_ptr<IgSpace>
{
  return ig_space_;
}

template<int dim,int space_dim>
const IgCoefficients &
IgGridFunction<dim,space_dim>::
get_coefficients() const
{
  return coeffs_;
}


IGA_NAMESPACE_CLOSE


#include <igatools/functions/ig_grid_function.inst>
