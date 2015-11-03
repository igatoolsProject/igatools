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

#include <igatools/geometry/formula_grid_function_handler.h>

IGA_NAMESPACE_OPEN

template<int dim, int space_dim>
FormulaGridFunctionHandler<dim, space_dim>::
FormulaGridFunctionHandler(const std::shared_ptr<GridFunctionType> &grid_function)
  :
  parent_t::GridFunctionHandler(grid_function)
{}



template<int dim_, int space_dim_>
auto
FormulaGridFunctionHandler<dim_, space_dim_>::
set_flags(const topology_variant &sdim,
          const Flags &flag) -> void
{
  this->get_grid_handler().set_flags(sdim, grid_element::Flags::point);
  parent_t::set_flags(sdim, flag);
}



template<int dim, int space_dim>
auto
FormulaGridFunctionHandler<dim, space_dim>::
fill_cache(const topology_variant &sdim,
           ElementAccessor &elem,
           const int s_id) const  -> void
{
  parent_t::fill_cache(sdim, elem, s_id);

//  GridFunctionType &grid_function = *std::dynamic_pointer_cast<GridFunctionType>(this->get_grid_function());
//  auto disp = FillCacheDispatcher(grid_function, *this, elem, s_id);

  auto disp = FillCacheDispatcher(*this, elem, s_id);
  boost::apply_visitor(disp, sdim);

}

IGA_NAMESPACE_CLOSE

#include <igatools/geometry/formula_grid_function_handler.inst>
