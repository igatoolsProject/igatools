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

template<int dim, int range>
FormulaGridFunctionHandler<dim, range>::
FormulaGridFunctionHandler(const std::shared_ptr<GridFunctionType> &grid_function)
  :
  parent_t::GridFunctionHandler(grid_function)
{}



template<int dim_, int range_>
auto
FormulaGridFunctionHandler<dim_, range_>::
set_flags(const topology_variant &sdim,
          const Flags &flag) -> void
{
  this->get_grid_handler().set_flags(sdim, grid_element::Flags::point);
  parent_t::set_flags(sdim, flag);
}


template<int dim, int range>
FormulaGridFunctionHandler<dim, range>::
FillCacheDispatcher::
FillCacheDispatcher(const self_t &grid_function_handler,
                    ElementAccessor &elem,
                    const int s_id)
  :
  grid_function_handler_(grid_function_handler),
  elem_(elem),
  s_id_(s_id)
{}


template<int dim, int range>
auto
FormulaGridFunctionHandler<dim, range>::
fill_cache(const topology_variant &sdim,
           ElementAccessor &elem,
           const int s_id) const  -> void
{
//  parent_t::fill_cache(sdim, elem, s_id);

  this->grid_handler_->fill_cache(sdim, elem.get_grid_element(), s_id);

//  GridFunctionType &grid_function = *std::dynamic_pointer_cast<GridFunctionType>(this->get_grid_function());
//  auto disp = FillCacheDispatcher(grid_function, *this, elem, s_id);

  auto disp = FillCacheDispatcher(*this, elem, s_id);
  boost::apply_visitor(disp, sdim);

}



template<int dim, int range>
template<int sdim>
void
FormulaGridFunctionHandler<dim, range>::
FillCacheDispatcher::
operator()(const Topology<sdim> &sub_elem)
{
  const auto &formula_grid_function =
    *std::dynamic_pointer_cast<GridFunctionType>(grid_function_handler_.get_grid_function());


  auto &local_cache = grid_function_handler_.get_element_cache(elem_);
  auto &cache = local_cache.template get_sub_elem_cache<sdim>(s_id_);

  if (!cache.fill_none())
  {
    const auto &grid_pts = elem_.get_grid_element().template get_points<sdim>(s_id_);
    if (cache.template status_fill<_D<0>>())
    {
      auto &F = cache.template get_data<_D<0>>();
      formula_grid_function.evaluate_0(grid_pts, F);
      F.set_status_filled(true);
    }

    if (cache.template status_fill<_D<1>>())
    {
      auto &DF = cache.template get_data<_D<1>>();
      formula_grid_function.evaluate_1(grid_pts, DF);
      DF.set_status_filled(true);
    }

    if (cache.template status_fill<_D<2>>())
    {
      auto &D2F = cache.template get_data<_D<2>>();
      formula_grid_function.evaluate_2(grid_pts, D2F);
      D2F.set_status_filled(true);
    }
//        if (cache.template status_fill<_Divergence>())
//          Assert(false,ExcNotImplemented());
  }

  cache.set_filled(true);
}

IGA_NAMESPACE_CLOSE

#include <igatools/geometry/formula_grid_function_handler.inst>
