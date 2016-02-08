//-+--------------------------------------------------------------------
// Igatools a general purpose Isogeometric analysis library.
// Copyright (C) 2012-2016  by the igatools authors (see authors.txt).
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

#include <igatools/functions/ig_grid_function_handler.h>
#include <igatools/functions/ig_grid_function.h>
#include <igatools/functions/grid_function_element.h>
#include <igatools/basis_functions/reference_basis_handler.h>

IGA_NAMESPACE_OPEN

template<int dim_, int range_>
IgGridFunctionHandler<dim_, range_>::
IgGridFunctionHandler(const std::shared_ptr<GridFunctionType> &ig_grid_function)
  :
  parent_t(ig_grid_function),
  ig_basis_handler_(ig_grid_function->get_basis()->create_cache_handler()),
  ig_grid_function_(ig_grid_function)
{}


template<int dim_, int range_>
auto
IgGridFunctionHandler<dim_, range_>::
get_ig_grid_function() const -> std::shared_ptr<GridFunctionType>
{
  return ig_grid_function_;
}



template<int dim_, int range_>
auto
IgGridFunctionHandler<dim_, range_>::
set_flags(const topology_variant &sdim,
          const Flags &flag) -> void
{
//  this->get_grid_handler().set_flags(sdim,flag);

  parent_t::set_flags(sdim,flag);


  using BsFlags = basis_element::Flags;
  BsFlags ig_basis_elem_flags = BsFlags::none;
  if (contains(flag,Flags::D0))
    ig_basis_elem_flags |= BsFlags::value;
  if (contains(flag,Flags::D1))
    ig_basis_elem_flags |= BsFlags::gradient;
  if (contains(flag,Flags::D2))
    ig_basis_elem_flags |= BsFlags::hessian;

  ig_basis_handler_->set_flags_impl(sdim,ig_basis_elem_flags);
  //*/
}


#if 0
template<int dim_, int range_>
void
IgGridFunctionHandler<dim_, range_>::
init_cache(ElementAccessor &elem,
           const eval_pts_variant &quad) const
{
  AssertThrow(false,ExcNotImplemented());
#if 0
  grid_handler_->init_cache(*(elem.grid_elem_), quad);

  auto disp = InitCacheDispatcher(*this, elem, flags_);
  boost::apply_visitor(disp, quad);
#endif
}
#endif


template<int dim_, int range_>
auto
IgGridFunctionHandler<dim_, range_>::
fill_cache(const topology_variant &sdim,
           ElementAccessor &elem,
           const int s_id) const-> void
{
  auto fill_cache_dispatcher =
    FillCacheDispatcher(
      (*this),
      elem,
      s_id);

  boost::apply_visitor(fill_cache_dispatcher, sdim);

//  this->get_grid_handler().fill_cache(sdim, elem.get_grid_element(), s_id);
}


template<int dim_, int range_>
template<int sdim>
void
IgGridFunctionHandler<dim_, range_>::
FillCacheDispatcher::
operator()(const Topology<sdim> &sub_elem)
{
  auto &grid_elem = ig_grid_function_elem_.get_grid_element();
  ig_grid_function_handler_.get_grid_handler().template fill_cache<sdim>(grid_elem,s_id_);


  const auto &ig_grid_function =
    *(std::dynamic_pointer_cast<const IgGridFunction<dim_,range_>>(ig_grid_function_handler_.get_grid_function()));

  auto &local_cache = ig_grid_function_handler_.get_element_cache(ig_grid_function_elem_);
  auto &cache = local_cache.template get_sub_elem_cache<sdim>(s_id_);

  if (!cache.fill_none())
  {
    const auto &grid_elem_id = grid_elem.get_index();

    const auto ig_basis = ig_grid_function.get_basis();
    const auto &ig_basis_handler = *ig_grid_function_handler_.ig_basis_handler_;
    auto ig_basis_elem = ig_basis->begin();
    ig_basis_elem->move_to(grid_elem_id);

    ig_basis_handler.template init_cache<sdim>(*ig_basis_elem,grid_elem.template get_quad<sdim>());
    ig_basis_handler.template fill_cache<sdim>(*ig_basis_elem,s_id_);


    const auto &dofs_property = ig_grid_function.get_dofs_property();

    const auto &ig_basis_elem_global_dofs = ig_basis_elem->get_local_to_global(dofs_property);
    const auto &ig_func_coeffs = ig_grid_function.get_coefficients();
    SafeSTLVector<Real> ig_func_elem_coeffs; // coefficients of the IgGridFunction restricted to the element
    for (const auto &global_dof : ig_basis_elem_global_dofs)
      ig_func_elem_coeffs.emplace_back(ig_func_coeffs[global_dof]);

    if (cache.template status_fill<_D<0>>())
    {
      using basis_element::_Value;
      auto &F = cache.template get_data<_D<0>>();
      F.fill(ig_basis_elem->template linear_combination<_Value,sdim>(ig_func_elem_coeffs,s_id_,dofs_property));
    }

    if (cache.template status_fill<_D<1>>())
    {
      using basis_element::_Gradient;
      auto &DF = cache.template get_data<_D<1>>();
      DF.fill(ig_basis_elem->template linear_combination<_Gradient,sdim>(ig_func_elem_coeffs,s_id_,dofs_property));
    }

    if (cache.template status_fill<_D<2>>())
    {
      using basis_element::_Hessian;
      auto &D2F = cache.template get_data<_D<2>>();
      D2F.fill(ig_basis_elem->template linear_combination<_Hessian,sdim>(ig_func_elem_coeffs,s_id_,dofs_property));
    }

//    Assert(cache.template status_fill<_D<3>>(),ExcNotImplemented());

  }

  cache.set_filled(true);
}


IGA_NAMESPACE_CLOSE

#include <igatools/functions/ig_grid_function_handler.inst>
