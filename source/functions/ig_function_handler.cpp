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

#include <igatools/functions/ig_function_handler.h>
#include <igatools/functions/ig_function.h>
#include <igatools/functions/function_element.h>
#include <igatools/basis_functions/physical_basis_element_handler.h>

IGA_NAMESPACE_OPEN


template<int dim,int codim,int range,int rank>
IgFunctionHandler<dim,codim,range,rank>::
IgFunctionHandler(const std::shared_ptr<FunctionType> &ig_function)
  :
  parent_t(ig_function),
  ig_basis_handler_(ig_function->get_basis()->create_cache_handler()),
  ig_function_(ig_function)
{}


#if 0
template<int dim_, int space_dim_>
auto
IgGridFunctionHandler<dim_, space_dim_>::
get_ig_grid_function() const -> std::shared_ptr<GridFunctionType>
{
  return ig_grid_function_;
}
#endif



template<int dim,int codim,int range,int rank>
auto
IgFunctionHandler<dim,codim,range,rank>::
set_flags(const topology_variant &sdim,
          const Flags &flag) -> void
{
//  this->get_grid_handler().set_flags(sdim,flag);

  parent_t::set_flags(sdim,flag);


  using SpFlags = space_element::Flags;
  SpFlags ig_space_elem_flags = SpFlags::none;
  if (contains(flag,Flags::D0))
    ig_space_elem_flags |= SpFlags::value;
  if (contains(flag,Flags::D1))
    ig_space_elem_flags |= SpFlags::gradient;
  if (contains(flag,Flags::D2))
    ig_space_elem_flags |= SpFlags::hessian;

  ig_basis_handler_->set_flags_impl(sdim,ig_space_elem_flags);
  //*/
}

#if 0
template<int dim_, int space_dim_>
void
IgGridFunctionHandler<dim_, space_dim_>::
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


template<int dim,int codim,int range,int rank>
auto
IgFunctionHandler<dim,codim,range,rank>::
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


template<int dim,int codim,int range,int rank>
template<int sdim>
void
IgFunctionHandler<dim,codim,range,rank>::
FillCacheDispatcher::
operator()(const Topology<sdim> &sub_elem)
{
  static_assert(sdim >= 0 && sdim <= dim,
                "Non-admissible sub-element dimensionality.");

  auto &domain_elem = ig_function_elem_.get_domain_element();
  ig_function_handler_.get_domain_handler()->template fill_cache<sdim>(domain_elem,s_id_);

  using IgFunc = IgFunction<dim,codim,range,rank>;
  const auto &ig_function =
    *(std::dynamic_pointer_cast<const IgFunc>(ig_function_handler_.get_function()));

  auto &local_cache = ig_function_handler_.get_element_cache(ig_function_elem_);
  auto &cache = local_cache.template get_sub_elem_cache<sdim>(s_id_);


  if (!cache.fill_none())
  {
    const auto &grid_elem = domain_elem.get_grid_function_element().get_grid_element();
    const auto &grid_elem_id = grid_elem.get_index();

    const auto ig_space = ig_function.get_basis();
    const auto &ig_basis_handler = *ig_function_handler_.ig_basis_handler_;
    auto ig_space_elem = ig_space->begin();
    ig_space_elem->move_to(grid_elem_id);

    ig_basis_handler.template init_cache<sdim>(*ig_space_elem,grid_elem.template get_quad<sdim>());
    ig_basis_handler.template fill_cache<sdim>(*ig_space_elem,s_id_);


    const auto &dofs_property = ig_function.get_dofs_property();

    const auto &ig_space_elem_global_dofs = ig_space_elem->get_local_to_global(dofs_property);
    const auto &ig_func_coeffs = ig_function.get_coefficients();
    SafeSTLVector<Real> ig_func_elem_coeffs; // coefficients of the IgGridFunction restricted to the element
    for (const auto &global_dof : ig_space_elem_global_dofs)
      ig_func_elem_coeffs.emplace_back(ig_func_coeffs[global_dof]);

    using _D0 = function_element::template _D<0>;
    if (cache.template status_fill<_D0>())
    {
      using space_element::_Value;
      auto &F = cache.template get_data<_D0>();
      F.fill(ig_space_elem->template linear_combination<_Value,sdim>(ig_func_elem_coeffs,s_id_,dofs_property));
    }

    using _D1 = function_element::template _D<1>;
    if (cache.template status_fill<_D1>())
    {
      using space_element::_Gradient;
      auto &DF = cache.template get_data<_D1>();
      DF.fill(ig_space_elem->template linear_combination<_Gradient,sdim>(ig_func_elem_coeffs,s_id_,dofs_property));
    }

    using _D2 = function_element::template _D<2>;
    if (cache.template status_fill<_D2>())
    {
      using space_element::_Hessian;
      auto &D2F = cache.template get_data<_D2>();
      D2F.fill(ig_space_elem->template linear_combination<_Hessian,sdim>(ig_func_elem_coeffs,s_id_,dofs_property));
    }

//    Assert(cache.template status_fill<_D<3>>(),ExcNotImplemented());

  }

  cache.set_filled(true);
}

IGA_NAMESPACE_CLOSE

#include <igatools/functions/ig_function_handler.inst>
