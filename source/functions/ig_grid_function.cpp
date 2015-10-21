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
IgGridFunction<dim,space_dim>::
IgGridFunction(const SharedPtrConstnessHandler<IgSpace> &space,
               const IgCoefficients &coeffs)
  :
  parent_t(
   space.data_is_const() ?
   SharedPtrConstnessHandler<GridType>(space.get_ptr_const_data()->get_ptr_const_grid()) :
   SharedPtrConstnessHandler<GridType>(space.get_ptr_data()->get_ptr_grid())),
  ig_space_(space)
{
#ifndef NDEBUG
  const auto &dof_distribution = *(ig_space_->get_ptr_const_dof_distribution());
  const auto &active_dofs = dof_distribution.get_dofs_id_same_property(DofProperties::active);

  for (const auto glob_dof : active_dofs)
    coeffs_[glob_dof] = coeffs.at(glob_dof);
#else
  coeffs_ = coeff;
#endif
}


template<int dim,int space_dim>
auto
IgGridFunction<dim,space_dim>::
create_cache_handler() const
-> std::unique_ptr<typename parent_t::ElementHandler>
{
  return std::make_unique<IgGridFunctionHandler<dim,space_dim>>(
    std::dynamic_pointer_cast<const self_t>(this->shared_from_this()));
}


template<int dim,int space_dim>
auto
IgGridFunction<dim,space_dim>::
const_create(const std::shared_ptr<const IgSpace> &space,
             const IgCoefficients &coeffs) -> std::shared_ptr<const self_t>
{
  return std::shared_ptr<const self_t>(new IgGridFunction(
    SharedPtrConstnessHandler<IgSpace>(space),coeffs));
}

template<int dim,int space_dim>
auto
IgGridFunction<dim,space_dim>::
create(const std::shared_ptr<IgSpace> &space,
       const IgCoefficients &coeffs) -> std::shared_ptr<self_t>
{
  return std::shared_ptr<self_t>(new IgGridFunction(
    SharedPtrConstnessHandler<IgSpace>(space),coeffs));
}



template<int dim,int space_dim>
auto
IgGridFunction<dim,space_dim>::
get_ig_space() const -> std::shared_ptr<const IgSpace>
{
  return ig_space_.get_ptr_const_data();
}

template<int dim,int space_dim>
const IgCoefficients &
IgGridFunction<dim,space_dim>::
get_coefficients() const
{
  return coeffs_;
}

template<int dim,int space_dim>
void
IgGridFunction<dim,space_dim>::
print_info(LogStream &out) const
{
  using std::to_string;
  out.begin_item("IgGridFunction<" +
                 to_string(dim) + "," +
                 to_string(space_dim) + ">");

  out.begin_item("ReferenceSpace<" +
                 to_string(dim) + ",1," +
                 to_string(space_dim) + ">:");
  ig_space_->print_info(out);
  out.end_item();

  out.begin_item("IgCoefficients:");
  coeffs_.print_info(out);
  out.end_item();

  out.end_item();
}


IGA_NAMESPACE_CLOSE


#include <igatools/functions/ig_grid_function.inst>
