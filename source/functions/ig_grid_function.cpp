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
#include <igatools/basis_functions/space_tools.h>
#include <igatools/base/quadrature_lib.h>


IGA_NAMESPACE_OPEN

template<int dim,int space_dim>
IgGridFunction<dim,space_dim>::
IgGridFunction(const SharedPtrConstnessHandler<IgSpace> &space,
               const EpetraTools::Vector &coeff)
  :
  parent_t(
   space.data_is_const() ?
   SharedPtrConstnessHandler<GridType>(space.get_ptr_const_data()->get_ptr_const_grid()) :
   SharedPtrConstnessHandler<GridType>(space.get_ptr_data()->get_ptr_grid())),
  ig_space_(space)
{
  const auto &dof_distribution = *(ig_space_->get_ptr_const_dof_distribution());
  const auto &active_dofs = dof_distribution.get_dofs_id_same_property(DofProperties::active);

  const auto &epetra_map = coeff.Map();

  for (const auto glob_dof : active_dofs)
  {
    auto loc_id = epetra_map.LID(glob_dof);
    Assert(loc_id >= 0,
           ExcMessage("Global dof " + std::to_string(glob_dof) + " not present in the input EpetraTools::Vector."));
    coeffs_[glob_dof] = coeff[loc_id];
  }
}


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




#ifdef MESH_REFINEMENT
template<int dim,int space_dim>
void
IgGridFunction<dim,space_dim>::
create_connection_for_insert_knots(const std::shared_ptr<self_t> &ig_grid_function)
{
  Assert(ig_grid_function != nullptr, ExcNullPtr());
  Assert(&(*ig_grid_function) == &(*this), ExcMessage("Different objects."));

  auto func_to_connect =
    std::bind(&self_t::rebuild_after_insert_knots,
              ig_grid_function.get(),
              std::placeholders::_1,
              std::placeholders::_2);

  using SlotType = typename Grid<dim>::SignalInsertKnotsSlot;
  this->get_grid()
  ->connect_insert_knots(SlotType(func_to_connect).track_foreign(ig_grid_function));
}


template<int dim,int space_dim>
void
IgGridFunction<dim,space_dim>::
rebuild_after_insert_knots(
  const SafeSTLArray<SafeSTLVector<Real>,dim> &knots_to_insert,
  const Grid<dim> &old_grid)
{
  auto ig_grid_function_pre_refinement =
    self_t::const_create(ig_space_->get_space_previous_refinement(),coeffs_);

  this->grid_function_previous_refinement_ = ig_grid_function_pre_refinement;


  const auto &ig_space = *(this->get_ig_space());

  const int max_degree = ig_space.get_max_degree();

  coeffs_ = space_tools::projection_l2_ig_grid_function<dim,space_dim>(
              *ig_grid_function_pre_refinement,
              ig_space,
              QGauss<dim>::create(max_degree+1),
              DofProperties::active);
//*/

//  Assert(false,ExcNotImplemented());
  /*
  this->grid_function_previous_refinement_ =
      self_t::const_create(
        this->grid_func_->get_grid_function_previous_refinement());
        //*/
}
#endif // MESH_REFINEMENT

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
  auto func = std::shared_ptr<self_t>(new IgGridFunction(
    SharedPtrConstnessHandler<IgSpace>(space),coeffs));

#ifdef MESH_REFINEMENT
  func->create_connection_for_insert_knots(func);
#endif

  return func;
}


template<int dim,int space_dim>
auto
IgGridFunction<dim,space_dim>::
const_create(const std::shared_ptr<const IgSpace> &space,
             const EpetraTools::Vector &coeffs) -> std::shared_ptr<const self_t>
{
  return std::shared_ptr<const self_t>(new IgGridFunction(
    SharedPtrConstnessHandler<IgSpace>(space),coeffs));
}

template<int dim,int space_dim>
auto
IgGridFunction<dim,space_dim>::
create(const std::shared_ptr<IgSpace> &space,
       const EpetraTools::Vector &coeffs) -> std::shared_ptr<self_t>
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
