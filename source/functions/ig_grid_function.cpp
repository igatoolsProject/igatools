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
IgGridFunction(const SharedPtrConstnessHandler<RefBasis> &ref_basis,
               const EpetraTools::Vector &coeff,
               const std::string &dofs_property,
               const std::string &name)
  :
  parent_t(
   ref_basis.data_is_const() ?
   SharedPtrConstnessHandler<GridType>(ref_basis->get_grid()) :
   SharedPtrConstnessHandler<GridType>(std::const_pointer_cast<Grid<dim>>(ref_basis->get_grid())),
   name),
  ref_basis_(ref_basis),
  dofs_property_(dofs_property)
{
  const auto &dof_distribution = *(ref_basis_->get_spline_space()->get_dof_distribution());
  const auto &active_dofs = dof_distribution.get_global_dofs(dofs_property_);

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
IgGridFunction(const SharedPtrConstnessHandler<RefBasis> &ref_basis,
               const IgCoefficients &coeffs,
               const std::string &dofs_property,
               const std::string &name)
  :
  parent_t(
   ref_basis.data_is_const() ?
   SharedPtrConstnessHandler<GridType>(ref_basis->get_grid()) :
   SharedPtrConstnessHandler<GridType>(std::const_pointer_cast<Grid<dim>>(ref_basis->get_grid())),
   name),
  ref_basis_(ref_basis),
  dofs_property_(dofs_property)
{
#ifndef NDEBUG
  const auto &dof_distribution = *(ref_basis_->get_spline_space()->get_dof_distribution());
  const auto &active_dofs = dof_distribution.get_global_dofs(dofs_property_);

  for (const auto glob_dof : active_dofs)
    coeffs_[glob_dof] = coeffs.at(glob_dof);
#else
  coeffs_ = coeffs;
#endif
}



#ifdef MESH_REFINEMENT

template<int dim,int space_dim>
void
IgGridFunction<dim,space_dim>::
rebuild_after_insert_knots(
  const SafeSTLArray<SafeSTLVector<Real>,dim> &knots_to_insert,
  const Grid<dim> &old_grid)
{
  auto ig_grid_function_pre_refinement =
    self_t::const_create(ref_basis_->get_basis_previous_refinement(),coeffs_);

  this->grid_function_previous_refinement_ = ig_grid_function_pre_refinement;


  const auto &ref_basis = *(this->get_basis());

  const int max_degree = ref_basis.get_spline_space()->get_max_degree();

  coeffs_ = space_tools::projection_l2_grid_function<dim,space_dim>(
              *ig_grid_function_pre_refinement,
              ref_basis,
              QGauss<dim>::create(max_degree+1),
              DofProperties::active);
}

#endif // MESH_REFINEMENT

template<int dim,int space_dim>
auto
IgGridFunction<dim,space_dim>::
create_cache_handler() const
-> std::unique_ptr<ElementHandler>
{
  using Handler = IgGridFunctionHandler<dim,space_dim>;
  return std::unique_ptr<Handler>(new Handler(
    std::dynamic_pointer_cast<const self_t>(this->shared_from_this())));
}


template<int dim,int space_dim>
auto
IgGridFunction<dim,space_dim>::
const_create(const std::shared_ptr<const RefBasis> &ref_basis,
             const IgCoefficients &coeffs,
             const std::string &dofs_property,
             const std::string &name) -> std::shared_ptr<const self_t>
{
  return std::shared_ptr<self_t>(new IgGridFunction(
    SharedPtrConstnessHandler<RefBasis>(ref_basis),coeffs,dofs_property,name));
}

template<int dim,int space_dim>
auto
IgGridFunction<dim,space_dim>::
create(const std::shared_ptr<RefBasis> &ref_basis,
       const IgCoefficients &coeffs,
       const std::string &dofs_property,
       const std::string &name) -> std::shared_ptr<self_t>
{
  auto func = std::shared_ptr<self_t>(new IgGridFunction(
    SharedPtrConstnessHandler<RefBasis>(ref_basis),coeffs,dofs_property,name));

#ifdef MESH_REFINEMENT
  func->create_connection_for_insert_knots(func);
#endif

  return func;
}


template<int dim,int space_dim>
auto
IgGridFunction<dim,space_dim>::
const_create(const std::shared_ptr<const RefBasis> &ref_basis,
             const EpetraTools::Vector &coeffs,
             const std::string &dofs_property,
             const std::string &name) -> std::shared_ptr<const self_t>
{
  return std::shared_ptr<self_t>(new IgGridFunction(
    SharedPtrConstnessHandler<RefBasis>(ref_basis),coeffs,dofs_property,name));
}

template<int dim,int space_dim>
auto
IgGridFunction<dim,space_dim>::
create(const std::shared_ptr<RefBasis> &ref_basis,
       const EpetraTools::Vector &coeffs,
       const std::string &dofs_property,
       const std::string &name) -> std::shared_ptr<self_t>
{
  return std::shared_ptr<self_t>(new IgGridFunction(
    SharedPtrConstnessHandler<RefBasis>(ref_basis),coeffs,dofs_property,name));
}



template<int dim,int space_dim>
auto
IgGridFunction<dim,space_dim>::
get_basis() const -> std::shared_ptr<const RefBasis>
{
  return ref_basis_.get_ptr_const_data();
}

template<int dim,int space_dim>
const IgCoefficients &
IgGridFunction<dim,space_dim>::
get_coefficients() const
{
  return coeffs_;
}

template<int dim,int space_dim>
const std::string &
IgGridFunction<dim,space_dim>::
get_dofs_property() const
{
  return dofs_property_;
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

  out.begin_item("ReferenceSpaceBasis<" +
                 to_string(dim) + ",1," +
                 to_string(space_dim) + ">:");
  ref_basis_->print_info(out);
  out.end_item();

  out.begin_item("IgCoefficients:");
  coeffs_.print_info(out);
  out.end_item();

  out << "Dofs property: " << dofs_property_ << std::endl;

  out.end_item();
}


IGA_NAMESPACE_CLOSE


#include <igatools/functions/ig_grid_function.inst>
