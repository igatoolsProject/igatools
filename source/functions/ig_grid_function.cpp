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

#include <igatools/functions/ig_grid_function.h>
#include <igatools/functions/ig_grid_function_handler.h>
#include <igatools/basis_functions/basis_tools.h>
#include <igatools/base/quadrature_lib.h>


IGA_NAMESPACE_OPEN

#ifdef IGATOOLS_USES_TRILINOS

template<int dim,int range>
IgGridFunction<dim,range>::
IgGridFunction(const SharedPtrConstnessHandler<RefBasis> &ref_basis,
               const EpetraTools::Vector &coeff,
               const std::string &dofs_property)
  :
  parent_t(
   ref_basis.data_is_const() ?
   SharedPtrConstnessHandler<GridType>(ref_basis->get_grid()) :
   SharedPtrConstnessHandler<GridType>(std::const_pointer_cast<Grid<dim>>(ref_basis->get_grid()))),
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

#endif // IGATOOLS_USES_TRILINOS


template<int dim,int range>
IgGridFunction<dim,range>::
IgGridFunction(const SharedPtrConstnessHandler<RefBasis> &ref_basis,
               const IgCoefficients &coeffs,
               const std::string &dofs_property)
  :
  parent_t(
   ref_basis.data_is_const() ?
   SharedPtrConstnessHandler<GridType>(ref_basis->get_grid()) :
   SharedPtrConstnessHandler<GridType>(std::const_pointer_cast<Grid<dim>>(ref_basis->get_grid()))),
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



#ifdef IGATOOLS_WITH_MESH_REFINEMENT

template<int dim,int range>
void
IgGridFunction<dim,range>::
rebuild_after_insert_knots(
  const SafeSTLArray<SafeSTLVector<Real>,dim> &knots_to_insert,
  const Grid<dim> &old_grid)
{
  auto ig_grid_function_pre_refinement =
    self_t::const_create(ref_basis_->get_basis_previous_refinement(),coeffs_);

  this->grid_function_previous_refinement_ = ig_grid_function_pre_refinement;


  const auto &ref_basis = *(this->get_basis());

  const int max_degree = ref_basis.get_spline_space()->get_max_degree();

  coeffs_ = basis_tools::projection_l2_grid_function<dim,range>(
              *ig_grid_function_pre_refinement,
              ref_basis,
              QGauss<dim>::create(max_degree+1),
              DofProperties::active);
}

#endif // IGATOOLS_WITH_MESH_REFINEMENT

template<int dim,int range>
auto
IgGridFunction<dim,range>::
create_cache_handler() const
-> std::unique_ptr<Handler>
{
  using Handler = IgGridFunctionHandler<dim,range>;
  return std::unique_ptr<Handler>(new Handler(
    std::dynamic_pointer_cast<const self_t>(this->shared_from_this())));
}


template<int dim,int range>
auto
IgGridFunction<dim,range>::
const_create(const std::shared_ptr<const RefBasis> &ref_basis,
             const IgCoefficients &coeffs,
             const std::string &dofs_property) -> std::shared_ptr<const self_t>
{
  return std::shared_ptr<self_t>(new IgGridFunction(
    SharedPtrConstnessHandler<RefBasis>(ref_basis),coeffs,dofs_property));
}

template<int dim,int range>
auto
IgGridFunction<dim,range>::
create(const std::shared_ptr<RefBasis> &ref_basis,
       const IgCoefficients &coeffs,
       const std::string &dofs_property) -> std::shared_ptr<self_t>
{
  auto func = std::shared_ptr<self_t>(new IgGridFunction(
    SharedPtrConstnessHandler<RefBasis>(ref_basis),coeffs,dofs_property));

#ifdef IGATOOLS_WITH_MESH_REFINEMENT
  func->create_connection_for_insert_knots(func);
#endif

  return func;
}


#ifdef IGATOOLS_USES_TRILINOS

template<int dim,int range>
auto
IgGridFunction<dim,range>::
const_create(const std::shared_ptr<const RefBasis> &ref_basis,
             const EpetraTools::Vector &coeffs,
             const std::string &dofs_property) -> std::shared_ptr<const self_t>
{
  return std::shared_ptr<self_t>(new IgGridFunction(
    SharedPtrConstnessHandler<RefBasis>(ref_basis),coeffs,dofs_property));
}

template<int dim,int range>
auto
IgGridFunction<dim,range>::
create(const std::shared_ptr<RefBasis> &ref_basis,
       const EpetraTools::Vector &coeffs,
       const std::string &dofs_property) -> std::shared_ptr<self_t>
{
  return std::shared_ptr<self_t>(new IgGridFunction(
    SharedPtrConstnessHandler<RefBasis>(ref_basis),coeffs,dofs_property));
}

#endif // IGATOOLS_USES_TRILINOS


template<int dim,int range>
auto
IgGridFunction<dim,range>::
get_basis() const -> std::shared_ptr<const RefBasis>
{
  return ref_basis_.get_ptr_const_data();
}

template<int dim,int range>
const IgCoefficients &
IgGridFunction<dim,range>::
get_coefficients() const
{
  return coeffs_;
}

template<int dim,int range>
const std::string &
IgGridFunction<dim,range>::
get_dofs_property() const
{
  return dofs_property_;
}

template<int dim,int range>
void
IgGridFunction<dim,range>::
print_info(LogStream &out) const
{
  using std::to_string;
  out.begin_item("IgGridFunction<" +
                 to_string(dim) + "," +
                 to_string(range) + ">");

  out.begin_item("ReferenceBasis<" +
                 to_string(dim) + ",1," +
                 to_string(range) + ">:");
  ref_basis_->print_info(out);
  out.end_item();

  out.begin_item("IgCoefficients:");
  coeffs_.print_info(out);
  out.end_item();

  out << "Dofs property: " << dofs_property_ << std::endl;

  out << "Name: " << this->name_ << std::endl;

  out.end_item();
}


#ifdef IGATOOLS_WITH_SERIALIZATION
template<int dim,int range>
template<class Archive>
void
IgGridFunction<dim,range>::
serialize(Archive &ar)
{
  using std::to_string;
  const std::string base_name =
    "GridFunction_" + to_string(dim) + "_" +
    to_string(range);

  ar &make_nvp(base_name,base_class<parent_t>(this));

  ar &make_nvp("ref_basis_",ref_basis_);

  ar &make_nvp("coeffs_",coeffs_);

  ar &make_nvp("dofs_property_",dofs_property_);
}
#endif


IGA_NAMESPACE_CLOSE


#include <igatools/functions/ig_grid_function.inst>
