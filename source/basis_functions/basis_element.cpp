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



#include <igatools/basis_functions/basis_element.h>
#include <igatools/basis_functions/reference_basis_element.h>
#include <igatools/basis_functions/physical_basis_element.h>


IGA_NAMESPACE_OPEN


template<int dim_,int codim_,int range_,int rank_>
BasisElement<dim_,codim_,range_,rank_>::
BasisElement(const std::shared_ptr<Sp> &basis)
  :
  basis_(basis)
{}



template<int dim_,int codim_,int range_,int rank_>
void
BasisElement<dim_,codim_,range_,rank_>::
print_info(LogStream &out) const
{
  out.begin_item("Element global connectivity (property=\"" + DofProperties::active + "\"):");
  const auto glob_dofs = this->get_local_to_global(DofProperties::active);
  glob_dofs.print_info(out);
  out.end_item();
}


template<int dim_,int codim_,int range_,int rank_>
void
BasisElement<dim_,codim_,range_,rank_>::
print_cache_info(LogStream &out) const
{
  out.begin_item("GridElement<" + std::to_string(dim) + "> cache:");
  this->get_grid_element().print_cache_info(out);
  out.end_item();

  //    Assert(all_sub_elems_cache_ != nullptr,ExcNullPtr());
  all_sub_elems_cache_.print_info(out);
}



template<int dim_,int codim_,int range_,int rank_>
auto
BasisElement<dim_,codim_,range_,rank_>::
get_index() const -> const IndexType &
{
  return this->get_grid_element().get_index();
}




template<int dim_,int codim_,int range_,int rank_>
SafeSTLVector<Index>
BasisElement<dim_,codim_,range_,rank_>::
get_local_to_global(const std::string &dofs_property) const
{
  SafeSTLVector<Index> dofs_global;
  SafeSTLVector<Index> dofs_loc_to_patch;
  SafeSTLVector<Index> dofs_loc_to_elem;
  this->basis_->get_spline_space()->get_element_dofs(
    this->get_index(),
    dofs_global,
    dofs_loc_to_patch,
    dofs_loc_to_elem,
    dofs_property);

  return dofs_global;
}

template<int dim_,int codim_,int range_,int rank_>
SafeSTLVector<Index>
BasisElement<dim_,codim_,range_,rank_>::
get_local_to_patch(const std::string &dofs_property) const
{
  SafeSTLVector<Index> dofs_global;
  SafeSTLVector<Index> dofs_loc_to_patch;
  SafeSTLVector<Index> dofs_loc_to_elem;
  this->basis_->get_spline_space()->get_element_dofs(
    this->get_index(),
    dofs_global,
    dofs_loc_to_patch,
    dofs_loc_to_elem,
    dofs_property);

  return dofs_loc_to_patch;
}

template<int dim_,int codim_,int range_,int rank_>
SafeSTLVector<Index>
BasisElement<dim_,codim_,range_,rank_>::
get_local_dofs(const std::string &dofs_property) const
{
  SafeSTLVector<Index> dofs_global;
  SafeSTLVector<Index> dofs_loc_to_patch;
  SafeSTLVector<Index> dofs_loc_to_elem;
  this->basis_->get_spline_space()->get_element_dofs(
    this->get_index(),
    dofs_global,
    dofs_loc_to_patch,
    dofs_loc_to_elem,
    dofs_property);

  return dofs_loc_to_elem;
}

template<int dim_,int codim_,int range_,int rank_>
Size
BasisElement<dim_,codim_,range_,rank_>::
get_num_basis(const std::string &dofs_property) const
{
  const auto dofs_global = this->get_local_to_global(dofs_property);
  return dofs_global.size();
}


template<int dim_,int codim_,int range_,int rank_>
bool
BasisElement<dim_,codim_,range_,rank_>::
operator==(const self_t &a) const
{
  Assert(this->has_same_basis_of(a),
         ExcMessage("Comparison between elements defined on different spaces"));
  return this->get_grid_element() == a.get_grid_element();
}

template<int dim_,int codim_,int range_,int rank_>
bool
BasisElement<dim_,codim_,range_,rank_>::
operator!=(const self_t &a) const
{
  Assert(this->has_same_basis_of(a),
         ExcMessage("Comparison between elements defined on different spaces"));
  return this->get_grid_element() != a.get_grid_element();
}



template<int dim_,int codim_,int range_,int rank_>
template <int k>
ValueVector<Real>
BasisElement<dim_,codim_,range_,rank_>::
get_w_measures(const int j) const
{
  ValueVector<Real> w_measures;

  using RefElem = const ReferenceBasisElement<dim_,range_,rank_>;
  RefElem *as_ref_elem = dynamic_cast<RefElem *>(this);
  if (as_ref_elem)
    w_measures = as_ref_elem->template get_w_measures<k>(j);

  using PhysElem = const PhysicalBasisElement<dim_,range_,rank_,codim_>;
  PhysElem *as_phys_elem = dynamic_cast<PhysElem *>(this);
  if (as_phys_elem)
    w_measures = as_phys_elem->template get_w_measures<k>(j);

  return w_measures;
}


template<int dim_,int codim_,int range_,int rank_>
ValueVector<Real>
BasisElement<dim_,codim_,range_,rank_>::
get_element_w_measures() const
{
  return this->template get_w_measures<dim>(0);
}

template<int dim_,int codim_,int range_,int rank_>
auto
BasisElement<dim_,codim_,range_,rank_>::
get_basis() const -> std::shared_ptr<Sp>
{
  return basis_;
}


template<int dim_,int codim_,int range_,int rank_>
bool
BasisElement<dim_,codim_,range_,rank_>::
has_same_basis_of(const self_t &elem) const
{
  return (basis_ == elem.basis_);
}


template<int dim_,int codim_,int range_,int rank_>
auto
BasisElement<dim_,codim_,range_,rank_>::
get_element_values(const std::string &dofs_property) const
-> ValueTable<Value>
{
  return this->template get_basis_data<space_element::_Value,dim_>(0,dofs_property);
}

template<int dim_,int codim_,int range_,int rank_>
auto
BasisElement<dim_,codim_,range_,rank_>::
get_element_gradients(const std::string &dofs_property) const
-> ValueTable<Derivative<1>>
{
  return this->template get_basis_data<space_element::_Gradient,dim_>(0,dofs_property);
}

template<int dim_,int codim_,int range_,int rank_>
auto
BasisElement<dim_,codim_,range_,rank_>::
get_element_hessians(const std::string &dofs_property) const
-> ValueTable<Derivative<2>>
{
  return this->template get_basis_data<space_element::_Hessian,dim_>(0,dofs_property);
}

template<int dim_,int codim_,int range_,int rank_>
auto
BasisElement<dim_,codim_,range_,rank_>::
get_element_divergences(const std::string &dofs_property) const
-> ValueTable<Div>
{
  return this->template get_basis_data<space_element::_Divergence,dim_>(0,dofs_property);
}

template<int dim_,int codim_,int range_,int rank_>
DenseMatrix
BasisElement<dim_,codim_,range_,rank_>::
integrate_u_v(const PropId &dofs_property)
{
  const auto &w_meas = this->get_element_w_measures();
  const auto &u = this->get_element_values(dofs_property);

  const int n_pts = u.get_num_points();

  const int n_basis = u.get_num_functions();
  DenseMatrix M(n_basis,n_basis);
  M = 0.0;

  for (int i = 0; i < n_basis; ++i)
  {
    const auto u_i = u.get_function_view(i);

    for (int j = i; j < n_basis; ++j)
    {
      const auto u_j = u.get_function_view(j);

      for (int pt = 0; pt < n_pts; ++pt)
        M(i,j) += scalar_product(u_i[pt],u_j[pt]) * w_meas[pt];
    } // end loop j

    for (int j = 0; j < i; ++j)
      M(i,j) = M(j,i);
  } // end loop i

  return M;
}

template<int dim_,int codim_,int range_,int rank_>
DenseMatrix
BasisElement<dim_,codim_,range_,rank_>::
integrate_gradu_gradv(const PropId &dofs_property)
{
  const auto &w_meas = this->get_element_w_measures();
  const auto &gradu = this->get_element_gradients(dofs_property);

  const int n_pts = gradu.get_num_points();

  const int n_basis = gradu.get_num_functions();
  DenseMatrix M(n_basis,n_basis);
  M = 0.0;

  for (int i = 0; i < n_basis; ++i)
  {
    const auto gradu_i = gradu.get_function_view(i);

    for (int j = i; j < n_basis; ++j)
    {
      const auto gradu_j = gradu.get_function_view(j);

      for (int pt = 0; pt < n_pts; ++pt)
        M(i,j) += scalar_product(gradu_i[pt],gradu_j[pt]) * w_meas[pt];
    } // end loop j

    for (int j = 0; j < i; ++j)
      M(i,j) = M(j,i);
  } // end loop i

  return M;
}

template<int dim_,int codim_,int range_,int rank_>
DenseVector
BasisElement<dim_,codim_,range_,rank_>::
integrate_u_func(const ValueVector<Value> &func_at_points,
                 const PropId &dofs_property)
{
  const auto &u = this->get_element_values(dofs_property);

  const int n_pts = u.get_num_points();
  Assert(n_pts == func_at_points.get_num_points(),
         ExcDimensionMismatch(n_pts,func_at_points.get_num_points()));

  const int n_basis = u.get_num_functions();
  DenseVector R(n_basis);
  R = 0.0;

  const auto &w_meas = this->get_element_w_measures();

  for (int i = 0; i < n_basis; ++i)
  {
    const auto u_i = u.get_function_view(i);

    for (int pt = 0; pt < n_pts; ++pt)
      R(i) += scalar_product(u_i[pt],func_at_points[pt]) * w_meas[pt];
  } // end loop i

  return R;
}


IGA_NAMESPACE_CLOSE


#include <igatools/basis_functions/basis_element.inst>


