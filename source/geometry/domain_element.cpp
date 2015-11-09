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

#include <igatools/geometry/domain_element.h>
#include <igatools/functions/function_element.h>

IGA_NAMESPACE_OPEN

template<int dim_,int codim_>
DomainElement<dim_,codim_>::
DomainElement(const std::shared_ptr<ContainerType> &domain,
              const ListIt &index,
              const PropId &prop)
  :
  domain_(domain),
  grid_func_elem_(domain_->get_grid_function()->create_element(index,prop))
{}


template<int dim_,int codim_>
auto
DomainElement<dim_,codim_>::
get_index() const -> const IndexType &
{
  return grid_func_elem_->get_index();
}

template<int dim_,int codim_>
void
DomainElement<dim_,codim_>::
operator++()
{
  ++(*grid_func_elem_);
}

template<int dim_,int codim_>
void
DomainElement<dim_,codim_>::
move_to(const IndexType &elem_id)
{
  grid_func_elem_->move_to(elem_id);
}

template<int dim_,int codim_>
auto
DomainElement<dim_,codim_>::
get_grid_function_element() const -> const GridFuncElem &
{
  return *grid_func_elem_;
}

template<int dim_,int codim_>
auto
DomainElement<dim_,codim_>::
get_grid_function_element() -> GridFuncElem &
{
  return *grid_func_elem_;
}


template<int dim_,int codim_>
void
DomainElement<dim_,codim_>::
print_info(LogStream &out) const
{
  out.begin_item("GridFunctionElement:");
  grid_func_elem_->print_info(out);
  out.end_item();
}

template<int dim_,int codim_>
void
DomainElement<dim_,codim_>::
print_cache_info(LogStream &out) const
{
  out.begin_item("GridFunctionElement's cache");
  grid_func_elem_->print_cache_info(out);
  out.end_item();

  local_cache_.print_info(out);
}


template<int dim_,int codim_>
bool
DomainElement<dim_,codim_>::
operator ==(const self_t &elem) const
{
  Assert(domain_ == elem.domain_,
         ExcMessage("Cannot compare elements on different grid."));
  return (*grid_func_elem_ == *(elem.grid_func_elem_));
}



template<int dim_,int codim_>
bool
DomainElement<dim_,codim_>::
operator !=(const self_t &elem) const
{
  Assert(domain_ == elem.domain_,
         ExcMessage("Cannot compare elements on different grid."));
  return (*grid_func_elem_ != *(elem.grid_func_elem_));
}



template<int dim_,int codim_>
bool
DomainElement<dim_,codim_>::
operator <(const self_t &elem) const
{
  Assert(domain_ == elem.domain_,
         ExcMessage("Cannot compare elements on different grid."));
  return (*grid_func_elem_ < *(elem.grid_func_elem_));
}



template<int dim_,int codim_>
bool
DomainElement<dim_,codim_>::
operator >(const self_t &elem) const
{
  Assert(domain_ == elem.domain_,
         ExcMessage("Cannot compare elements on different grid."));
  return (*grid_func_elem_ > *(elem.grid_func_elem_));
}


#if 0
template<int dim_,int codim_>
template<int sdim>
auto
DomainElement<dim_,codim_>::
get_w_measures(const int s_id) const -> ValueVector<Real>
{
  const auto &meas = get_values_from_cache<_Measure, sdim>(s_id);
  const auto &w = grid_func_elem_->get_grid_element().template get_weights<sdim>(s_id);
  auto w_meas = meas;
  auto it_w = w.begin();
  for (auto &w_m : w_meas)
  {
    w_m *= *(it_w);
    ++it_w;
  }
  return w_meas;
}
#endif

template<int dim_,int codim_>
auto
DomainElement<dim_,codim_>::
get_element_measures() const -> const ValueVector<Real> &
{
  return get_values_from_cache<_Measure,dim_>(0);
}

template<int dim_,int codim_>
auto
DomainElement<dim_,codim_>::
get_element_w_measures() const -> const ValueVector<Real> &
{
  return get_values_from_cache<_W_Measure,dim_>(0);
}


template<int dim_,int codim_>
auto
DomainElement<dim_,codim_>::
get_element_points() const -> const ValueVector<Point> &
{
  return grid_func_elem_->template
  get_values_from_cache<grid_function_element::_D<0>,dim_>(0);
}

template<int dim_,int codim_>
auto
DomainElement<dim_,codim_>::
get_element_jacobians() const -> const ValueVector<Jacobian> &
{
  return grid_func_elem_->template
  get_values_from_cache<grid_function_element::_D<1>,dim_>(0);
}

template<int dim_,int codim_>
auto
DomainElement<dim_,codim_>::
get_element_hessians() const -> const ValueVector<Hessian> &
{
  return grid_func_elem_->template
  get_values_from_cache<grid_function_element::_D<2>,dim_>(0);
}

template<int dim_,int codim_>
auto
DomainElement<dim_,codim_>::
get_exterior_normals() const -> const ValueVector<SafeSTLArray<Point, codim_> > &
{
  const int sdim = dim_;
  const int s_id = 0;
  AssertThrow(codim_ == 1, ExcNotImplemented());
#if 0
  ValueVector<SafeSTLArray<Point, codim_>> res;

  const auto &DF = grid_func_elem_->template
  get_values<grid_function_element::_D<1>, sdim>(s_id);
  const auto n_points = DF.get_num_points();
  res.resize(n_points);

  for (int pt = 0; pt < n_points; ++pt)
  {
    res[pt][0] = cross_product<dim_, codim_>(DF[pt]);
    res[pt][0] /= res[pt][0].norm();
  }

  return res;
#endif
  return get_values_from_cache<domain_element::_ExtNormal,sdim>(s_id);
}



#if 0
template<int dim_,int codim_>
auto
DomainElement<dim_,codim_>::
compute_inv_first_fundamental_form() const -> ValueVector<MetricTensor>
{
  ValueVector<MetricTensor> res;
  const auto &DF = this->template get_values<_Gradient, dim>(0);
  const auto n_points = DF.get_num_points();

  res.resize(n_points);
  Real det;
  for (int i = 0; i< n_points; ++i)
  {
    const auto &A = DF[i];
    const auto A_t   = co_tensor(transpose(A));
    const auto G     = compose(A_t, A);
    res[i] = inverse(G, det);
  }

  return res;
}



template<int dim_,int codim_>
auto
DomainElement<dim_,codim_>::
compute_second_fundamental_form() const -> ValueVector<MetricTensor>
{
  Assert(codim==1, ExcNotImplemented());

  const auto &D2_F  = this->template get_values<_Hessian, dim>(0);
  const auto normal = this->get_external_normals();

  const auto n_points = D2_F.get_num_points();

  // const auto G_inv = compute_inv_first_fundamental_form();

  ValueVector<MetricTensor> res;
  res.resize(n_points);

  MetricTensor A;
  for (int pt = 0; pt < n_points; ++pt)
  {
    const auto &D2_F_pt = D2_F[pt];
    for (int u=0; u<dim; ++u)
    {
      const auto B = co_tensor(transpose(D2_F_pt[u]));
      A[u] = action(B, normal[pt]);
    }
    res[pt] = -A;

    // res[pt] = -compose(A, G_inv[pt]);

  }

  return res;
}


template<int dim_,int codim_>
auto
DomainElement<dim_,codim_>::
get_principal_curvatures() const -> const ValueVector<SafeSTLVector<Real>> &
{
#if 0
  Assert(codim==1, ExcNotImplemented());

  const auto H = compute_second_fundamental_form();
  const auto G_inv = compute_inv_first_fundamental_form();

  const auto n_points = H.get_num_points();

  ValueVector<SafeSTLVector<Real>> res(n_points);

  for (int pt = 0; pt < n_points; ++pt)
  {
    const MetricTensor B = compose(H[pt], G_inv[pt]);
    const auto A = unroll_to_matrix(B);
    res[pt] = A.eigen_values();
  }
  return res;
#endif
  Assert(codim==1, ExcNotImplemented());
  return get_values_from_cache<_Curvature,dim_>(0);
}







template<int dim_,int codim_>
auto
DomainElement<dim_,codim_>::
get_D_external_normals() const -> ValueVector< Derivative<1> >
{
  Assert(codim==1, ExcNotImplemented());

  const auto H = compute_second_fundamental_form();
  const auto &DF = this->template get_values<_Gradient, dim>(0);
  const auto G_inv = compute_inv_first_fundamental_form();

  const auto n_points = H.get_num_points();
  ValueVector< Derivative<1> > Dn(n_points);

  for (int pt = 0; pt< n_points; ++pt)
  {
    auto L = compose(DF[pt], G_inv[pt]);
    Dn[pt] = compose(L, H[pt]);
  }
  return Dn;
}
#endif


//template<int dim_,int codim_>
//auto
//DomainElement<dim_,codim_>::
//clone() const -> std::shared_ptr<self_t>
//{
//    auto elem = std::make_shared<self_t>(*this,CopyPolicy::deep);
//    Assert(elem != nullptr, ExcNullPtr());
//    return elem;
//}



//  const auto &D2_F  = this->template get_values<2, dim>(0);
//  const auto normal = this->get_external_normals();
//
//  const auto n_points = D2_F.get_num_points();
//  ValueVector<SafeSTLVector<Real>> res(n_points);
//  const auto G_inv = compute_inv_first_fundamental_form();
//  DenseMatrix A(dim, dim);
//  for (int pt = 0; pt < n_points; ++pt)
//  {
//      const auto B = unroll_to_matrix(G_inv[pt]);
//      for (int i = 0; i<dim; ++i)
//          for (int j = 0; j<dim; ++j)
//              A(i,j) = -scalar_product(D2_F[pt][i][j], normal[pt]);
//      DenseMatrix C(dim, dim);
//      boost::numeric::ublas::axpy_prod(B, A, C, true);
//
//      res[pt] = C.eigen_values();
//  }



//    template<int sub_dim>
//    ValueVector<Points<space_dim> >
//    get_boundary_normals(const int s_id) const
//    {
//        Assert(dim==sub_dim+1, ExcNotImplemented());
//        ValueVector<Points<space_dim>> res;
//        const auto &DF_inv = get_inverse_values<1, sub_dim>(s_id);
//        const auto n_hat  = get_grid()->template get_boundary_normals<sub_dim>(s_id)[0];
//
//        const auto n_points = DF_inv.get_num_points();
//        res.resize(n_points);
//        for (int i = 0; i< n_points; ++i)
//        {
//            const auto DF_inv_t = co_tensor(transpose(DF_inv[i]));
//            res[i] = action(DF_inv_t, n_hat);
//            res[i] /= res[i].norm();
//        }
//        return res;
//    }




IGA_NAMESPACE_CLOSE

#include <igatools/geometry/domain_element.inst>

