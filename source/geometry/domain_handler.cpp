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

#include <igatools/geometry/domain_handler.h>
#include <igatools/geometry/domain.h>
#include <igatools/geometry/domain_element.h>

IGA_NAMESPACE_OPEN

template<int dim_, int codim_>
DomainHandler<dim_, codim_>::
DomainHandler(std::shared_ptr<DomainType> domain)
  :
  domain_(domain),
  grid_func_handler_(domain->get_grid_function()->create_cache_handler()),
  flags_(Flags::none)
{
  Assert(domain_ != nullptr, ExcNullPtr());
}



template<int dim_, int codim_>
auto
DomainHandler<dim_, codim_>::
get_domain() const -> std::shared_ptr<DomainType>
{
  return domain_;
}

template<int dim_, int codim_>
auto
DomainHandler<dim_, codim_>::
get_element_cache(ElementAccessor &elem) const
-> typename ElementAccessor::CacheType &
{
  return  elem.local_cache_;
}


template<int dim_, int codim_>
void
DomainHandler<dim_, codim_>::
set_flags(const topology_variant &sdim,
          const Flags &flag)
{
  using GridFuncFlags = typename GridFuncType::Handler::Flags;
  using GridFlags = typename GridFuncType::GridType::Handler::Flags;

  GridFlags  grid_flag = GridFlags::none;
  GridFuncFlags  grid_func_flag = GridFuncFlags::none;
  Flags dom_flag = Flags::none;

  for (auto &fl : domain_element::all_flags)
    if (contains(flag, fl))
      dom_flag  |= domain_element::activate::domain[fl];


  for (auto &fl : domain_element::all_flags)
    if (contains(dom_flag, fl))
    {
      grid_func_flag |= domain_element::activate::grid_func[fl];
      grid_flag |= domain_element::activate::grid[fl];
    }

  grid_func_handler_->set_flags(sdim, grid_func_flag);
  grid_func_handler_->get_grid_handler().set_flags(sdim, grid_flag);
  auto disp = SetFlagsDispatcher(dom_flag, flags_);
  boost::apply_visitor(disp, sdim);
}

template<int dim_, int codim_>
template <int sdim>
void
DomainHandler<dim_, codim_>::
set_flags(const Flags &flag)
{
  this->set_flags(Topology<sdim>(), flag);
}



template<int dim_, int codim_>
void
DomainHandler<dim_, codim_>::
set_element_flags(const Flags &flag)
{
  this->set_flags(Topology<dim_>(), flag);
}

template<int dim_, int codim_>
void
DomainHandler<dim_, codim_>::
init_cache(ElementAccessor &elem,
           const eval_pts_variant &quad) const
{
  grid_func_handler_->init_cache(*(elem.grid_func_elem_), quad);

  auto disp = InitCacheDispatcher(*this, elem, flags_);
  boost::apply_visitor(disp, quad);
}

template<int dim_, int codim_>
void
DomainHandler<dim_, codim_>::
init_cache(ElementIterator &elem,
           const eval_pts_variant &quad) const
{
  this->init_cache(*elem, quad);
}

template<int dim_, int codim_>
void
DomainHandler<dim_, codim_>::
init_element_cache(ElementIterator &elem,
                   const std::shared_ptr<const Quadrature<dim_>> &quad) const
{
  this->init_cache(elem,quad);
}

template<int dim_, int codim_>
void
DomainHandler<dim_, codim_>::
init_element_cache(ElementAccessor &elem,
                   const std::shared_ptr<const Quadrature<dim_>> &quad) const
{
  this->init_cache(elem,quad);
}


template<int dim_, int codim_>
template <int sdim>
void
DomainHandler<dim_, codim_>::
fill_cache(ElementAccessor &elem,
           const int s_id) const
{
  this->fill_cache(Topology<sdim>(), elem, s_id);
}

template<int dim_, int codim_>
template <int sdim>
void
DomainHandler<dim_, codim_>::
fill_cache(ElementIterator &elem,
           const int s_id) const
{
  this->fill_cache(Topology<sdim>(), elem, s_id);
}


template<int dim_, int codim_>
void
DomainHandler<dim_, codim_>::
fill_cache(const topology_variant &sdim,
           ElementAccessor &elem,
           const int s_id) const
{
  grid_func_handler_->fill_cache(sdim, *(elem.grid_func_elem_), s_id);

  auto disp = FillCacheDispatcher(elem, s_id);
  boost::apply_visitor(disp, sdim);
}


template<int dim_, int codim_>
void
DomainHandler<dim_, codim_>::
fill_cache(const topology_variant &sdim,
           ElementIterator &elem,
           const int s_id) const
{
  this->fill_cache(sdim, *elem, s_id);
}



template<int dim_, int codim_>
void
DomainHandler<dim_, codim_>::
fill_element_cache(ElementAccessor &elem) const
{
  this->fill_cache(Topology<dim_>(), elem, 0);
}

template<int dim_, int codim_>
void
DomainHandler<dim_, codim_>::
fill_element_cache(ElementIterator &elem) const
{
  this->fill_cache(Topology<dim_>(), elem, 0);
}


template<int dim_, int codim_>
DomainHandler<dim_, codim_>::
SetFlagsDispatcher::
SetFlagsDispatcher(const Flags cache_flag, FlagsArray &flags)
  :
  cache_flag_(cache_flag),
  flags_(flags)
{}


template<int dim_, int codim_>
template<int sdim>
void
DomainHandler<dim_, codim_>::
SetFlagsDispatcher::
operator()(const Topology<sdim> &)
{
  flags_[sdim] |= cache_flag_;
}

template<int dim_, int codim_>
DomainHandler<dim_, codim_>::
InitCacheDispatcher::
InitCacheDispatcher(const self_t &domain_handler,
                    ElementAccessor &elem,
                    const FlagsArray &flags)
  :
  domain_handler_(domain_handler),
  elem_(elem),
  flags_(flags)
{}


template<int dim_, int codim_>
template<int sdim>
void
DomainHandler<dim_, codim_>::
InitCacheDispatcher::
operator()(const std::shared_ptr<const Quadrature<sdim>> &quad)
{
  auto &cache = domain_handler_.get_element_cache(elem_);

  const auto n_points = elem_.get_grid_function_element().get_grid_element().template get_quad<sdim>()
                        ->get_num_points();
  for (auto &s_id: UnitElement<dim_>::template elems_ids<sdim>())
  {
    auto &s_cache = cache.template get_sub_elem_cache<sdim>(s_id);
    s_cache.resize(flags_[sdim], n_points);
  }
}

template<int dim_, int codim_>
DomainHandler<dim_, codim_>::
FillCacheDispatcher::
FillCacheDispatcher(ElementAccessor &elem,
                    const int s_id)
  :
  elem_(elem),
  s_id_(s_id)
{}



template<int dim_, int codim_>
template<int sdim>
void
DomainHandler<dim_, codim_>::
FillCacheDispatcher::
operator()(const Topology<sdim> &)
{
  auto &cache = elem_.local_cache_.template get_sub_elem_cache<sdim>(s_id_);

  using _Measure = domain_element::_Measure;
  if (cache.template status_fill<_Measure>())
  {
    auto &s_elem = UnitElement<dim_>::template get_elem<sdim>(s_id_);

    const auto &DF = elem_.grid_func_elem_->template get_values_from_cache<grid_function_element::_D<1>, sdim>(s_id_);

    const auto n_points = DF.get_num_points();

    typename GridFunction<sdim, space_dim>::Gradient DF1;

    auto &measures = cache.template get_data<_Measure>();
    for (int pt = 0 ; pt < n_points; ++pt)
    {
      for (int l=0; l<sdim; ++l)
        DF1[l] = DF[pt][s_elem.active_directions[l]];

      measures[pt] = fabs(determinant<sdim,space_dim>(DF1));
    }
    measures.set_status_filled(true);
  }

  using _W_Measure = domain_element::_W_Measure;
  if (cache.template status_fill<_W_Measure>())
  {
    const auto &meas = cache.template get_data<_Measure>();
    const auto &w = elem_.grid_func_elem_->get_grid_element().template get_weights<sdim>(s_id_);

    auto &w_measures = cache.template get_data<_W_Measure>();

    auto it_w = w.cbegin();
    auto it_meas = meas.cbegin();
    for (auto &w_m : w_measures)
    {
      w_m = (*it_w) * (*it_meas);
      ++it_w;
      ++it_meas;
    }
    w_measures.set_status_filled(true);
  }

  using _InvJacobian = domain_element::_InvJacobian;
  if (cache.template status_fill<_InvJacobian>())
  {
    const auto &DF = elem_.grid_func_elem_->
                     template get_values_from_cache<grid_function_element::_D<1>,sdim>(s_id_);

    const auto n_points = DF.get_num_points();

    auto &D_invF = cache.template get_data<_InvJacobian>();
    Real det;
    for (int pt = 0 ; pt < n_points; ++pt)
      D_invF[pt] = inverse(DF[pt], det);

    D_invF.set_status_filled(true);
  }

  using _InvHessian = domain_element::_InvHessian;
  if (cache.template status_fill<_InvHessian>())
  {
    const auto &D2_F = elem_.grid_func_elem_->
                       template get_values_from_cache<grid_function_element::_D<2>,sdim>(s_id_);
    const auto &D1_invF = cache.template get_data<_InvJacobian>();
    auto &D2_invF       = cache.template get_data<_InvHessian>();

    const auto n_points = D2_F.get_num_points();

    for (int pt = 0 ; pt < n_points; ++pt)
    {
      for (int u=0; u<dim_; ++u)
      {
        const auto tmp_u = action(D2_F[pt], D1_invF[pt][u]);
        for (int v=0; v<dim_; ++v)
        {
          const auto tmp_u_v = action(tmp_u, D1_invF[pt][v]);
          D2_invF[pt][u][v] = - action(D1_invF[pt], tmp_u_v);
        } //end loop v
      } // end loop u
    } // end loop pt
    D2_invF.set_status_filled(true);
  }

  using _BoundaryNormal = domain_element::_BoundaryNormal;
  if (cache.template status_fill<_BoundaryNormal>())
  {
    Assert(dim_ == sdim+1, ExcMessage("The boundary normal is defined only if sdim == dim-1"));

    const auto n_hat  = UnitElement<dim_>::template get_elem<sdim>(s_id_).get_boundary_normal(0);

    const auto &D1_invF = cache.template get_data<_InvJacobian>();
    auto &bndry_normals = cache.template get_data<_BoundaryNormal>();

    int pt = 0;
    for (auto &bndry_normal_pt : bndry_normals)
    {
      const auto D1_invF_t = co_tensor(transpose(D1_invF[pt]));
      bndry_normal_pt = action(D1_invF_t, n_hat);
      bndry_normal_pt /= bndry_normal_pt.norm();
      ++pt;
    }

    bndry_normals.set_status_filled(true);
  }

  using _ExtNormal = domain_element::_ExtNormal;
  if (cache.template status_fill<_ExtNormal>())
  {
    Assert(sdim == dim_,  ExcMessage("The boundary normal is defined only if sdim == dim"));
    Assert(codim_ == 1, ExcNotImplemented());
    Assert(s_id_ == 0,ExcDimensionMismatch(s_id_,0));

    const auto &DF = elem_.grid_func_elem_->
                     template get_values_from_cache<grid_function_element::_D<1>,sdim>(s_id_);

    auto &ext_normals = cache.template get_data<_ExtNormal>();

    int pt = 0;
    for (auto &ext_normal_pt : ext_normals)
    {
      ext_normal_pt[0] = cross_product<dim_,codim_>(DF[pt]);
      ext_normal_pt[0] /= ext_normal_pt[0].norm();
      ++pt;
    }

    ext_normals.set_status_filled(true);
  }


  using  _FirstFundamentalForm = domain_element::_FirstFundamentalForm;
  if (cache.template status_fill<_FirstFundamentalForm>())
  {
    const auto &DF = elem_.grid_func_elem_->
                     template get_values_from_cache<grid_function_element::_D<1>,sdim>(s_id_);
    const auto n_points = DF.get_num_points();

    auto &first_fundamental_form = cache.template get_data<_FirstFundamentalForm>();

    for (int pt = 0; pt < n_points; ++pt)
    {
      const auto &A = DF[pt];
      const auto A_t   = co_tensor(transpose(A));
      first_fundamental_form[pt] = compose(A_t, A);
    }

    first_fundamental_form.set_status_filled(true);
  }


  using _SecondFundamentalForm = domain_element::_SecondFundamentalForm;
  if (cache.template status_fill<_SecondFundamentalForm>())
  {
    Assert(codim_==1, ExcNotImplemented());
    const auto &D2_F = elem_.grid_func_elem_->
                       template get_values_from_cache<grid_function_element::_D<2>,sdim>(s_id_);
    const auto n_points = D2_F.get_num_points();

    const auto &normal = cache.template get_data<_ExtNormal>();

    auto &second_fundamental_form = cache.template get_data<_SecondFundamentalForm>();

    for (int pt = 0; pt < n_points; ++pt)
    {
      auto &A = second_fundamental_form[pt];
      const auto &D2_F_pt = D2_F[pt];

      auto inner_normal = normal[pt][0];  //valid only for codim_ == 1
      inner_normal.operator-();

      for (int u = 0; u < dim_ ; ++u)
      {
        const auto B = co_tensor(transpose(D2_F_pt[u]));
        A[u] = action(B, inner_normal);
      }
    }
    second_fundamental_form.set_status_filled(true);
  }


  using _Curvature = domain_element::_Curvature;
  if (cache.template status_fill<_Curvature>())
  {
    Assert(sdim == dim_, ExcNotImplemented());
    Assert(codim_ == 1, ExcNotImplemented());

//        const auto H = elem.compute_second_fundamental_form();
//        const auto G_inv = elem.compute_inv_first_fundamental_form();

    const auto &form_1 = cache.template get_data<_FirstFundamentalForm>();
    const auto &form_2 = cache.template get_data<_SecondFundamentalForm>();

    const auto n_points = form_1.get_num_points();

    auto &curvatures = cache.template get_data<_Curvature>();

    Real det;
    for (int pt = 0; pt < n_points; ++pt)
    {
      //          const MetricTensor B = compose(H[pt], G_inv[pt]);
//          const auto B = compose(form_2[pt], inverse(form_1[pt],det));
      const auto A = unroll_to_matrix(
                       compose(form_2[pt], inverse(form_1[pt],det)));
      curvatures[pt] = A.eigen_values();
    }

    curvatures.set_status_filled(true);
  }

  using _ExtNormalD1 = domain_element::_ExtNormalD1;
  if (cache.template status_fill<_ExtNormalD1>())
  {
    Assert(sdim == dim_, ExcNotImplemented());
    Assert(codim_==1, ExcNotImplemented());

    const auto &DF = elem_.grid_func_elem_->
                     template get_values_from_cache<grid_function_element::_D<1>,sdim>(s_id_);

    const auto form_1 = cache.template get_data<_FirstFundamentalForm>();
    const auto form_2 = cache.template get_data<_SecondFundamentalForm>();

    const auto n_points = form_1.get_num_points();

    auto &Dn = cache.template get_data<_ExtNormalD1>();

    Real det;
    for (int pt = 0; pt< n_points; ++pt)
    {
//          const auto L = compose(DF[pt], inverse(form_1[pt],det));
      Dn[pt] = compose(
                 compose(DF[pt], inverse(form_1[pt],det)), form_2[pt]);
    }
    Dn.set_status_filled(true);
  }

  cache.set_filled(true);
} // end operator()


//    if (flag_.fill_inv_hessians())
//    {
//        const auto &D1_F = elem.get_gradients();
//        const auto &D2_F = elem.get_hessians();
//        const auto &D1_invF = std::get<1>(cache->inv_derivatives_);
//        auto &D2_invF = std::get<2>(cache->inv_derivatives_);
//
//        for (int i=0; i<n_points; ++i)
//            for (int u=0; u<dim_; ++u)
//            {
//                const auto tmp_u = action(D2_F[i], D1_invF[i][u]);
//                for (int v=0; v<dim_; ++v)
//                {
//                    const auto tmp_u_v = action(tmp_u, D1_invF[i][v]);
//                    D2_invF[i][u][v] = - action(D1_invF[i], tmp_u_v);
//                }
//            }
//    }
//}

IGA_NAMESPACE_CLOSE

#include <igatools/geometry/domain_handler.inst>

