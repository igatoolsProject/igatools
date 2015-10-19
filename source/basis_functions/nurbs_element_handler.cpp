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

#include <igatools/basis_functions/nurbs_element_handler.h>
#include <igatools/basis_functions/nurbs_element.h>

#include <algorithm>

#ifdef NURBS
using std::shared_ptr;
IGA_NAMESPACE_OPEN

template<int dim_, int range_ , int rank_>
NURBSElementHandler<dim_, range_, rank_>::
NURBSElementHandler(shared_ptr<const Space> space)
  :
  base_t(space),
  bsp_elem_handler_(space->get_spline_space()->create_cache_handler()),
  w_func_elem_handler_(
   dynamic_cast<IgGridFunctionHandler<dim_,1> *>(
     space->get_weight_func()->create_cache_handler().release()))
{
  Assert(w_func_elem_handler_ != nullptr, ExcNullPtr());
}

template<int dim_, int range_ , int rank_>
void
NURBSElementHandler<dim_, range_, rank_>::
set_flags_impl(const topology_variant &sdim,
               const typename space_element::Flags &flag)
{
  using BasisFlags = space_element::Flags;
  auto nrb_elem_flags = flag;
  if (contains(flag,space_element::Flags::divergence))
    nrb_elem_flags |= space_element::Flags::gradient;


  auto bsp_elem_flags = BasisFlags::none;

  using WFuncFlags = grid_function_element::Flags;
  auto w_func_elem_flags = WFuncFlags::none;

  if (contains(nrb_elem_flags, BasisFlags::value))
  {
    bsp_elem_flags |= BasisFlags::value;
    w_func_elem_flags |= WFuncFlags::D0;
  }
  if (contains(nrb_elem_flags, BasisFlags::gradient))
  {
    bsp_elem_flags |= BasisFlags::value | BasisFlags::gradient ;
    w_func_elem_flags |= WFuncFlags::D0 | WFuncFlags::D1;
  }
  if (contains(nrb_elem_flags, BasisFlags::hessian))
  {
    bsp_elem_flags |= BasisFlags::value | BasisFlags::gradient | BasisFlags::hessian;
    w_func_elem_flags |= WFuncFlags::D0 | WFuncFlags::D1 | WFuncFlags::D2;
  }


  bsp_elem_handler_->set_flags_impl(sdim,bsp_elem_flags);
  w_func_elem_handler_->set_flags(sdim,w_func_elem_flags);

  auto set_flag_dispatcher = SetFlagDispatcher(nrb_elem_flags,*this);
  boost::apply_visitor(set_flag_dispatcher,sdim);
}


template<int dim_, int range_ , int rank_>
void
NURBSElementHandler<dim_, range_, rank_>::
init_cache_impl(BaseElem &elem,
                const eval_pts_variant &quad) const
{
  auto init_cache_dispatcher = InitCacheDispatcher(*this,elem);
  boost::apply_visitor(init_cache_dispatcher,quad);
}


template<int dim_, int range_ , int rank_>
void
NURBSElementHandler<dim_, range_, rank_>::
fill_cache_impl(const topology_variant &sdim,
                BaseElem &elem,
                const int s_id) const
{
  auto fill_cache_dispatcher = FillCacheDispatcher(*this,elem,s_id);
  boost::apply_visitor(fill_cache_dispatcher,sdim);
}







template<int dim_, int range_ , int rank_>
auto
NURBSElementHandler<dim_, range_, rank_>::
get_nurbs_space() const -> std::shared_ptr<const Space>
{
  auto nrb_space = std::dynamic_pointer_cast<const Space>(this->get_space());
  Assert(nrb_space != nullptr,ExcNullPtr());
  return nrb_space;
}

template<int dim_, int range_ , int rank_>
void
NURBSElementHandler<dim_, range_, rank_>::
print_info(LogStream &out) const
{
  Assert(false,ExcNotImplemented());
  AssertThrow(false,ExcNotImplemented());
#if 0

  out.begin_item("Grid Cache:");
  base_t::print_info(out);
  out.end_item();

  cacheutils::print_caches(splines1d_, out);
  //out.begin_item("One dimensional splines cache:");
  //splines1d_.print_info(out);
  //out.end_item();
#endif
}


template<int dim_, int range_ , int rank_>
NURBSElementHandler<dim_, range_, rank_>::
SetFlagDispatcher::
SetFlagDispatcher(const typename space_element::Flags nrb_flag,
                  self_t &nrb_handler)
  :
  nrb_flag_(nrb_flag),
  nrb_handler_(nrb_handler)
{}


template<int dim_, int range_ , int rank_>
NURBSElementHandler<dim_, range_, rank_>::
InitCacheDispatcher::
InitCacheDispatcher(const self_t &nrb_handler,
                    SpaceElement<dim_,0,range_,rank_> &elem)
  :
  nrb_handler_(nrb_handler),
  elem_(elem)
{}


template<int dim_, int range_ , int rank_>
NURBSElementHandler<dim_, range_, rank_>::
FillCacheDispatcher::
FillCacheDispatcher(const self_t &nrb_handler,
                    SpaceElement<dim_,0,range_,rank_> &elem,
                    const int s_id)
  :
  nrb_handler_(nrb_handler),
  elem_(elem),
  s_id_(s_id)
{}

template<int dim_, int range_ , int rank_>
void
NURBSElementHandler<dim_, range_, rank_>::
FillCacheDispatcher::
evaluate_nurbs_values_from_bspline(
  const BSplineElem &bspline_elem,
  const SafeSTLVector<Index> &bsp_local_to_patch,
  const ValueTable<typename BSplineElem::Value> &P,
  const IgCoefficients &w_coefs,
  const ValueVector<typename WeightFuncElem::Value> &Q,
  DataWithFlagStatus<ValueTable<Value>> &phi) const
{
  /*
   * This function evaluates the values of the n+1 NURBS basis function R_0,...,R_n
   * from the set of BSpline basis function N_0,...,N_n
   * where the i-th NURBS basis function is defined as
   *
   *         w_i P_i
   * R_i = -----------
   *            Q
   *
   * with P_i a basis function of a BSplineSpace
   * and Q an IgGridFunction built over a scalar BSpline space
   *
   */
  Assert(!phi.empty(), ExcEmptyObject());

  const auto n_pts = P.get_num_points();

  const auto &bsp_space = *bspline_elem.get_bspline_space();

  SafeSTLVector<Real> invQ(n_pts);

  int bsp_fn_id = 0;
  int offset = 0;
  for (int comp = 0 ; comp < n_components ; ++comp)
  {
//    Assert(bsp_space->get_num_basis(comp) == w_coefs.size(),
//           ExcDimensionMismatch(bsp_space->get_num_basis(comp),w_coefs.size()));

    for (int pt = 0 ; pt < n_pts ; ++pt)
      invQ[pt] = 1.0 / Q[pt](0);

    const int n_funcs_comp = bspline_elem.get_num_basis_comp(comp);
    for (int w_fn_id = 0 ; w_fn_id < n_funcs_comp ; ++w_fn_id, ++bsp_fn_id)
    {
      const auto &P_fn = P.get_function_view(bsp_fn_id);

      auto R_fn = phi.get_function_view(bsp_fn_id);

      const Real w = w_coefs[bsp_local_to_patch[bsp_fn_id]-offset];

      for (int pt = 0 ; pt < n_pts ; ++pt)
        R_fn[pt](comp) = P_fn[pt](comp) * invQ[pt] * w ;
    } // end loop w_fn_id
    offset += bsp_space.get_num_basis(comp);
  } // end loop comp

  phi.set_status_filled(true);
}


template<int dim_, int range_ , int rank_>
void
NURBSElementHandler<dim_, range_, rank_>::
FillCacheDispatcher::
evaluate_nurbs_gradients_from_bspline(
  const BSplineElem &bspline_elem,
  const SafeSTLVector<Index> &bsp_local_to_patch,
  const ValueTable<typename BSplineElem::Value> &P,
  const ValueTable<typename BSplineElem::template Derivative<1>> &dP,
  const IgCoefficients &w_coefs,
  const ValueVector<typename WeightFuncElem::Value> &Q,
  const ValueVector<typename WeightFuncElem::template Derivative<1>> &dQ,
  DataWithFlagStatus<ValueTable<Derivative<1>>> &D1_phi) const
{
  /*
   * This function evaluates the gradients of the n+1 NURBS basis function R_0,...,R_n
   * from the set of BSpline basis function P_0,...,P_n
   * where the i-th NURBS basis function is defined as
   *
   *         w_i P_i
   * R_i = ----------
   *            Q
   *
   * with P_i a basis function of a BSplineSpace
   * and Q an IgFunction built over a scalar space.
   *
   * Then the gradient dR_i is:
   *               _                     _
   *              |  dP_i       P_i * dQ  |
   * dR_i = w_i * | ------- -  ---------- |
   *              |    Q           Q^2    |
   *              |_                     _|
   */

  Assert(!D1_phi.empty(), ExcEmptyObject());

  const auto n_pts = P.get_num_points();

  const auto &bsp_space = *bspline_elem.get_bspline_space();

  SafeSTLVector<Real> invQ(n_pts);
  SafeSTLVector<SafeSTLArray<Real,dim>> dQ_invQ2(n_pts);

  int bsp_fn_id = 0;
  int offset = 0;
  for (int comp = 0 ; comp < n_components ; ++comp)
  {
//    Assert(bsp_space->get_num_basis(comp) == w_coefs.size(),
//           ExcDimensionMismatch(bsp_space->get_num_basis(comp),w_coefs.size()));

    for (int pt = 0 ; pt < n_pts ; ++pt)
    {
      invQ[pt] = 1.0 / Q[pt](0);

      const auto &dQ_pt = dQ[pt];
      auto &dQ_invQ2_pt = dQ_invQ2[pt];

      for (int i = 0 ; i < dim ; ++i)
        dQ_invQ2_pt[i] = invQ[pt] * invQ[pt] * dQ_pt(i)(0);
    }

    const int n_funcs_comp = bspline_elem.get_num_basis_comp(comp);
    for (int w_fn_id = 0 ; w_fn_id < n_funcs_comp ; ++w_fn_id, ++bsp_fn_id)
    {
      const auto &P_fn  =  P.get_function_view(bsp_fn_id);
      const auto &dP_fn = dP.get_function_view(bsp_fn_id);

      auto dR_fn = D1_phi.get_function_view(bsp_fn_id);

      const Real w = w_coefs[bsp_local_to_patch[bsp_fn_id]-offset];

      for (int pt = 0 ; pt < n_pts ; ++pt)
      {
        auto &dR_fn_pt = dR_fn[pt];

        const auto  &P_fn_pt =  P_fn[pt];
        const auto &dP_fn_pt = dP_fn[pt];

        const Real invQ_pt = invQ[pt];
        const auto &dQ_invQ2_pt = dQ_invQ2[pt];

        const auto &P_fn_pt_comp = P_fn_pt(comp);
        for (int i = 0 ; i < dim ; ++i)
          dR_fn_pt(i)(comp) = (dP_fn_pt(i)(comp)*invQ_pt - P_fn_pt_comp*dQ_invQ2_pt[i]) * w;
      } // end loop pt
    } // end loop w_fn_id
    offset += bsp_space.get_num_basis(comp);
  } // end loop comp

  D1_phi.set_status_filled(true);
}

template<int dim_, int range_ , int rank_>
void
NURBSElementHandler<dim_, range_, rank_>::
FillCacheDispatcher::
evaluate_nurbs_hessians_from_bspline(
  const BSplineElem &bspline_elem,
  const SafeSTLVector<Index> &bsp_local_to_patch,
  const ValueTable<typename BSplineElem::Value> &P,
  const ValueTable<typename BSplineElem::template Derivative<1>> &dP,
  const ValueTable<typename BSplineElem::template Derivative<2>> &d2P,
  const IgCoefficients &w_coefs,
  const ValueVector<typename WeightFuncElem::Value> &Q,
  const ValueVector<typename WeightFuncElem::template Derivative<1>> &dQ,
  const ValueVector<typename WeightFuncElem::template Derivative<2>> &d2Q,
  DataWithFlagStatus<ValueTable<Derivative<2>>> &D2_phi) const
{
  /*
   * This function evaluates the gradients of the n+1 NURBS basis function R_0,...,R_n
   * from the set of BSpline basis function N_0,...,N_n
   * where the k-th NURBS basis function is defined as
   *
   *        Pk
   * Rk = -------
   *         Q
   *
   * with Pk a basis function of a BSplineSpace
   * and Q an IgFunction built over a scalar space.
   *
   * Then the gradient dRk is defined by the partial derivatives:
   *
   *          dPk_i     Pk * dQ_i
   * dRk_i = ------- - -----------
   *            Q          Q^2
   *
   * And the hessian d2Rk is:
   *                                _                                         _
   *            d2Pk_ij      1     |                                           |    2 * Pk
   * d2Rk_ij = --------- - ----- * | dPk_i * dQ_j + dPk_j * dQ_i + Pk * d2Q_ij | + -------- * dQ_i * dQ_j
   *               Q        Q^2    |_                                         _|      Q^3
   *
   */
  Assert(!D2_phi.empty(), ExcEmptyObject());
  const auto n_pts = P.get_num_points();

  const auto &bsp_space = *bspline_elem.get_bspline_space();

  SafeSTLVector<Real> invQ(n_pts);
  SafeSTLVector<Real> invQ2(n_pts);
  SafeSTLVector<SafeSTLArray<Real,dim> > dQ_invQ2(n_pts);
  SafeSTLVector<SafeSTLArray<SafeSTLArray<Real,dim>,dim> > Q_terms_2nd_order(n_pts);

  int bsp_fn_id = 0;
  int offset = 0;
  for (int comp = 0 ; comp < n_components ; ++comp)
  {
//    Assert(bsp_space->get_num_basis(comp) == bsp_space.size(),
//           ExcDimensionMismatch(nrb_space->get_num_basis(comp),w_coefs.size()));

    for (int pt = 0 ; pt < n_pts ; ++pt)
    {
      auto &invQ_pt  = invQ[pt];
      auto &invQ2_pt = invQ2[pt];

      const auto &dQ_pt = dQ[pt];
      const auto &d2Q_pt = d2Q[pt];

      auto &dQ_invQ2_pt = dQ_invQ2[pt];
      auto &Q_terms_2nd_order_pt = Q_terms_2nd_order[pt];

      invQ_pt = 1.0 / Q[pt](0);
      invQ2_pt = invQ_pt * invQ_pt;

      const Real two_invQ_pt = 2.0 * invQ_pt;

      int hessian_entry_fid = 0;
      for (int i = 0 ; i < dim ; ++i)
      {
        dQ_invQ2_pt[i] = invQ2_pt * dQ_pt(i)(0);

        for (int j = 0 ; j < dim ; ++j, ++hessian_entry_fid)
          Q_terms_2nd_order_pt[i][j] = dQ_invQ2_pt[i] * dQ_pt(j)(0) * two_invQ_pt
                                       - d2Q_pt(hessian_entry_fid)(0) * invQ2_pt;
      } // end loop i

    } // end loop pt


    const int n_funcs_comp = bspline_elem.get_num_basis_comp(comp);
    for (int w_fn_id = 0 ; w_fn_id < n_funcs_comp ; ++w_fn_id, ++bsp_fn_id)
    {
      const auto &P_fn   =   P.get_function_view(bsp_fn_id);
      const auto &dP_fn  =  dP.get_function_view(bsp_fn_id);
      const auto &d2P_fn = d2P.get_function_view(bsp_fn_id);

      auto d2R_fn = D2_phi.get_function_view(bsp_fn_id);

      const Real w = w_coefs[bsp_local_to_patch[bsp_fn_id]-offset];

      for (int pt = 0 ; pt < n_pts ; ++pt)
      {
        auto &d2R_fn_pt = d2R_fn[pt];

        const auto   &P_fn_pt_comp = P_fn[pt](comp);
        const auto  &dP_fn_pt =  dP_fn[pt];
        const auto &d2P_fn_pt = d2P_fn[pt];

        const Real invQ_pt = invQ[pt];
        const auto &dQ_invQ2_pt = dQ_invQ2[pt];
        const auto &Q_terms_2nd_order_pt = Q_terms_2nd_order[pt];

        int hessian_entry_fid = 0;
        for (int i = 0 ; i < dim ; ++i)
        {
          const auto &dP_fn_pt_i_comp = dP_fn_pt(i)(comp);
          const auto &dQ_invQ2_pt_i = dQ_invQ2_pt[i];
          const auto &Q_terms_2nd_order_pt_i = Q_terms_2nd_order_pt[i];

          for (int j = 0 ; j < dim ; ++j, ++hessian_entry_fid)
          {
            const auto &dP_fn_pt_j_comp = dP_fn_pt(j)(comp);
            const auto &dQ_invQ2_pt_j = dQ_invQ2_pt[j];
            const auto &Q_terms_2nd_order_pt_i_j = Q_terms_2nd_order_pt_i[j];

            auto &d2R_fn_pt_i_j_comp = d2R_fn_pt(hessian_entry_fid)(comp);
            const auto &d2P_fn_pt_i_j_comp = d2P_fn_pt(hessian_entry_fid)(comp);

            d2R_fn_pt_i_j_comp = (d2P_fn_pt_i_j_comp * invQ_pt
                                  - dP_fn_pt_i_comp  * dQ_invQ2_pt_j
                                  - dP_fn_pt_j_comp * dQ_invQ2_pt_i
                                  +  P_fn_pt_comp * Q_terms_2nd_order_pt_i_j) * w;
          } // end loop j
        } // end loop i
      } // end loop pt
    } // end loop w_fn_id
    offset += bsp_space.get_num_basis(comp);
  } // end loop comp
  D2_phi.set_status_filled(true);
}


IGA_NAMESPACE_CLOSE

#include <igatools/basis_functions/nurbs_element_handler.inst>

#endif // #ifdef NURBS
