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

#include <igatools/basis_functions/bspline_element.h>
#include <igatools/basis_functions/bspline_element_handler.h>
//#include <igatools/basis_functions/bernstein_basis.h>


#include <algorithm>
using std::shared_ptr;

using std::set;

IGA_NAMESPACE_OPEN

namespace
{
//TODO(pauletti, Sep 9, 2014): this class seems to have a more general use put in own
// file under group tensor product values utilities
/**
 * @class DerivativeEvaluationSymmetryManager
 *
 * This class is used to manage the symmetries in the evaluation of a k-th order
 *  derivative of a scalar function \f$ f \mathbb{R}^d -to \mathbb{R} \f$
 *  according with the Schwarz theorem (equality of mixed partials).
 *
 * The total number of derivatives is given by \f$ d^k \f$ where the number of
 * different values
 * is \f$ \binorm{d-1+k}{d-1} = \binorm{d-1+k}{k} \f$.
 *
 */
template<int size>
SafeSTLVector<TensorIndex<size> >
partition(const int n)
{
  SafeSTLVector<TensorIndex<size>> v;
  TensorIndex<size> arr(0);

  arr[0] = n;
  v.push_back(arr);

  for (int j=1; j<n+1; ++j)
  {
    auto w = partition<size-1>(j);
    for (auto a : w)
    {
      arr[0] = n-j;
      std::copy(a.begin(), a.end(), arr.begin()+1);
      v.push_back(arr);
    }
  }
  return v;
}

template<>
SafeSTLVector<TensorIndex<1> >
partition<1>(const int n)
{
  TensorIndex<1> arr(n);
  return SafeSTLVector<TensorIndex<1>>(1,arr);
}



template<>
SafeSTLVector<TensorIndex<0> >
partition<0>(const int n)
{
  return SafeSTLVector<TensorIndex<0>>();
}



template<int dim, int order>
class TensorFunctionDerivativesSymmetry
{
public:
//    static const int num_entries_total = pow(dim,order);
  static const int num_entries_eval = constexpr_binomial_coefficient(dim-1+order,order);

  TensorFunctionDerivativesSymmetry()
  {
    auto uni_indices = partition<dim>(order);
    std::copy(uni_indices.begin(), uni_indices.end(), univariate_order.begin());



    for (int j=0; j<num_entries_eval; ++j)
    {
      auto &der_ind = eval_indices[j];
      int s=0;
      for (int dir=0; dir<dim; ++dir)
      {
        for (int l=0; l<uni_indices[j][dir]; ++l)
          der_ind[s+l] = dir;
        s += uni_indices[j][dir];
      }

      auto ind = sequence<order>();
      SafeSTLVector<TensorIndex<order>> v;
      do
      {
        TensorIndex<order> ti;
        for (int i=0; i<order; ++i)
          ti[i] = eval_indices[j][ind[i]];
        v.push_back(ti);
      }
      while (std::next_permutation(ind.begin(),ind.end()));

      auto it = std::unique(v.begin(), v.end());
      v.resize(std::distance(v.begin(),it));

      copy_indices[j] = v;
    }
  }

  void print_info(LogStream &out) const
  {
    out.begin_item("univariate derivative orders:");
    univariate_order.print_info(out);
    out.end_item();

    out.begin_item("Assigment indices:");
    eval_indices.print_info(out);
    out.end_item();

    out.begin_item("all equal indices indices:");
    copy_indices.print_info(out);
    out.end_item();
  }
  SafeSTLArray<TensorIndex<dim>, num_entries_eval> univariate_order;

  SafeSTLArray<TensorIndex<order>, num_entries_eval> eval_indices;

  SafeSTLArray<SafeSTLVector<TensorIndex<order>>, num_entries_eval> copy_indices;

};


}; // of the namespace




template<int dim_, int range_ , int rank_>
BSplineElementHandler<dim_, range_, rank_>::
BSplineElementHandler(shared_ptr<const Basis> space)
  :
  base_t(space)
{}


#if 0
template<int dim_, int range_ , int rank_>
auto
BSplineElementHandler<dim_, range_, rank_>::
create(std::shared_ptr<const Basis> space) -> std::unique_ptr<self_t>
{
  auto handler = std::unique_ptr<self_t>(new self_t(space));
  Assert(handler != nullptr, ExcNullPtr());
  return handler;
}
#endif


template<int dim_, int range_ , int rank_>
template<int sdim>
void
BSplineElementHandler<dim_, range_, rank_>::
SetFlagDispatcher::
operator()(const Topology<sdim> &topology)
{
  using GridFlags = grid_element::Flags;
  //TODO (martinelli, Aug 27, 2015): select the proper grid flags depending on the BSpline element flags
  const auto grid_flags = GridFlags::point |
                          GridFlags::weight;
  grid_handler_.template set_flags<sdim>(grid_flags);

  flags_[sdim] = flag_in_;
}


template<int dim_, int range_ , int rank_>
void
BSplineElementHandler<dim_, range_, rank_>::
set_flags_impl(const topology_variant &topology,
               const typename space_element::Flags &flag)
{
  auto elem_flags = flag;
  if (contains(flag,space_element::Flags::divergence))
    elem_flags |= space_element::Flags::gradient;


  auto set_flag_dispatcher = SetFlagDispatcher(elem_flags,this->grid_handler_,this->flags_);
  boost::apply_visitor(set_flag_dispatcher,topology);
}


template<int dim_, int range_ , int rank_>
template<int sdim>
void
BSplineElementHandler<dim_, range_, rank_>::
InitCacheDispatcher::
init_cache_1D()
{
  const auto &quad = *bsp_elem_.get_grid_element().template get_quad<sdim>();
  const auto &bsp_basis = dynamic_cast<const Basis &>(*bsp_elem_.get_space_basis());

  const auto &spline_space = *bsp_basis.spline_space_;

  const auto &degree = spline_space.get_degree_table();
  const auto &active_components_id = spline_space.get_active_components_id();

  const auto n_pts = quad.get_num_coords_direction();

  auto &splines_1D_table = bsp_elem_.all_splines_1D_table_[sdim];

  const int n_sub_elems = UnitElement<dim>::template num_elem<sdim>();
  splines_1D_table.resize(n_sub_elems);

  for (auto s_id = 0 ; s_id < n_sub_elems ; ++s_id)
  {
    auto &splines_1D_table_sub_elem = splines_1D_table[s_id];
    splines_1D_table_sub_elem = typename BSplineElem::Splines1DTable(spline_space.get_components_map());

    const auto &sub_elem = UnitElement<dim>::template get_elem<sdim>(s_id);
    TensorSize<dim> n_coords(1);
    for (int dir = 0 ; dir < sdim ; ++dir)
      n_coords[sub_elem.active_directions[dir]] = n_pts[dir];


    for (auto comp : active_components_id)
    {
      auto &splines_1D_comp = splines_1D_table_sub_elem[comp];

      const auto &deg_comp = degree[comp];

      for (const int dir : UnitElement<dim>::active_directions)
        splines_1D_comp[dir].resize(deg_comp[dir]+1,n_coords[dir]);
    } // end loop comp
  } // end loop sub_elem
}


template<int dim_, int range_ , int rank_>
template<int sdim>
void
BSplineElementHandler<dim_, range_, rank_>::
InitCacheDispatcher::
init_cache_multiD()
{
  auto &cache = bsp_elem_.all_sub_elems_cache_;

  const auto n_basis = bsp_elem_.get_basis_offset()[BaseSpace::n_components];

  const auto n_pts = bsp_elem_.get_grid_element().template get_quad<sdim>()->get_num_points();

  const auto flag = flags_[sdim];

  for (auto &s_id: UnitElement<dim_>::template elems_ids<sdim>())
  {
    auto &s_cache = cache.template get_sub_elem_cache<sdim>(s_id);
    s_cache.resize(flag, n_pts, n_basis);
  }
}

template<int dim_, int range_ , int rank_>
template<int sdim>
void
BSplineElementHandler<dim_, range_, rank_>::
InitCacheDispatcher::
operator()(const std::shared_ptr<const Quadrature<sdim>> &quad)
{
  grid_handler_.template init_cache<sdim>(bsp_elem_.get_grid_element(),quad);

  init_cache_1D<sdim>();

  init_cache_multiD<sdim>();
}


template<int dim_, int range_ , int rank_>
void
BSplineElementHandler<dim_, range_, rank_>::
init_cache_impl(BaseElem &elem,
                const eval_pts_variant &quad) const
{
  auto init_cache_dispatcher =
    InitCacheDispatcher(
      this->grid_handler_,
      this->flags_,
      dynamic_cast<BSplineElem &>(elem));
  boost::apply_visitor(init_cache_dispatcher,quad);
}



template <int dim, int range, int rank>
void
BSplineElementHandler<dim, range, rank>::
FillCacheDispatcherNoGlobalCache::
copy_to_inactive_components_values(const SafeSTLVector<Index> &inactive_comp,
                                   const SafeSTLArray<Index, n_components> &active_map,
                                   ValueTable<Value> &D_phi) const
{
  const auto &comp_offset = bsp_elem_.get_basis_offset();

  Assert(D_phi.get_num_functions() == comp_offset[BaseSpace::n_components],
         ExcDimensionMismatch(D_phi.get_num_functions(),
                              comp_offset[BaseSpace::n_components]));


  const Size n_points = D_phi.get_num_points();
  for (int comp : inactive_comp)
  {
    const auto act_comp = active_map[comp];
    const auto n_basis_comp = bsp_elem_.get_num_basis_comp(comp);
    const Size act_offset = comp_offset[act_comp];
    const Size offset     = comp_offset[comp];
    for (Size basis_i = 0; basis_i < n_basis_comp;  ++basis_i)
    {
      const auto act_D_phi = D_phi.get_function_view(act_offset+basis_i);

      auto inact_D_phi = D_phi.get_function_view(offset+basis_i);

      for (int qp = 0; qp < n_points; ++qp)
        inact_D_phi[qp](comp) = act_D_phi[qp](act_comp);

    } // end loop basis_i
  } // end loop comp
}



template <int dim, int range, int rank>
template <int order>
void
BSplineElementHandler<dim, range, rank>::
FillCacheDispatcherNoGlobalCache::
copy_to_inactive_components(const SafeSTLVector<Index> &inactive_comp,
                            const SafeSTLArray<Index, n_components> &active_map,
                            ValueTable<Derivative<order>> &D_phi) const
{
  const auto &comp_offset = bsp_elem_.get_basis_offset();

  Assert(D_phi.get_num_functions() == comp_offset[BaseSpace::n_components],
         ExcDimensionMismatch(D_phi.get_num_functions(),
                              comp_offset[BaseSpace::n_components]));

  const Size n_points = D_phi.get_num_points();
  const Size n_ders = Derivative<order>::size;
  for (int comp : inactive_comp)
  {
    const auto act_comp = active_map[comp];
    const auto n_basis_comp = bsp_elem_.get_num_basis_comp(comp);
    const Size act_offset = comp_offset[act_comp];
    const Size offset     = comp_offset[comp];
    for (Size basis_i = 0; basis_i < n_basis_comp;  ++basis_i)
    {
      const auto act_D_phi = D_phi.get_function_view(act_offset+basis_i);
      auto     inact_D_phi = D_phi.get_function_view(offset+basis_i);

      for (int pt = 0; pt < n_points; ++pt)
      {
        const auto &act_D_phi_pt = act_D_phi[pt];
        auto &inact_D_phi_pt = inact_D_phi[pt];

        for (int der = 0; der < n_ders; ++der)
          inact_D_phi_pt(der)(comp) = act_D_phi_pt(der)(act_comp);
      } // end loop pt
    } // end loop basis_i
  } // end loop comp
}







template <int dim, int range, int rank>
void
BSplineElementHandler<dim, range, rank>::
FillCacheDispatcherNoGlobalCache::
evaluate_bspline_values(
  const ComponentContainer<std::unique_ptr<const TensorProductFunctionEvaluator<dim>>> &elem_values,
  ValueTable<Value> &phi) const
{
  const auto &comp_offset = bsp_elem_.get_basis_offset();

  Assert(phi.get_num_functions() == comp_offset[BaseSpace::n_components],
         ExcDimensionMismatch(phi.get_num_functions(),
                              comp_offset[BaseSpace::n_components]));

  const Size n_points = phi.get_num_points();
  const TensorIndex<dim> der_tensor_id; // [0,0,..,0] tensor index
  for (int comp : elem_values.get_active_components_id())
  {
    const auto &values = *elem_values[comp];
    const int n_basis_comp = bsp_elem_.get_num_basis_comp(comp);
    const Size offset = comp_offset[comp];

    for (int func_id = 0; func_id < n_basis_comp; ++func_id)
    {
      auto phi_i = phi.get_function_view(offset + func_id);
      auto const &func_tensor_id = values.func_flat_to_tensor(func_id);

      for (int pt = 0; pt < n_points; ++pt)
        phi_i[pt](comp) = values.evaluate(der_tensor_id, func_tensor_id, pt);
    } // end func_id loop
  } // end comp loop


  copy_to_inactive_components_values(elem_values.get_inactive_components_id(),
                                     elem_values.get_comp_map(),
                                     phi);
}




template <int dim, int range, int rank>
template <int order>
void
BSplineElementHandler<dim, range, rank>::
FillCacheDispatcherNoGlobalCache::
evaluate_bspline_derivatives(
  const ComponentContainer<std::unique_ptr<const TensorProductFunctionEvaluator<dim>>> &elem_values,
  ValueTable<Derivative<order>> &D_phi) const
{
  /*
   * This code computes any order of derivatives for a multivariate
   * B-spline on the current element
   * We use the formula
   * \partial_(\alpha_1,...,\alpha_n) B(qp) = \Pi d^{\alpha_i} B_i(qp_i)
   */

  Assert(D_phi.size() > 0, ExcEmptyObject());
  const int n_points = D_phi.get_num_points();


  const auto &comp_offset = bsp_elem_.get_basis_offset();

  Assert(D_phi.get_num_functions() == comp_offset[BaseSpace::n_components],
         ExcDimensionMismatch(D_phi.get_num_functions(),
                              comp_offset[BaseSpace::n_components]));

  TensorFunctionDerivativesSymmetry<dim,order> sym;
  const auto n_der = TensorFunctionDerivativesSymmetry<dim,order>::num_entries_eval;

  const auto &univariate_order = sym.univariate_order ;
  const auto &copy_indices = sym.copy_indices;

  for (int comp : elem_values.get_active_components_id())
  {
    const auto &values = *elem_values[comp];
    const int n_basis_comp = bsp_elem_.get_num_basis_comp(comp);
    const int offset = comp_offset[comp];

    for (int func_id = 0 ; func_id < n_basis_comp; ++func_id)
    {
      auto D_phi_i = D_phi.get_function_view(offset + func_id);

      auto const &func_tensor_id = values.func_flat_to_tensor(func_id);

      for (int der_id = 0 ; der_id < n_der ; ++der_id)
      {
        auto const &der_tensor_id = univariate_order[der_id];

        const auto &copy_indices_der = copy_indices[der_id];
//                const auto copy_indices_der_size = copy_indices_der.size();

        const auto &copy_indices_der_0 = copy_indices_der[0];

        const auto copy_indices_der_end = copy_indices_der.end();

        for (int pt = 0; pt < n_points; ++pt)
        {
          auto &der = D_phi_i[pt];
          Real &der_copy_indices_der_0_comp = der(copy_indices_der_0)(comp);

          der_copy_indices_der_0_comp = values.evaluate(der_tensor_id, func_tensor_id, pt);

//                    for (int k = 1 ; k < copy_indices_der_size ; ++k)
//                        der(copy_indices_der[k])(comp) = der_copy_indices_der_0_comp;

          auto copy_indices_der_it = copy_indices_der.begin()+1;
          for (; copy_indices_der_it != copy_indices_der_end ; ++copy_indices_der_it)
            der(*copy_indices_der_it)(comp) = der_copy_indices_der_0_comp;

        } // end loop pt

      } // end loop der_id

    } // end loop func_id

  } // end comp loop

  copy_to_inactive_components<order>(elem_values.get_inactive_components_id(),
                                     elem_values.get_comp_map(), D_phi);
}


template<int dim_, int range_ , int rank_>
template<int sdim>
void
BSplineElementHandler<dim_, range_, rank_>::
FillCacheDispatcherNoGlobalCache::
fill_cache_1D(const Quadrature<dim> &extended_sub_elem_quad)
{
  auto &grid_elem = bsp_elem_.get_grid_element();



  //--------------------------------------------------------------------------------------
  // filling the 1D cache --- begin

  const auto &grid = *grid_elem.get_grid();
  const auto n_inter = grid.get_num_intervals();

  const auto elem_size = grid_elem.template get_side_lengths<dim>(0);
  const auto &elem_tensor_id = grid_elem.get_index().get_tensor_index();

  const auto &bsp_basis = dynamic_cast<const Basis &>(*bsp_elem_.get_space_basis());

  const auto &spline_space = *bsp_basis.spline_space_;

  const auto &degree = spline_space.get_degree_table();

  const auto &active_components_id = spline_space.get_active_components_id();

  const auto &n_coords = extended_sub_elem_quad.get_num_coords_direction();

  const auto &bezier_op   = bsp_basis.operators_;
  const auto &end_interval = bsp_basis.end_interval_;

  auto &splines_1D_table_subelems = bsp_elem_.all_splines_1D_table_[sdim];
  auto &splines_1D_table = splines_1D_table_subelems[s_id_];

  for (const int dir : UnitElement<dim>::active_directions)
  {
    const auto &pt_coords_internal = extended_sub_elem_quad.get_coords_direction(dir);

    const auto len = elem_size[dir];

    const auto interval_id = elem_tensor_id[dir];

    Real alpha;

    const int n_pts_1D = n_coords[dir];

    SafeSTLVector<Real> pt_coords_boundary(n_pts_1D);

    const SafeSTLVector<Real> *pt_coords_ptr = nullptr;

    for (auto comp : active_components_id)
    {
      if (interval_id == 0) // processing the leftmost interval
      {
        // first interval (i.e. left-most interval)

        alpha = end_interval[comp][dir].first;
        const Real one_minus_alpha = 1. - alpha;

        for (int ipt = 0 ; ipt < n_pts_1D ; ++ipt)
          pt_coords_boundary[ipt] = one_minus_alpha +
                                    pt_coords_internal[ipt] * alpha;

        pt_coords_ptr = &pt_coords_boundary;
      } // end process_interval_left
      else if (interval_id == n_inter[dir]-1) // processing the rightmost interval
      {
        // last interval (i.e. right-most interval)

        alpha = end_interval[comp][dir].second;

        for (int ipt = 0 ; ipt < n_pts_1D ; ++ipt)
          pt_coords_boundary[ipt] = pt_coords_internal[ipt] *
                                    alpha;

        pt_coords_ptr = &pt_coords_boundary;
      } // end process_interval_right
      else
      {
        // internal interval

        alpha = 1.0;

        pt_coords_ptr = &pt_coords_internal;
      } // end process_interval_internal


      const Real alpha_div_interval_length = alpha / len;

      const auto &oper = bezier_op.get_operator(dir,interval_id,comp);

      //------------------------------------------------------------
      //resize_and_fill_bernstein_values
      const int deg = degree[comp][dir];

      auto &splines_1D_comp = splines_1D_table[comp];
      auto &splines_1D = splines_1D_comp[dir];
      for (int order = 0; order < MAX_NUM_DERIVATIVES; ++order)
      {
        const Real scaling_factor = std::pow(alpha_div_interval_length, order);
        const auto &bernstein_values = BernsteinBasis::derivative(order, deg,*pt_coords_ptr);

        auto &splines = splines_1D.get_derivative(order);
        splines = oper.scale_action(scaling_factor,bernstein_values);
      } // end loop order
      //------------------------------------------------------------

    } // end loop comp

  } // end loop dir
  //
  // filling the 1D cache --- end
  //-------------------------------------------------------------------------------
}

template<int dim_, int range_ , int rank_>
template<int sdim>
void
BSplineElementHandler<dim_, range_, rank_>::
FillCacheDispatcherNoGlobalCache::
fill_cache_multiD(const Quadrature<dim> &extended_sub_elem_quad)
{
  //-------------------------------------------------------------------------------
  // Multi-variate spline evaluation from 1D values --- begin
  const auto &splines_1D_table_subelems = bsp_elem_.all_splines_1D_table_[sdim];
  const auto &splines_1D_table = splines_1D_table_subelems[s_id_];

  using TPFE = const TensorProductFunctionEvaluator<dim>;
  ComponentContainer<std::unique_ptr<TPFE>> val_1d(splines_1D_table.get_comp_map());

  for (auto c : val_1d.get_active_components_id())
    val_1d[c] = std::make_unique<TPFE>(extended_sub_elem_quad,splines_1D_table[c]);

  // Multi-variate spline evaluation from 1D values --- end
  //-------------------------------------------------------------------------------




  //-------------------------------------------------------------------------------
  auto &sub_elem_cache =
    bsp_elem_.all_sub_elems_cache_.template get_sub_elem_cache<sdim>(s_id_);


  using Elem = SpaceElement<dim_,0,range_,rank_>;
  using _Value      = typename Elem::_Value;
  using _Gradient   = typename Elem::_Gradient;
  using _Hessian    = typename Elem::_Hessian;
  using _Divergence = typename Elem::_Divergence;

  if (sub_elem_cache.template status_fill<_Value>())
  {
    auto &values = sub_elem_cache.template get_data<_Value>();
    evaluate_bspline_values(val_1d, values);
    values.set_status_filled(true);
  }
  if (sub_elem_cache.template status_fill<_Gradient>())
  {
    auto &gradients = sub_elem_cache.template get_data<_Gradient>();
    evaluate_bspline_derivatives<1>(val_1d, gradients);
    gradients.set_status_filled(true);
  }
  if (sub_elem_cache.template status_fill<_Hessian>())
  {
    auto &hessians = sub_elem_cache.template get_data<_Hessian>();
    evaluate_bspline_derivatives<2>(val_1d, hessians);
    hessians.set_status_filled(true);
  }
  if (sub_elem_cache.template status_fill<_Divergence>())
  {
    auto &divergences = sub_elem_cache.template get_data<_Divergence>();
    eval_divergences_from_gradients(
      sub_elem_cache.template get_data<_Gradient>(),divergences);
    divergences.set_status_filled(true);
  }

  sub_elem_cache.set_filled(true);
  //-------------------------------------------------------------------------------
}

template<int dim_, int range_ , int rank_>
template<int sdim>
void
BSplineElementHandler<dim_, range_, rank_>::
FillCacheDispatcherNoGlobalCache::
operator()(const Topology<sdim> &topology)
{
  auto &grid_elem = bsp_elem_.get_grid_element();
  grid_handler_.template fill_cache<sdim>(grid_elem,s_id_);

  const auto extended_sub_elem_quad =
    extend_sub_elem_quad<sdim,dim>(
      *grid_elem.template get_quad<sdim>(),
      s_id_);

  fill_cache_1D<sdim>(extended_sub_elem_quad);

  fill_cache_multiD<sdim>(extended_sub_elem_quad);
}



template<int dim_, int range_ , int rank_>
void
BSplineElementHandler<dim_, range_, rank_>::
fill_cache_impl(const topology_variant &topology,
                BaseElem &elem,
                const int s_id) const
{
  auto fill_cache_dispatcher =
    FillCacheDispatcherNoGlobalCache(
      s_id,
      this->grid_handler_,
      dynamic_cast<BSplineElem &>(elem));
  boost::apply_visitor(fill_cache_dispatcher,topology);
}

template<int dim_, int range_ , int rank_>
auto
BSplineElementHandler<dim_, range_, rank_>::
get_bspline_basis() const -> std::shared_ptr<const Basis>
{
  auto bsp_basis = std::dynamic_pointer_cast<const Basis>(this->get_space());
  Assert(bsp_basis != nullptr,ExcNullPtr());
  return bsp_basis;
}





template<int dim_, int range_ , int rank_>
void
BSplineElementHandler<dim_, range_, rank_>::
print_info(LogStream &out) const
{
  /*
  out.begin_item("Grid Cache:");
  this->grid_handler_.print_info(out);
  out.end_item();
  //*/

#if 0
  out.begin_item("Splines 1D Cache:");
  cacheutils::print_caches(splines1d_, out);
//    splines1d_.print_info(out);
  out.end_item();
#endif
}



#if 0
template<int dim_, int range_ , int rank_>
BSplineElementHandler<dim_, range_, rank_>::
GlobalCache::
GlobalCache(const std::shared_ptr<const Quadrature<dim>> &quad, const ComponentMap &component_map)
  :
  quad_(quad),
  basis_values_1d_table_(BasisValues1dTable(component_map))
{}

template<int dim_, int range_ , int rank_>
auto
BSplineElementHandler<dim_, range_, rank_>::
GlobalCache::
entry(const int comp, const int dir, const Index interval_id) -> BasisValues1d &
{
  return basis_values_1d_table_[comp][dir][interval_id];
}


template<int dim_, int range_ , int rank_>
void
BSplineElementHandler<dim_, range_, rank_>::
GlobalCache::
print_info(LogStream &out) const
{
  using std::to_string;
  for (const auto comp : basis_values_1d_table_.get_active_components_id())
  {
    out.begin_item("Active Component ID: " + to_string(comp));

    for (const int dir : UnitElement<dim_>::active_directions)
    {
      out.begin_item("Direction : " + to_string(dir));

      for (const auto &interv_id_and_basis : basis_values_1d_table_[comp][dir])
      {
        const auto interval_id = interv_id_and_basis.first;
        const auto &basis = interv_id_and_basis.second;

        out.begin_item("Interval ID: " + to_string(interval_id));
        basis.print_info(out);
        out.end_item();
      }
      out.end_item();
    } // end loop dir
    out.end_item();
  } // end loop comp
}


template<int dim_, int range_ , int rank_>
auto
BSplineElementHandler<dim_, range_, rank_>::
GlobalCache::
get_element_values(const TensorIndex<dim> &elem_tensor_id) const
-> ComponentContainer<std::unique_ptr<const TensorProductFunctionEvaluator<dim>>>
{
  ComponentContainer<std::unique_ptr<const TensorProductFunctionEvaluator<dim>> >
  result(basis_values_1d_table_.get_comp_map());

  SafeSTLArray<BasisValues1dConstView, dim> values_1D;
  for (auto c : result.get_active_components_id())
  {
    const auto &value = basis_values_1d_table_[c];
    for (int i = 0 ; i < dim_ ; ++i)
      values_1D[i] = BasisValues1dConstView(value[i].at(elem_tensor_id[i]));

    result[c] = std::make_unique<const TensorProductFunctionEvaluator<dim>>(*this->quad_,values_1D);
  }
  return result;
}
#endif


IGA_NAMESPACE_CLOSE

#include <igatools/basis_functions/bspline_element_handler.inst>
