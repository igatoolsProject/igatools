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


#ifndef NURBS_ELEMENT_HANDLER_H_
#define NURBS_ELEMENT_HANDLER_H_

#include <igatools/base/config.h>

#ifdef USE_NURBS

#include <igatools/base/cache_status.h>
#include <igatools/base/flags_handler.h>
#include <igatools/base/quadrature.h>


#include <igatools/utils/value_table.h>
#include <igatools/geometry/grid_handler.h>
#include <igatools/basis_functions/bspline_element_handler.h>
#include <igatools/functions/ig_grid_function_handler.h>




IGA_NAMESPACE_OPEN


template<int,int,int> class NURBS;
template<int,int,int> class NURBSElement;

/**
 * Global NURBS uniform quadrature
 * computational optimization cache.
 *
 * @ingroup handlers
 */
template<int dim_, int range_, int rank_>
class NURBSElementHandler
  : public ReferenceElementHandler<dim_,range_,rank_>
{
  using base_t = ReferenceElementHandler<dim_,range_,rank_>;
  using self_t = NURBSElementHandler<dim_,range_,rank_>;
  using Basis = NURBS<dim_,range_,rank_>;
  static const Size n_components =  Basis::n_components;

  using IndexType = typename Grid<dim_>::IndexType;

  template<class T>
  using ComponentContainer = typename Basis::template ComponentContainer<T>;

  template <int order>
  using Derivative = typename Basis::template Derivative<order>;

  using Value = typename Basis::Value;

protected:
  using ElementIterator = typename Basis::ElementIterator;
  using ElementAccessor = typename Basis::ElementAccessor;

  using BaseSpace = ReferenceSpaceBasis<dim_,range_,rank_>;
  using RefElementIterator = typename BaseSpace::ElementIterator;
  using RefElementAccessor = typename BaseSpace::ElementAccessor;


public:
  static const int dim = dim_;



public:
  /** @name Constructors.*/
  ///@{
  /**
   * Default constructor. Not allowed to be used.
   */
  NURBSElementHandler() = delete;

  NURBSElementHandler(std::shared_ptr<const Basis> space);

  /**
   * Copy constructor. Not allowed to be used.
   */
  NURBSElementHandler(const self_t &) = delete;

  /**
   * Move constructor. Not allowed to be used.
   */
  NURBSElementHandler(self_t &&) = delete;
  ///@}

  /**
   * Assignment operators.
   */
  ///@{
  /**
   * Copy assignment operator. Not allowed to be used.
   */
  self_t &operator=(const self_t &) = delete;

  /**
   * Move assignment operator. Not allowed to be used.
   */
  self_t &operator=(self_t &&) = delete;
  ///@}

  /**
   * Destructor.
   */
  virtual ~NURBSElementHandler() = default;


  using topology_variant = typename base_t::topology_variant;
  using eval_pts_variant = typename base_t::eval_pts_variant;



  virtual void print_info(LogStream &out) const override final;

private:


  std::unique_ptr<SpaceElementHandler<dim_,0,range_,rank_>> bsp_elem_handler_;

  std::unique_ptr<IgGridFunctionHandler<dim_,1>> w_func_elem_handler_;


  using WeightElem = typename Basis::WeightFunction::ElementAccessor;
//  using WeightElemTable = typename Basis::template ComponentContainer<std::shared_ptr<WeightElem>>;



  /**
   * Returns the NURBS used to define the NURBSElementHandler object.
   */
  std::shared_ptr<const Basis> get_nurbs_space() const;



private:

  using BaseElem = SpaceElement<dim_,0,range_,rank_>;

  virtual void set_flags_impl(const topology_variant &topology,
                              const typename space_element::Flags &flag) override final;

  virtual void init_cache_impl(BaseElem &elem,
                               const eval_pts_variant &quad) const override final;

  virtual void fill_cache_impl(const topology_variant &topology,
                               BaseElem &elem,
                               const int s_id) const override final;



  struct SetFlagDispatcher : boost::static_visitor<void>
  {
    SetFlagDispatcher(const typename space_element::Flags nrb_flag,
                      self_t &nrb_handler);

    template<int sdim>
    void operator()(const Topology<sdim> &topology)
    {
      nrb_handler_.flags_[sdim] = nrb_flag_;
    }


  private:
    const typename space_element::Flags nrb_flag_;
    self_t &nrb_handler_;
  };


  struct InitCacheDispatcher : boost::static_visitor<void>
  {
    InitCacheDispatcher(const self_t &nrb_handler,
                        SpaceElement<dim_,0,range_,rank_> &elem);

    template<int sdim>
    void operator()(const std::shared_ptr<const Quadrature<sdim>> &quad)
    {
//      using NURBSElem = NURBSElement<dim_,range_,rank_>;
//      auto &nrb_elem = dynamic_cast<NURBSElem &>(elem_);

      auto &bsp_elem = *nrb_elem_.bspline_elem_;
      nrb_handler_.bsp_elem_handler_->template init_cache<sdim>(bsp_elem,quad);

      auto &w_func_elem = *(nrb_elem_.weight_elem_);
      nrb_handler_.w_func_elem_handler_->init_cache(w_func_elem,quad);

      auto &cache = nrb_handler_.get_element_cache(nrb_elem_);

      const auto n_basis = nrb_elem_.get_num_basis(DofProperties::active);

      Assert(quad == nrb_elem_.get_grid_element().template get_quad<sdim>(),
             ExcMessage("Different quadratures."));
      const auto n_pts = quad->get_num_points();

      const auto flag = nrb_handler_.flags_[sdim];

      for (auto &s_id: UnitElement<dim_>::template elems_ids<sdim>())
      {
        auto &s_cache = cache.template get_sub_elem_cache<sdim>(s_id);
        s_cache.resize(flag, n_pts, n_basis);
      }
    }


  private:
    const self_t &nrb_handler_;
    NURBSElement<dim_,range_,rank_> &nrb_elem_;
  };


  struct FillCacheDispatcher : boost::static_visitor<void>
  {
    FillCacheDispatcher(const self_t &nrb_handler,
                        SpaceElement<dim_,0,range_,rank_> &elem,
                        const int s_id);

    template<int sdim>
    void operator()(const Topology<sdim> &topology)
    {
//      static_assert(sdim == dim,"The case with sdim != dim is not implemented!");
      Assert(sdim == dim,ExcNotImplemented());

      auto &bsp_elem = *nrb_elem_.bspline_elem_;
      nrb_handler_.bsp_elem_handler_->template fill_cache<sdim>(bsp_elem,s_id_);

      auto &w_func_elem = *(nrb_elem_.weight_elem_);
      nrb_handler_.w_func_elem_handler_->fill_cache(topology,w_func_elem,s_id_);

      auto &cache =
        nrb_handler_.get_element_cache(nrb_elem_).template get_sub_elem_cache<sdim>(s_id_);

      using space_element::_Value;
      using space_element::_Gradient;
      using space_element::_Hessian;
      using _D0 = grid_function_element::template _D<0>;
      using _D1 = grid_function_element::template _D<1>;
      using _D2 = grid_function_element::template _D<2>;

      const auto bsp_local_to_patch = bsp_elem.get_local_to_patch(DofProperties::active);

      const auto &w_coefs =
        nrb_handler_.w_func_elem_handler_->get_ig_grid_function()->get_coefficients();

      if (cache.template status_fill<_Value>())
      {
        const auto &P = bsp_elem.template get_basis_data<_Value,sdim>(s_id_,DofProperties::active);
        const auto &Q = w_func_elem.template get_values_from_cache<_D0,sdim>(s_id_);

        auto &values = cache.template get_data<_Value>();
        evaluate_nurbs_values_from_bspline(bsp_elem,bsp_local_to_patch,P,w_coefs,Q,values);
      }

      if (cache.template status_fill<_Gradient>())
      {
        const auto &P = bsp_elem.template get_basis_data<_Value,sdim>(s_id_,DofProperties::active);
        const auto &dP = bsp_elem.template get_basis_data<_Gradient,sdim>(s_id_,DofProperties::active);
        const auto &Q = w_func_elem.template get_values_from_cache<_D0,sdim>(s_id_);
        const auto &dQ = w_func_elem.template get_values_from_cache<_D1,sdim>(s_id_);

        auto &gradients = cache.template get_data<_Gradient>();
        evaluate_nurbs_gradients_from_bspline(bsp_elem,bsp_local_to_patch,P,dP,w_coefs,Q,dQ,gradients);
      }

      if (cache.template status_fill<_Hessian>())
      {
        const auto &P = bsp_elem.template get_basis_data<_Value,sdim>(s_id_,DofProperties::active);
        const auto &dP = bsp_elem.template get_basis_data<_Gradient,sdim>(s_id_,DofProperties::active);
        const auto &d2P = bsp_elem.template get_basis_data<_Hessian,sdim>(s_id_,DofProperties::active);
        const auto &Q = w_func_elem.template get_values_from_cache<_D0,sdim>(s_id_);
        const auto &dQ = w_func_elem.template get_values_from_cache<_D1,sdim>(s_id_);
        const auto &d2Q = w_func_elem.template get_values_from_cache<_D2,sdim>(s_id_);
        auto &hessians = cache.template get_data<_Hessian>();
        evaluate_nurbs_hessians_from_bspline(bsp_elem,bsp_local_to_patch,P,dP,d2P,w_coefs,Q,dQ,d2Q, hessians);
      }

      using space_element::_Divergence;
      if (cache.template status_fill<_Divergence>())
      {
        const auto &gradient = cache.template get_data<_Gradient>();
        auto &divergence = cache.template get_data<_Divergence>();
        eval_divergences_from_gradients(gradient,divergence);
        divergence.set_status_filled(true);
      }
      cache.set_filled(true);
    }

  private:
    const self_t &nrb_handler_;
    NURBSElement<dim_,range_,rank_> &nrb_elem_;
    const int s_id_;


    using BSplineElem = BSplineElement<dim_,range_,rank_>;
    using WeightFuncElem = GridFunctionElement<dim_,1>;

    /**
     * Computes the value of the non-zero NURBS basis
     * functions over the current element,
     *   at the evaluation points pre-allocated in the cache.
     *
     * \warning If the output result @p derivatives_phi_hat is not correctly pre-allocated,
     * an exception will be raised.
     */
    void evaluate_nurbs_values_from_bspline(
      const BSplineElem &bspline_elem,
      const SafeSTLVector<Index> &bsp_local_to_patch,
      const ValueTable<typename BSplineElem::Value> &N,
      const IgCoefficients &w_coeffs,
      const ValueVector<typename WeightFuncElem::Value> &Q,
      DataWithFlagStatus<ValueTable<Value>> &phi) const;

    /**
     * Computes the 1st order derivative of the non-zero NURBS basis
     * functions over the current element,
     *   at the evaluation points pre-allocated in the cache.
     *
     * \warning If the output result @p derivatives_phi_hat is not correctly pre-allocated,
     * an exception will be raised.
     */
    void evaluate_nurbs_gradients_from_bspline(
      const BSplineElem &bspline_elem,
      const SafeSTLVector<Index> &bsp_local_to_patch,
      const ValueTable<typename BSplineElem::Value> &N,
      const ValueTable<typename BSplineElem::template Derivative<1>> &dN,
      const IgCoefficients &w_coeffs,
      const ValueVector<typename WeightFuncElem::Value> &Q,
      const ValueVector<typename WeightFuncElem::template Derivative<1>> &dQ,
      DataWithFlagStatus<ValueTable<Derivative<1>>> &D1_phi) const;

    /**
     * Computes the 2nd order derivative of the non-zero NURBS basis
     * functions over the current element,
     *   at the evaluation points pre-allocated in the cache.
     *
     * \warning If the output result @p derivatives_phi_hat is not correctly pre-allocated,
     * an exception will be raised.
     */
    void evaluate_nurbs_hessians_from_bspline(
      const BSplineElem &bspline_elem,
      const SafeSTLVector<Index> &bsp_local_to_patch,
      const ValueTable<typename BSplineElem::Value> &N,
      const ValueTable<typename BSplineElem::template Derivative<1>> &dN,
      const ValueTable<typename BSplineElem::template Derivative<2>> &d2N,
      const IgCoefficients &w_coeffs,
      const ValueVector<typename WeightFuncElem::Value> &Q,
      const ValueVector<typename WeightFuncElem::template Derivative<1>> &dQ,
      const ValueVector<typename WeightFuncElem::template Derivative<2>> &d2Q,
      DataWithFlagStatus<ValueTable<Derivative<2>>> &D2_phi) const;
  };

};




IGA_NAMESPACE_CLOSE

#endif // #ifdef USE_NURBS

#endif // #ifndef NURBS_ELEMENT_HANDLER_H_

