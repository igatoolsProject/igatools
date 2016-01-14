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

#ifndef BSPLINE_ELEMENT_HANDLER_H_
#define BSPLINE_ELEMENT_HANDLER_H_

#include <igatools/base/config.h>

//TODO(pauletti, Sep 9, 2014): should we instantiate the cartesian product instead
#include <igatools/utils/cartesian_product_array-template.h>

#include <igatools/basis_functions/reference_element_handler.h>
#include <igatools/basis_functions/bspline.h>
#include <igatools/basis_functions/bernstein_basis.h>
#include <igatools/basis_functions/values1d_const_view.h>


IGA_NAMESPACE_OPEN

/**
 *
 * @ingroup handlers
 */
template<int dim_, int range_, int rank_>
class BSplineElementHandler
  : public ReferenceElementHandler<dim_,range_,rank_>
{
  using base_t = ReferenceElementHandler<dim_,range_,rank_>;
  using self_t = BSplineElementHandler<dim_,range_,rank_>;
  using Basis = BSpline<dim_,range_,rank_>;
  static const Size n_components =  SplineSpace<dim_,range_,rank_>::n_components;

  using IndexType = typename Grid<dim_>::IndexType;

  template<class T>
  using ComponentContainer = typename Basis::template ComponentContainer<T>;


  template <int order>
  using Derivative = typename Basis::template Derivative<order>;

  using Value = typename Basis::Value;


protected:

  using BaseSpace = ReferenceBasis<dim_,range_,rank_>;
  using RefElementIterator = typename BaseSpace::ElementIterator;
  using RefElementAccessor = typename BaseSpace::ElementAccessor;


public:
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
   * @name Constructors.
   */
  ///@{

  BSplineElementHandler() = delete;

  BSplineElementHandler(std::shared_ptr<const Basis> space);

  /**
   * Copy constructor. Not allowed to be used.
   */
  BSplineElementHandler(const self_t &) = delete;

  /**
   * Move constructor. Not allowed to be used.
   */
  BSplineElementHandler(self_t &&) = delete;
  ///@}

public:
  static const int dim = dim_;

  /**
   * Destructor.
   */
  virtual ~BSplineElementHandler() = default;

//  static std::unique_ptr<self_t> create(std::shared_ptr<const Basis> space);

  using topology_variant = typename base_t::topology_variant;
  using eval_pts_variant = typename base_t::eval_pts_variant;




public:
  virtual void print_info(LogStream &out) const override final ;



  virtual void set_flags_impl(const topology_variant &topology,
                              const typename space_element::Flags &flag) override final;

private:
  struct SetFlagsDispatcher : boost::static_visitor<void>
  {
    SetFlagsDispatcher(const typename space_element::Flags flag_in,
                       GridHandler<dim_> &grid_handler,
                       SafeSTLArray<typename space_element::Flags, dim+1> &flags);

    template<int sdim>
    void operator()(const Topology<sdim> &topology);


  private:
    const typename space_element::Flags flag_in_;
    GridHandler<dim_> &grid_handler_;
    SafeSTLArray<typename space_element::Flags, dim+1> &flags_;
  };

public:
  using BaseElem = SpaceElement<dim_,0,range_,rank_>;
  using BSplineElem = BSplineElement<dim_,range_,rank_>;

  virtual void init_cache_impl(BaseElem &elem,
                               const eval_pts_variant &quad) const override final;

private:
  struct InitCacheDispatcher : boost::static_visitor<void>
  {
    InitCacheDispatcher(const GridHandler<dim_> &grid_handler,
                        const SafeSTLArray<typename space_element::Flags, dim+1> &flags,
                        BSplineElem &elem);


    template<int sdim>
    void operator()(const std::shared_ptr<const Quadrature<sdim>> &quad);

  private:
    const GridHandler<dim_> &grid_handler_;
    const SafeSTLArray<typename space_element::Flags, dim+1> &flags_;
    BSplineElem &bsp_elem_;

    template<int sdim>
    void init_cache_1D();

    template<int sdim>
    void init_cache_multiD();
  };

public:
  virtual void fill_cache_impl(const topology_variant &topology,
                               BaseElem &elem,
                               const int s_id) const override final;

private:
  struct FillCacheDispatcherNoGlobalCache : boost::static_visitor<void>
  {
    FillCacheDispatcherNoGlobalCache(const int s_id,
                                     const GridHandler<dim_> &grid_handler,
                                     BSplineElem &elem);

    template<int sdim>
    void operator()(const Topology<sdim> &topology);

  private:

    template<int sdim>
    void fill_cache_1D(const Quadrature<dim> &extended_sub_elem_quad);

    template<int sdim>
    void fill_cache_multiD(const Quadrature<dim> &extended_sub_elem_quad);

    /**
     * Computes the values (i.e. the 0-th order derivative) of the non-zero
     *  B-spline basis
     * functions over the current element,
     *   at the evaluation points pre-allocated in the cache.
     *
     * \warning If the output result @p D_phi is not correctly
     * pre-allocated,
     * an exception will be raised.
     */
    void evaluate_bspline_values(
      const ComponentContainer<std::unique_ptr<const TensorProductFunctionEvaluator<dim>>> &elem_values,
      ValueTable<Value> &D_phi) const;

    /**
     * Computes the k-th order derivative of the non-zero B-spline basis
     * functions over the current element,
     *   at the evaluation points pre-allocated in the cache.
     *
     * \warning If the output result @p D_phi is not correctly pre-allocated,
     * an exception will be raised.
     */
    template <int order>
    void evaluate_bspline_derivatives(
      const ComponentContainer<std::unique_ptr<const TensorProductFunctionEvaluator<dim>>> &elem_values,
      ValueTable<Derivative<order>> &D_phi) const;


    void
    copy_to_inactive_components_values(const SafeSTLVector<Index> &inactive_comp,
                                       const SafeSTLArray<Index, n_components> &active_map,
                                       ValueTable<Value> &D_phi) const;

    template <int order>
    void
    copy_to_inactive_components(const SafeSTLVector<Index> &inactive_comp,
                                const SafeSTLArray<Index, n_components> &active_map,
                                ValueTable<Derivative<order>> &D_phi) const;


    const int s_id_;
    const GridHandler<dim_> &grid_handler_;
    BSplineElem &bsp_elem_;
  };



#if 0
  /**
   * One-dimensional B-splines values and derivatives at quadrature points.
   * The values are stored with the following index ordering:
   *
   * splines1d_[comp][dir][interval][order][function][point]
   *
   * @ingroup serializable
   */
  class GlobalCache
  {
  private:


    /**
     * Quadrature points used for the 1D basis evaluation.
     */
    std::shared_ptr<const Quadrature<dim>> quad_;


    using BasisValues1dTable = ComponentContainer<SafeSTLArray<std::map<Index,BasisValues1d>,dim>>;

    /**
     * Values (and derivatives) of 1D basis precomputed in the initalized
     * interval of a given direction.
     *
     * @note The map's key is the interval id. In Debug mode, it will be
     * raised an assertion if
     * the requested values are not initialized for the interval.
     *
     */
    BasisValues1dTable basis_values_1d_table_;


  public:
    using ComponentMap = typename BasisValues1dTable::ComponentMap;

    GlobalCache() = default;

    GlobalCache(const std::shared_ptr<const Quadrature<dim>> &quad, const ComponentMap &component_map);


    ComponentContainer<std::unique_ptr<const TensorProductFunctionEvaluator<dim>> >
        get_element_values(const TensorIndex<dim> &elem_tensor_id) const;

    BasisValues1d &entry(const int comp, const int dir, const Index interval_id);

    void print_info(LogStream &out) const;

  private:
#ifdef SERIALIZATION
    /**
     * @name Functions needed for boost::serialization
     * @see <a href="http://www.boost.org/doc/libs/release/libs/serialization/">boost::serialization</a>
     */
    ///@{
    friend class boost::serialization::access;

    template<class Archive>
    void
    serialize(Archive &ar, const unsigned int version)
    {
      ar &boost::serialization::make_nvp("quad_",quad_);
      ar &boost::serialization::make_nvp("basis_values_1d_table_",basis_values_1d_table_);
    }
    ///@}
#endif // SERIALIZATION
  };
#endif



  /**
   * Returns the BSpline used to define the BSplineElementHandler object.
   */
  std::shared_ptr<const Basis> get_bspline_basis() const;
};

IGA_NAMESPACE_CLOSE


#endif
