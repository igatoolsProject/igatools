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

#ifndef PHYS_SPACE_ELEMENT_HANDLER_H_
#define PHYS_SPACE_ELEMENT_HANDLER_H_

#include <igatools/base/config.h>
#include <igatools/basis_functions/bspline_element_handler.h>
#include <igatools/basis_functions/nurbs_element_handler.h>

IGA_NAMESPACE_OPEN




inline
auto
phys_space_to_reference_space_flag(
  const Transformation transformation_type,
  const typename space_element::Flags phys_space_flag)
-> typename space_element::Flags
{
  using space_element::Flags;

  Flags ref_flag = Flags::none;

  bool fill_values = false;
  bool fill_gradients = false;
  bool fill_hessians = false;
  bool fill_divergences = false;


  if (contains(phys_space_flag,Flags::value))
    fill_values = true;
  if (contains(phys_space_flag,Flags::gradient))
    fill_gradients = true;
  if (contains(phys_space_flag,Flags::hessian))
    fill_hessians = true;
  if (contains(phys_space_flag,Flags::divergence))
    fill_divergences = true;

  bool fill_D0_phi_hat = false;
  bool fill_D1_phi_hat = false;
  bool fill_D2_phi_hat = false;

  if (transformation_type == Transformation::h_grad)
  {
    fill_D0_phi_hat = fill_values;
    fill_D1_phi_hat = fill_gradients || fill_hessians || fill_divergences;
    fill_D2_phi_hat = fill_hessians;
  }
  else if (transformation_type == Transformation::h_div  ||
  transformation_type == Transformation::h_curl ||
  transformation_type == Transformation::l_2)
  {
    fill_D0_phi_hat = fill_values || fill_gradients || fill_hessians  || fill_divergences;
    fill_D1_phi_hat = fill_gradients || fill_hessians  || fill_divergences;
    fill_D2_phi_hat = fill_hessians;
  }

  if (fill_D0_phi_hat)
    ref_flag |= Flags::value;
  if (fill_D1_phi_hat)
    ref_flag |= Flags::gradient;
  if (fill_D2_phi_hat)
    ref_flag |= Flags::hessian;

  return ref_flag;
}


inline
domain_element::Flags
phys_space_to_domain_flag(
  const Transformation &transformation_type,
  const typename space_element::Flags phys_space_flag)
{
  /*
  ValueFlags transfer_flag =
    ValueFlags::measure |
    ValueFlags::w_measure |
    ValueFlags::outer_normal |
    ValueFlags::boundary_normal |
    ValueFlags::point;

  map_flags = ValueFlags::none;
  map_flags = flags & transfer_flag;
  //*/

  using SpaceFlags = space_element::Flags;

  using DomainFlags = domain_element::Flags;
  DomainFlags domain_flag = DomainFlags::none;

  if (contains(phys_space_flag, SpaceFlags::point))
    domain_flag |= DomainFlags::point;
  if (contains(phys_space_flag, SpaceFlags::w_measure))
    domain_flag |= DomainFlags::w_measure;


  if (transformation_type == Transformation::h_grad)
  {
    if (contains(phys_space_flag, SpaceFlags::value))
    {}

    if (contains(phys_space_flag, SpaceFlags::gradient))
      domain_flag |= (DomainFlags::inv_jacobian);

    if (contains(phys_space_flag, SpaceFlags::hessian))
      domain_flag |= (DomainFlags::hessian | DomainFlags::inv_jacobian);

//    if (contains(phys_space_flag, SpaceFlags::divergence))
//      AssertThrow(false,ExcNotImplemented());
//      domain_flag |= (DomainFlags::inv_gradient);

  }
  else
  {
    AssertThrow(false,ExcNotImplemented());
  }

  return domain_flag;
}





template<int dim, int range, int rank, int codim,Transformation type_>
class PhysicalSpace;

/**
 * Element handler for an isogeometric space
 *
 */
template<int dim_,int range_,int rank_,int codim_,Transformation type_>
class PhysSpaceElementHandler
  :
  public SpaceElementHandler<dim_,codim_,range_,rank_,type_>
{

  using PhysSpace = PhysicalSpace<dim_,range_,rank_,codim_,type_>;
  using RefSpace =  typename PhysSpace::RefSpace;
  using RefPhysSpaceElementHandler = typename PhysSpace::RefSpace::ElementHandler;
//    using PFCache = typename PhysSpace::PushForwardType;

  using ElementIterator = typename PhysSpace::ElementIterator;
  using ElementAccessor = typename PhysSpace::ElementAccessor;
  using PushFwd = typename PhysSpace::PushFwd;

  using base_t = SpaceElementHandler<dim_,codim_,range_,rank_,type_>;
  using self_t = PhysSpaceElementHandler<dim_,range_,rank_,codim_,type_>;

  using eval_pts_variant = QuadVariants<dim_>;
  using topology_variant = TopologyVariants<dim_>;

public:
  static const int dim = dim_;

//    using PhysSpace::PushForwardType::type;

  /**
   * @name Constructors
   */
  ///@{
#if 0
  /**
   * Default constructor. It does nothing but it is needed for the
   * <a href="http://www.boost.org/doc/libs/release/libs/serialization/">boost::serialization</a>
   * mechanism.
   */
  PhysSpaceElementHandler() = default;
#endif

protected:

  /**
   * Default constructor. Not allowed to be used.
   */
  PhysSpaceElementHandler() = delete;

  PhysSpaceElementHandler(std::shared_ptr<const PhysSpace> space);
  /**
   * Copy constructor. Not allowed to be used.
   */
  PhysSpaceElementHandler(const self_t &) = delete;

  /**
   * Move constructor. Not allowed to be used.
   */
  PhysSpaceElementHandler(self_t &&) = delete;

public:
  /**
   * Destructor.
   */
  virtual ~PhysSpaceElementHandler() = default;
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
   * @name Creators.
   */
  ///@{
  static std::unique_ptr<self_t> create(std::shared_ptr<const PhysSpace> space);
  ///@}

#if 0
  /**
   * Resets all the internal data in order to use the
   * quadrature scheme for the elements of the space with ID specified by
   * the input parameter <tt>elements_flat_id</tt>.
   */
  virtual void reset_selected_elements(
    const ValueFlags &flag,
    const eval_pts_variant &eval_points,
    const SafeSTLVector<int> &elements_flat_id) override final;


  virtual void init_cache(SpaceElement<dim_,codim_,range_,rank_,type_> &sp_elem,
                          const topology_variant &topology) override final;


  virtual void fill_cache(SpaceElement<dim_,codim_,range_,rank_,type_> &sp_elem,
                          const topology_variant &topology,
                          const int sub_elem_id) override final;
#endif

  void print_info(LogStream &out) const override final;

private:

  using RefElemHandler = SpaceElementHandler<RefSpace::dim,0,RefSpace::range,RefSpace::rank,Transformation::h_grad>;
  std::shared_ptr<RefElemHandler> ref_space_handler_;


  using PhysDomainHandler = typename PhysSpace::PhysDomain::ElementHandler;
  std::shared_ptr<PhysDomainHandler> phys_domain_handler_;




  virtual void set_flags_impl(
    const typename space_element::Flags &flag,
    const topology_variant &topology) override final
  {
    auto set_flag_dispatcher = SetFlagDispatcher(
                                 flag,
                                 *ref_space_handler_,
                                 *phys_domain_handler_,
                                 this->flags_);
    boost::apply_visitor(set_flag_dispatcher,topology);
  }

  struct SetFlagDispatcher : boost::static_visitor<void>
  {
    SetFlagDispatcher(const typename space_element::Flags phys_elem_flag,
                      RefElemHandler &ref_space_handler,
                      PhysDomainHandler &phys_domain_handler,
                      SafeSTLArray<typename space_element::Flags, dim+1> &flags)
      :
      phys_elem_flag_(phys_elem_flag),
      ref_space_handler_(ref_space_handler),
      phys_domain_handler_(phys_domain_handler),
      flags_(flags)
    {
    }

    template<int sdim>
    void operator()(const Topology<sdim> &topology)
    {
      ref_space_handler_.template set_flags<sdim>(
        phys_space_to_reference_space_flag(type_,phys_elem_flag_));

      phys_domain_handler_.template set_flags<sdim>(
        phys_space_to_domain_flag(type_,phys_elem_flag_));

      flags_[sdim] = phys_elem_flag_;

//      Assert(false,ExcNotImplemented());
    }


  private:
    const typename  space_element::Flags   phys_elem_flag_;
    RefElemHandler &ref_space_handler_;
    PhysDomainHandler &phys_domain_handler_;
    SafeSTLArray<typename space_element::Flags, dim+1> &flags_;
  };


  using BaseElem = SpaceElement<dim_,codim_,range_,rank_,type_>;

  virtual void init_cache_impl(BaseElem &elem,
                               const eval_pts_variant &quad) const override final
  {
    auto init_cache_dispatcher =
      InitCacheDispatcher(*ref_space_handler_,*phys_domain_handler_,this->flags_,elem);
    boost::apply_visitor(init_cache_dispatcher,quad);
  }

  struct InitCacheDispatcher : boost::static_visitor<void>
  {
    InitCacheDispatcher(const RefElemHandler &ref_space_handler,
                        const PhysDomainHandler &phys_domain_handler,
                        const SafeSTLArray<typename space_element::Flags, dim+1> &flags,
                        BaseElem &elem)
      :
      ref_space_handler_(ref_space_handler),
      phys_domain_handler_(phys_domain_handler),
      flags_(flags),
      elem_(elem)
    {}


    template<int sdim>
    void operator()(const std::shared_ptr<const Quadrature<sdim>> &quad)
    {
      using PhysSpaceElem = PhysicalSpaceElement<dim_,range_,rank_,codim_,type_>;
      auto &phys_space_elem  = dynamic_cast<PhysSpaceElem &>(elem_);

      ref_space_handler_.template init_cache<sdim>(
        *phys_space_elem.ref_space_element_,quad);

      phys_domain_handler_.init_cache(
        *phys_space_elem.phys_domain_element_,quad);


      auto &cache = phys_space_elem.all_sub_elems_cache_;

      const auto n_basis = phys_space_elem.get_num_basis(DofProperties::active);

      const auto n_pts = quad->get_num_points();

      const auto flag = flags_[sdim];

      for (auto &s_id: UnitElement<dim_>::template elems_ids<sdim>())
      {
        auto &s_cache = cache.template get_sub_elem_cache<sdim>(s_id);
        s_cache.resize(flag, n_pts, n_basis);
      }
    }

  private:
    const RefElemHandler &ref_space_handler_;
    const PhysDomainHandler &phys_domain_handler_;
    const SafeSTLArray<typename space_element::Flags, dim+1> &flags_;
    BaseElem &elem_;
  };


  virtual void fill_cache_impl(BaseElem &elem,
                               const topology_variant &topology,
                               const int s_id) const override final
  {
    auto fill_cache_dispatcher =
      FillCacheDispatcher(s_id,*ref_space_handler_,*phys_domain_handler_,elem);
    boost::apply_visitor(fill_cache_dispatcher,topology);
  }


  struct FillCacheDispatcher : boost::static_visitor<void>
  {
    FillCacheDispatcher(const int s_id,
                        const RefElemHandler &ref_space_handler,
                        const PhysDomainHandler &phys_domain_handler,
                        BaseElem &elem)
      :
      s_id_(s_id),
      ref_space_handler_(ref_space_handler),
      phys_domain_handler_(phys_domain_handler),
      elem_(elem)
    {}

    template<int sdim>
    void operator()(const Topology<sdim> &topology)
    {
      using PhysSpaceElem = PhysicalSpaceElement<dim_,range_,rank_,codim_,type_>;
      auto &phys_space_elem  = dynamic_cast<PhysSpaceElem &>(elem_);

      auto &ref_space_elem = *phys_space_elem.ref_space_element_;
      auto &phys_domain_elem = *phys_space_elem.phys_domain_element_;

      ref_space_handler_.template fill_cache<sdim>(ref_space_elem,s_id_);

      phys_domain_handler_.template fill_cache<sdim>(phys_domain_elem,s_id_);

      auto &all_sub_elems_cache = phys_space_elem.all_sub_elems_cache_;
      auto &sub_elem_cache = all_sub_elems_cache.template get_sub_elem_cache<sdim>(s_id_);

      using _Value = typename BaseElem::_Value;
      using _Gradient = typename BaseElem::_Gradient;
      using _Hessian = typename BaseElem::_Hessian;
      using _Divergence = typename BaseElem::_Divergence;


      if (sub_elem_cache.template status_fill<_Value>())
      {
        PushFwd::template
        transform_0<RefSpace::range,RefSpace::rank,sdim>(
          s_id_,
          ref_space_elem,
          phys_domain_elem,
          sub_elem_cache);
      }
      if (sub_elem_cache.template status_fill<_Gradient>())
      {
        PushFwd::template
        transform_1<RefSpace::range,RefSpace::rank,sdim>(
          s_id_,ref_space_elem,phys_domain_elem,sub_elem_cache);
      }
      if (sub_elem_cache.template status_fill<_Hessian>())
      {
        PushFwd::template
        transform_2<RefSpace::range,RefSpace::rank,sdim>(
          s_id_,ref_space_elem,phys_domain_elem,sub_elem_cache);
      }
      if (sub_elem_cache.template status_fill<_Divergence>())
      {
        auto &divergences = sub_elem_cache.template get_data<_Divergence>();
        eval_divergences_from_gradients(
          sub_elem_cache.template get_data<_Gradient>(),
          divergences);

        divergences.set_status_filled(true);
      }

      sub_elem_cache.set_filled(true);

    }

  private:

    const int s_id_;
    const RefElemHandler &ref_space_handler_;
    const PhysDomainHandler &phys_domain_handler_;
    BaseElem &elem_;
  };

#if 0
  struct ResetDispatcher : boost::static_visitor<void>
  {
    ResetDispatcher(
      const ValueFlags flag_in,
      const SafeSTLVector<Index> &elements_flat_id,
      RefElemHandler &ref_space_handler,
      Map &mapping,
      SafeSTLArray<ValueFlags, dim+1> &flags)
      :
      flag_in_(flag_in),
      elements_flat_id_(elements_flat_id),
      ref_space_handler_(ref_space_handler),
      phys_domain_handler_(mapping),
      flags_(flags)
    {};

    template<int sub_elem_dim>
    void operator()(const Quadrature<sub_elem_dim> &quad);

    const ValueFlags flag_in_;
    const SafeSTLVector<Index> &elements_flat_id_;
    RefElemHandler &ref_space_handler_;
    Map &phys_domain_handler_;
    SafeSTLArray<ValueFlags, dim+1> &flags_;
  };


  struct InitCacheDispatcher : boost::static_visitor<void>
  {
    InitCacheDispatcher(
      const SafeSTLArray<ValueFlags, dim+1> &flags,
      RefElemHandler &ref_space_handler,
      Map &mapping,
      PhysicalSpaceElement<dim_,range_,rank_,codim_,type_> &phys_elem)
      :
      flags_(flags),
      ref_space_handler_(ref_space_handler),
      phys_domain_handler_(mapping),
      phys_elem_(phys_elem)
    {};

    template<int sub_elem_dim>
    void operator()(const Topology<sub_elem_dim> &topology);

    const SafeSTLArray<ValueFlags, dim+1> &flags_;
    RefElemHandler &ref_space_handler_;
    Map &phys_domain_handler_;
    PhysicalSpaceElement<dim_,range_,rank_,codim_,type_> &phys_elem_;
  };



  struct FillCacheDispatcher : boost::static_visitor<void>
  {
    FillCacheDispatcher(
      const int sub_elem_id,
      RefElemHandler &ref_space_handler,
      Map &mapping,
      PhysicalSpaceElement<dim_,range_,rank_,codim_,type_> &phys_elem)
      :
      sub_elem_id_(sub_elem_id),
      ref_space_handler_(ref_space_handler),
      phys_domain_handler_(mapping),
      phys_elem_(phys_elem)
    {};

    template<int sub_elem_dim>
    void operator()(const Topology<sub_elem_dim> &topology);

    const int sub_elem_id_;
    RefElemHandler &ref_space_handler_;
    Map &phys_domain_handler_;
    PhysicalSpaceElement<dim_,range_,rank_,codim_,type_> &phys_elem_;
  };

#endif

#if 0
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
    ar &boost::serialization::make_nvp("PhysSpaceElementHandler_base_t",
                                       boost::serialization::base_object<base_t>(*this));

    ar.template register_type<BSplineElementHandler<dim_,range_,rank_> >();
#ifdef NURBS
    ar.template register_type<NURBSElementHandler<dim_,range_,rank_> >();
#endif // NURBS
    ar &boost::serialization::make_nvp("ref_space_handler_",ref_space_handler_);
    ar &boost::serialization::make_nvp("phys_domain_handler_",phys_domain_handler_);


//        ar &boost::serialization::make_nvp("push_fwd_",push_fwd_);


    ar &boost::serialization::make_nvp("flags_",flags_);
  }
  ///@}
#endif
#endif

};


IGA_NAMESPACE_CLOSE

#endif
