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





template<int dim, int range, int rank, int codim>
class PhysicalBasis;

/**
 * Element handler for an isogeometric space
 *
 * @ingroup handlers
 */
template<int dim_,int range_,int rank_,int codim_>
class PhysicalBasisElementHandler
  :
  public BasisElementHandler<dim_,codim_,range_,rank_>
{

  using PhysSpace = PhysicalBasis<dim_,range_,rank_,codim_>;
  using RefBasis =  typename PhysSpace::RefBasis;
  using RefPhysicalBasisElementHandler = typename PhysSpace::RefBasis::ElementHandler;
//    using PFCache = typename PhysSpace::PushForwardType;

  using ElementIterator = typename PhysSpace::ElementIterator;
  using ElementAccessor = typename PhysSpace::ElementAccessor;

  using base_t = BasisElementHandler<dim_,codim_,range_,rank_>;
  using self_t = PhysicalBasisElementHandler<dim_,range_,rank_,codim_>;

  using eval_pts_variant = QuadVariants<dim_>;
  using topology_variant = TopologyVariants<dim_>;

public:
  static const int dim = dim_;

//    using PhysSpace::PushForwardType::type;

  /**
   * @name Constructors
   */
  ///@{

public:

  /**
   * Default constructor. Not allowed to be used.
   */
  PhysicalBasisElementHandler() = delete;

  PhysicalBasisElementHandler(std::shared_ptr<const PhysSpace> space);
  /**
   * Copy constructor. Not allowed to be used.
   */
  PhysicalBasisElementHandler(const self_t &) = delete;

  /**
   * Move constructor. Not allowed to be used.
   */
  PhysicalBasisElementHandler(self_t &&) = delete;

public:
  /**
   * Destructor.
   */
  virtual ~PhysicalBasisElementHandler() = default;
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


  void print_info(LogStream &out) const override final;

private:

  using RefElemHandler = BasisElementHandler<RefBasis::dim,0,RefBasis::range,RefBasis::rank>;
  std::unique_ptr<RefElemHandler> ref_space_handler_;


  using PhysDomainHandler = typename PhysSpace::PhysDomain::ElementHandler;
  std::unique_ptr<PhysDomainHandler> phys_domain_handler_;




  virtual void set_flags_impl(const topology_variant &topology,
                              const typename space_element::Flags &flag) override final;

  struct SetFlagsDispatcher : boost::static_visitor<void>
  {
    SetFlagsDispatcher(const typename space_element::Flags phys_elem_flag,
                       const Transformation &transformation_type,
                       RefElemHandler &ref_space_handler,
                       PhysDomainHandler &phys_domain_handler,
                       SafeSTLArray<typename space_element::Flags, dim+1> &flags);

    template<int sdim>
    void operator()(const Topology<sdim> &topology);


  private:
    const typename  space_element::Flags   phys_elem_flag_;
    const Transformation transformation_type_;
    RefElemHandler &ref_space_handler_;
    PhysDomainHandler &phys_domain_handler_;
    SafeSTLArray<typename space_element::Flags, dim+1> &flags_;
  };


  using BaseElem = BasisElement<dim_,codim_,range_,rank_>;

  virtual void init_cache_impl(BaseElem &elem,
                               const eval_pts_variant &quad) const override final;

  struct InitCacheDispatcher : boost::static_visitor<void>
  {
    InitCacheDispatcher(const RefElemHandler &ref_space_handler,
                        const PhysDomainHandler &phys_domain_handler,
                        const SafeSTLArray<typename space_element::Flags, dim+1> &flags,
                        BaseElem &elem);


    template<int sdim>
    void operator()(const std::shared_ptr<const Quadrature<sdim>> &quad);

  private:
    const RefElemHandler &ref_space_handler_;
    const PhysDomainHandler &phys_domain_handler_;
    const SafeSTLArray<typename space_element::Flags, dim+1> &flags_;
    BaseElem &elem_;
  };


  virtual void fill_cache_impl(const topology_variant &topology,
                               BaseElem &elem,
                               const int s_id) const override final;


  struct FillCacheDispatcher : boost::static_visitor<void>
  {
    FillCacheDispatcher(const int s_id,
                        const RefElemHandler &ref_space_handler,
                        const PhysDomainHandler &phys_domain_handler,
                        const self_t &phys_space_handler,
                        BaseElem &elem);

    template<int sdim>
    void operator()(const Topology<sdim> &topology);

  private:

    const int s_id_;
    const RefElemHandler &ref_space_handler_;
    const PhysDomainHandler &phys_domain_handler_;
    const self_t &phys_space_handler_;
    BaseElem &elem_;
  };


  std::shared_ptr<const PhysSpace> phys_space_;
};


IGA_NAMESPACE_CLOSE

#endif
