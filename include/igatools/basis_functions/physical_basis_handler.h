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

#ifndef __PHYS_BASIS_HANDLER_H_
#define __PHYS_BASIS_HANDLER_H_

#include <igatools/base/config.h>
#include <igatools/basis_functions/bspline_handler.h>
#include <igatools/basis_functions/nurbs_handler.h>

IGA_NAMESPACE_OPEN




inline
auto
phys_basis_to_reference_basis_flag(
  const Transformation transformation_type,
  const typename basis_element::Flags phys_basis_flag)
-> typename basis_element::Flags
{
  using basis_element::Flags;

  Flags ref_flag = Flags::none;

  bool fill_values = false;
  bool fill_gradients = false;
  bool fill_hessians = false;
  bool fill_divergences = false;


  if (contains(phys_basis_flag,Flags::value))
    fill_values = true;
  if (contains(phys_basis_flag,Flags::gradient))
    fill_gradients = true;
  if (contains(phys_basis_flag,Flags::hessian))
    fill_hessians = true;
  if (contains(phys_basis_flag,Flags::divergence))
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
phys_basis_to_domain_flag(
  const Transformation &transformation_type,
  const typename basis_element::Flags phys_basis_flag)
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

  using BasisFlags = basis_element::Flags;

  using DomainFlags = domain_element::Flags;
  DomainFlags domain_flag = DomainFlags::none;

  if (contains(phys_basis_flag, BasisFlags::point))
    domain_flag |= DomainFlags::point;
  if (contains(phys_basis_flag, BasisFlags::w_measure))
    domain_flag |= DomainFlags::w_measure;


  if (transformation_type == Transformation::h_grad)
  {
    if (contains(phys_basis_flag, BasisFlags::value))
    {}

    if (contains(phys_basis_flag, BasisFlags::gradient))
      domain_flag |= (DomainFlags::inv_jacobian);

    if (contains(phys_basis_flag, BasisFlags::hessian))
      domain_flag |= (DomainFlags::hessian | DomainFlags::inv_jacobian);

//    if (contains(phys_basis_flag, BasisFlags::divergence))
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
 * Element handler for an isogeometric basis
 *
 * @ingroup handlers
 */
template<int dim_,int range_,int rank_,int codim_>
class PhysicalBasisHandler
  :
  public BasisHandler<dim_,codim_,range_,rank_>
{

  using PhysBasis = PhysicalBasis<dim_,range_,rank_,codim_>;
  using RefBasis =  typename PhysBasis::RefBasis;
  using RefPhysicalBasisHandler = typename PhysBasis::RefBasis::Handler;
//    using PFCache = typename PhysBasis::PushForwardType;

  using ElementIterator = typename PhysBasis::ElementIterator;
  using ElementAccessor = typename PhysBasis::ElementAccessor;

  using base_t = BasisHandler<dim_,codim_,range_,rank_>;
  using self_t = PhysicalBasisHandler<dim_,range_,rank_,codim_>;

  using eval_pts_variant = QuadVariants<dim_>;
  using topology_variant = TopologyVariants<dim_>;

public:
  static const int dim = dim_;

//    using PhysBasis::PushForwardType::type;

  /**
   * @name Constructors
   */
  ///@{

public:

  /**
   * Default constructor. Not allowed to be used.
   */
  PhysicalBasisHandler() = delete;

  PhysicalBasisHandler(std::shared_ptr<const PhysBasis> basis);
  /**
   * Copy constructor. Not allowed to be used.
   */
  PhysicalBasisHandler(const self_t &) = delete;

  /**
   * Move constructor. Not allowed to be used.
   */
  PhysicalBasisHandler(self_t &&) = delete;

public:
  /**
   * Destructor.
   */
  virtual ~PhysicalBasisHandler() = default;
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

  using RefElemHandler = BasisHandler<RefBasis::dim,0,RefBasis::range,RefBasis::rank>;
  std::unique_ptr<RefElemHandler> ref_basis_handler_;


  using PhysDomainHandler = typename PhysBasis::PhysDomain::Handler;
  std::unique_ptr<PhysDomainHandler> phys_domain_handler_;




  virtual void set_flags_impl(const topology_variant &topology,
                              const typename basis_element::Flags &flag) override final;

  struct SetFlagsDispatcher : boost::static_visitor<void>
  {
    SetFlagsDispatcher(const typename basis_element::Flags phys_elem_flag,
                       const Transformation &transformation_type,
                       RefElemHandler &ref_basis_handler,
                       PhysDomainHandler &phys_domain_handler,
                       SafeSTLArray<typename basis_element::Flags, dim+1> &flags);

    template<int sdim>
    void operator()(const Topology<sdim> &topology);


  private:
    const typename  basis_element::Flags   phys_elem_flag_;
    const Transformation transformation_type_;
    RefElemHandler &ref_basis_handler_;
    PhysDomainHandler &phys_domain_handler_;
    SafeSTLArray<typename basis_element::Flags, dim+1> &flags_;
  };


  using BaseElem = BasisElement<dim_,codim_,range_,rank_>;

  virtual void init_cache_impl(BaseElem &elem,
                               const eval_pts_variant &quad) const override final;

  struct InitCacheDispatcher : boost::static_visitor<void>
  {
    InitCacheDispatcher(const RefElemHandler &ref_basis_handler,
                        const PhysDomainHandler &phys_domain_handler,
                        const SafeSTLArray<typename basis_element::Flags, dim+1> &flags,
                        BaseElem &elem);


    template<int sdim>
    void operator()(const std::shared_ptr<const Quadrature<sdim>> &quad);

  private:
    const RefElemHandler &ref_basis_handler_;
    const PhysDomainHandler &phys_domain_handler_;
    const SafeSTLArray<typename basis_element::Flags, dim+1> &flags_;
    BaseElem &elem_;
  };


  virtual void fill_cache_impl(const topology_variant &topology,
                               BaseElem &elem,
                               const int s_id) const override final;


  struct FillCacheDispatcher : boost::static_visitor<void>
  {
    FillCacheDispatcher(const int s_id,
                        const RefElemHandler &ref_basis_handler,
                        const PhysDomainHandler &phys_domain_handler,
                        const self_t &phys_basis_handler,
                        BaseElem &elem);

    template<int sdim>
    void operator()(const Topology<sdim> &topology);

  private:

    const int s_id_;
    const RefElemHandler &ref_basis_handler_;
    const PhysDomainHandler &phys_domain_handler_;
    const self_t &phys_basis_handler_;
    BaseElem &elem_;
  };


  std::shared_ptr<const PhysBasis> phys_basis_;
};


IGA_NAMESPACE_CLOSE

#endif // __PHYS_BASIS_HANDLER_H_
