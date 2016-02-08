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

#ifndef NEW_PUSH_FORWARD_ELEMENT_ACCESSOR_H_
#define NEW_PUSH_FORWARD_ELEMENT_ACCESSOR_H_

#include <igatools/geometry/domain_element.h>
#include <igatools/basis_functions/physical_basis_element.h>

IGA_NAMESPACE_OPEN


template <int,int,int,int>
class PhysicalBasisElement;

constexpr
int physical_range(const int ref_range, const int space_dim, const Transformation type)
{
  return type == Transformation::h_grad ? ref_range : space_dim;
}



/**
 *
 * @ingroup elements
 */
template<int dim_, int codim_ = 0>
class PushForward
{
private:
  using self_t  = PushForward<dim_, codim_>;
  using PhysDomainElem = DomainElement<dim_, codim_>;



public:

  PushForward(const Transformation &transformation_type)
    :
    transformation_type_(transformation_type)
  {};


  static const int dim = dim_;
  static const int codim = codim_;
  static const int space_dim = dim + codim;

#if 0
  template<int ref_range>
  struct PhysRange
  {
    static const int value = physical_range(ref_range, space_dim, type_);
  };

  template <int range, int rank>
  using RefValue = Values<dim, range, rank>;

  template <int range, int rank, int order>
  using RefDerivative = Derivatives<dim, range, rank, order>;

  template <int range, int rank>
  using PhysValue = Values<space_dim, PhysRange<range>::value,rank>;

  template <int range, int rank, int order>
  using PhysDerivative = Derivatives<space_dim, PhysRange<range>::value, rank, order>;
#endif

  template <int range,int rank>
  using RefBasisElem = ReferenceBasisElement<dim_,range,rank>;


  template <int range,int rank>
  using PhysBasisElem = PhysicalBasisElement<dim_,range,rank,codim>;

public:

  template <int range, int rank, int sdim>
  void
  transform_0(const int s_id,
              const RefBasisElem<range,rank> &ref_elem,
              const PhysDomainElem &phys_domain_elem,
              typename PhysBasisElem<range,rank>::Cache &phys_sub_elem_cache) const
  {
    if (transformation_type_ == Transformation::h_grad)
    {
      this->template hgrad_transform_0<range,rank,sdim>(
        s_id,ref_elem,phys_sub_elem_cache);
    }
    else if (transformation_type_ == Transformation::h_div)
    {
      Assert(false,ExcNotImplemented());
    }
    else if (transformation_type_ == Transformation::h_curl)
    {
      Assert(false,ExcNotImplemented());
    }
    else if (transformation_type_ == Transformation::l_2)
    {
      Assert(false,ExcNotImplemented());
    }
    else
    {
      Assert(false,ExcMessage("InvalidTransformationType"));
    }
  }

  template <int range, int rank, int sdim>
  void
  transform_1(const int s_id,
              const RefBasisElem<range,rank> &ref_elem,
              const PhysDomainElem &phys_domain_elem,
              typename PhysBasisElem<range,rank>::Cache &phys_sub_elem_cache) const
  {
    if (transformation_type_ == Transformation::h_grad)
    {
      this->template hgrad_transform_1<range,rank,sdim>(
        s_id,ref_elem,phys_domain_elem,phys_sub_elem_cache);
    }
    else if (transformation_type_ == Transformation::h_div)
    {
      Assert(false,ExcNotImplemented());
    }
    else if (transformation_type_ == Transformation::h_curl)
    {
      Assert(false,ExcNotImplemented());
    }
    else if (transformation_type_ == Transformation::l_2)
    {
      Assert(false,ExcNotImplemented());
    }
    else
    {
      Assert(false,ExcMessage("InvalidTransformationType"));
    }
  }

  template <int range, int rank, int sdim>
  void
  transform_2(const int s_id,
              const RefBasisElem<range,rank> &ref_elem,
              const PhysDomainElem &phys_domain_elem,
              typename PhysBasisElem<range,rank>::Cache &phys_sub_elem_cache) const
  {
    if (transformation_type_ == Transformation::h_grad)
    {
      this->template hgrad_transform_2<range,rank,sdim>(
        s_id,ref_elem,phys_domain_elem,phys_sub_elem_cache);
    }
    else if (transformation_type_ == Transformation::h_div)
    {
      Assert(false,ExcNotImplemented());
    }
    else if (transformation_type_ == Transformation::h_curl)
    {
      Assert(false,ExcNotImplemented());
    }
    else if (transformation_type_ == Transformation::l_2)
    {
      Assert(false,ExcNotImplemented());
    }
    else
    {
      Assert(false,ExcMessage("InvalidTransformationType"));
    }
  }

private:
  template <int range, int rank, int sdim>
  void
  hgrad_transform_0(const int s_id,
                    const RefBasisElem<range,rank> &ref_elem,
                    typename PhysBasisElem<range,rank>::Cache &phys_sub_elem_cache) const
  {
    using PhysElem = PhysBasisElem<range,rank>;
    using _Value = typename PhysElem::_Value;
    auto &v = phys_sub_elem_cache.template get_data<_Value>();
    const auto &v_hat = ref_elem.template get_basis_data<_Value,sdim>(s_id,DofProperties::active);

    v.fill(v_hat);
  }


  template <int range, int rank, int sdim>
  void
  hgrad_transform_1(const int s_id,
                    const RefBasisElem<range,rank> &ref_elem,
                    const PhysDomainElem &phys_domain_elem,
                    typename PhysBasisElem<range,rank>::Cache &phys_sub_elem_cache) const
  {
    using PhysElem = PhysBasisElem<range,rank>;
    using _Gradient = typename PhysElem::_Gradient;

    using _InvJacobian = typename PhysDomainElem::_InvJacobian;

    const auto &Dv_hat = ref_elem.template get_basis_data<_Gradient,sdim>(s_id,DofProperties::active);

    auto &Dv = phys_sub_elem_cache.template get_data<_Gradient>();

    const int n_func   = Dv_hat.get_num_functions();
    const int n_points = Dv_hat.get_num_points();
    auto Dv_it     = Dv.begin();
    auto Dv_hat_it = Dv_hat.cbegin();

    const auto &DF_inv = phys_domain_elem.template get_values_from_cache<_InvJacobian,sdim>(s_id);
    for (int fn = 0; fn < n_func; ++fn)
      for (int pt = 0; pt < n_points; ++pt, ++Dv_hat_it, ++Dv_it)
        (*Dv_it) = compose((*Dv_hat_it), DF_inv[pt]);

    Dv.set_status_filled(true);
  }


  template <int range, int rank, int sdim>
  void
  hgrad_transform_2(const int s_id,
                    const RefBasisElem<range,rank> &ref_elem,
                    const PhysDomainElem &phys_domain_elem,
                    typename PhysBasisElem<range,rank>::Cache &phys_sub_elem_cache) const
  {
    using PhysElem = PhysBasisElem<range,rank>;
    using _Gradient = typename PhysElem::_Gradient;
    using _Hessian = typename PhysElem::_Hessian;

    using _InvJacobian = typename PhysDomainElem::_InvJacobian;

    const auto &D2v_hat  = ref_elem.template get_basis_data< _Hessian,sdim>(s_id,DofProperties::active);

    const auto &D1v  = phys_sub_elem_cache.template get_data<_Gradient>();
    auto &D2v  = phys_sub_elem_cache.template get_data<_Hessian>();


    const int n_func   = D2v_hat.get_num_functions();
    const int n_points = D2v_hat.get_num_points();
    auto D2v_it     = D2v.begin();
    auto D1v_it     = D1v.cbegin();
    auto D2v_hat_it = D2v_hat.cbegin();
    const auto D2F = phys_domain_elem.get_grid_function_element().
                     template get_values_from_cache<grid_function_element::_D<2>,sdim>(s_id);
    const auto &DF_inv =
      phys_domain_elem.template get_values_from_cache<_InvJacobian,sdim>(s_id);

    for (int fn = 0; fn < n_func; ++fn)
      for (Index pt = 0; pt < n_points; ++pt)
      {
        const auto &D2F_pt = D2F[pt];
        const auto &DF_inv_pt = DF_inv[pt];
        for (int u = 0 ; u < dim ; ++u)
        {
          const auto &w = DF_inv_pt[u];
          (*D2v_it)[u] = compose(
                           action(*D2v_hat_it, w) - compose((*D1v_it),action(D2F_pt,w)),
                           DF_inv_pt);
        }
        ++D2v_hat_it;
        ++D1v_it;
        ++D2v_it;
      }

    D2v.set_status_filled(true);
  }

private:
  const Transformation transformation_type_;
//  const int phys_range_;
};

IGA_NAMESPACE_CLOSE

#endif
