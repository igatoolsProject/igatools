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

#ifndef __DOMAIN_HANDLER_H_
#define __DOMAIN_HANDLER_H_

#include <igatools/base/config.h>
#include <igatools/geometry/grid.h>
#include <igatools/geometry/domain.h>
#include <igatools/geometry/grid_handler.h>

IGA_NAMESPACE_OPEN


template <int, int> class DomainElement;

/**
 * @brief The mapping is a deformation \f$ F : \hat\Omega \to \Omega\f$
 * which maps the reference domain \f$\hat\Omega \in \mathbb{R}^{dim}\f$ to the
 * physical domain \f$\Omega \in \mathbb{R}^{dim+codim}\f$.
 *
 * Domain is the physical domain, wether of a function or a space.
 *
 * It is a function with special properties: it codim is 0 and the map is always
 * the identity.
 *
 * @todo we should thing about renaming mapping to physical domain
 *
 * @ingroup handlers
 *
 * @author pauletti 2014, 2015
 * @author M. Martinelli, 2015
 */
template<int dim_, int codim_ = 0>
class DomainHandler :
  public std::enable_shared_from_this<DomainHandler<dim_,codim_> >
{
private:
  using self_t = DomainHandler<dim_, codim_>;

public:
  static const int space_dim = dim_ + codim_;
  static const int dim = dim_;

  using DomainType = const Domain<dim_, codim_>;
  using GridFuncType =  typename DomainType::GridFuncType;
  using GridFuncHandler = typename GridFuncType::ElementHandler;

  using ElementAccessor = DomainElement<dim_, codim_>;
  using ElementIterator = GridIterator<ElementAccessor>;

  using List = typename GridFuncType::List;
  using ListIt = typename GridFuncType::ListIt;
  using Flags = domain_element::Flags;
  using CacheFlags = domain_element::CacheFlags;
protected:
  using FlagsArray = SafeSTLArray<CacheFlags, dim+1>;


#if 0
  /** Type for the given order derivatives of the
   *  the mapping. */
  template<int order>
  using Derivative = typename FuncType::template Derivative<order>;

  /** Type for the diferent order derivatives of the inverse of
   * the mapping
   */
  template<int order>
  using InvDerivative = Derivatives<space_dim, dim_, 1, order>;

  /** Type of the mapping evaluation point. */
  using Point = typename FuncType::Point;
  /** Type of the mapping return value. */
  using Value = typename FuncType::Value;
  /** Type of the mapping gradient. */
  using Gradient = typename FuncType::Gradient;

  /** Typedef for the mapping hessian. */
  using Hessian = typename FuncType::Hessian;

#endif

  using topology_variant = TopologyVariants<dim_>;

  template<int k>
  using ConstQuad = const Quadrature<k>;
  using eval_pts_variant = SubElemPtrVariants<ConstQuad,dim_>;

private:

  DomainHandler() = delete;

public:
  DomainHandler(std::shared_ptr<DomainType> domain);


  ~DomainHandler() = default;

#if 0
  static std::shared_ptr<self_t>
  create(std::shared_ptr<DomainType> domain)
  {
    return std::shared_ptr<self_t>(new self_t(domain));
  }


  static std::shared_ptr<const self_t>
  const_create(std::shared_ptr<DomainType> domain)
  {
    return create(domain);
  }
#endif

  std::shared_ptr<DomainType> get_domain() const;


public:
  void set_flags(const topology_variant &sdim,
                 const Flags &flag);

  template <int sdim>
  void set_flags(const Flags &flag)
  {
    this->set_flags(Topology<sdim>(), flag);
  }

  void set_element_flags(const Flags &flag);

  void init_cache(ElementAccessor &elem,
                  const eval_pts_variant &quad) const;

  void init_cache(ElementIterator &elem,
                  const eval_pts_variant &quad) const;

  void init_element_cache(ElementIterator &elem,
                          const std::shared_ptr<Quadrature<dim_>> &quad) const;

  void init_element_cache(ElementAccessor &elem,
                          const std::shared_ptr<Quadrature<dim_>> &quad) const;

  void fill_cache(const topology_variant &sdim,
                  ElementAccessor &elem,
                  const int s_id) const;

  template <int sdim>
  void fill_cache(ElementAccessor &elem,
                  const int s_id) const
  {
    this->fill_cache(Topology<sdim>(), elem, s_id);
  }

  void fill_cache(const topology_variant &sdim,
                  ElementIterator &elem,
                  const int s_id) const;

  template <int sdim>
  void fill_cache(ElementIterator &elem,
                  const int s_id) const
  {
    this->fill_cache(Topology<sdim>(), elem, s_id);
  }

  void fill_element_cache(ElementAccessor &elem) const;

  void fill_element_cache(ElementIterator &elem) const;


protected:
//  std::shared_ptr<typename ElementAccessor::CacheType>
  typename ElementAccessor::CacheType
  &get_element_cache(ElementAccessor &elem) const
  {
    return  elem.local_cache_;
  }


private:
  /**
   * Alternative to
   * template <int sdim> set_flags()
   */
  struct SetFlagsDispatcher : boost::static_visitor<void>
  {
    SetFlagsDispatcher(const CacheFlags cache_flag, FlagsArray &flags)
      :
      cache_flag_(cache_flag),
      flags_(flags)
    {}

    template<int sdim>
    void operator()(const Topology<sdim> &)
    {
      flags_[sdim] |= cache_flag_;
    }

    const CacheFlags cache_flag_;
    FlagsArray &flags_;
  };



  struct InitCacheDispatcher : boost::static_visitor<void>
  {
    InitCacheDispatcher(const self_t &domain_handler,
                        ElementAccessor &elem,
                        const FlagsArray &flags)
      :
      domain_handler_(domain_handler),
      elem_(elem),
      flags_(flags)
    {}


    template<int sdim>
    void operator()(const std::shared_ptr<const Quadrature<sdim>> &quad)
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

    const self_t &domain_handler_;
    ElementAccessor &elem_;
    const FlagsArray &flags_;
  };



  struct FillCacheDispatcher : boost::static_visitor<void>
  {
    FillCacheDispatcher(ElementAccessor &elem,
                        const int s_id)
      :
      elem_(elem),
      s_id_(s_id)
    {}


    template<int sdim>
    void operator()(const Topology<sdim> &)
    {
      auto &cache = elem_.local_cache_.template get_sub_elem_cache<sdim>(s_id_);

      using _Measure = typename ElementAccessor::_Measure;
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


      using _InvJacobian = typename ElementAccessor::_InvJacobian;
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

#if 0
      if (cache.template status_fill<_InvHessian>())
      {
        //        const auto &D1_F = elem.template get_values<_Gradient, sdim>(j);
        const auto &D2_F = elem.template get_values<_Hessian, sdim>(s_id);
        const auto &D1_invF = cache.template get_data<_InvGradient>();
        auto &D2_invF       = cache.template get_data<_InvHessian>();

        for (int pt = 0 ; pt < n_points; ++pt)
          for (int u=0; u<dim_; ++u)
          {
            const auto tmp_u = action(D2_F[pt], D1_invF[pt][u]);
            for (int v=0; v<dim_; ++v)
            {
              const auto tmp_u_v = action(tmp_u, D1_invF[pt][v]);
              D2_invF[pt][u][v] = - action(D1_invF[pt], tmp_u_v);
            }
          }

        cache.template set_status_filled<_InvHessian>(true);
      }
#endif

      using _BoundaryNormal = typename ElementAccessor::_BoundaryNormal;
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

      using _ExtNormal = typename ElementAccessor::_ExtNormal;
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


#if 0
      if (cache.template status_fill<_Curvature>())
      {
        Assert(sdim == dim_, ExcNotImplemented());
        Assert(codim_ == 1, ExcNotImplemented());

        const auto H = elem.compute_second_fundamental_form();
        const auto G_inv = elem.compute_inv_first_fundamental_form();

        auto &curvatures = cache.template get_data<_Curvature>();

        for (int pt = 0; pt < n_points; ++pt)
        {
          //          const MetricTensor B = compose(H[pt], G_inv[pt]);
          const auto B = compose(H[pt], G_inv[pt]);
          const auto A = unroll_to_matrix(B);
          curvatures[pt] = A.eigen_values();
        }

        cache.template set_status_filled<_Curvature>(true);
      }
#endif
      cache.set_filled(true);
    }

    ElementAccessor &elem_;
    const int s_id_;
  };


private:
  std::shared_ptr<DomainType> domain_;

  std::unique_ptr<GridFuncHandler> grid_func_handler_;

  FlagsArray flags_;

//  friend ElementAccessor;
};

IGA_NAMESPACE_CLOSE

#endif

