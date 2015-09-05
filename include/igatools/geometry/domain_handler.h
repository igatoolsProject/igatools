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

#ifndef __PHYSICAL_DOMAIN_CACHE_HANDLER_H_
#define __PHYSICAL_DOMAIN_CACHE_HANDLER_H_

#include <igatools/base/config.h>
#include <igatools/geometry/grid.h>
#include <igatools/geometry/domain.h>
#include <igatools/geometry/grid_handler.h>

IGA_NAMESPACE_OPEN


template <int, int, class> class DomainElementBase;
template <int, int> class DomainElement;
template <int, int> class ConstDomainElement;
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
 * @ingroup containers
 * @ingroup serializable
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
  using GridType = const Grid<dim_>;
  using GridHandler = typename GridType::ElementHandler;

  using ElementAccessor = DomainElement<dim_, codim_>;
  using ElementIterator = GridIterator<ElementAccessor>;
  using ConstElementAccessor = ConstDomainElement<dim_, codim_>;
  using ElementConstIterator = GridIterator<ConstElementAccessor>;

  using List = typename GridType::List;
  using ListIt = typename GridType::ListIt;
  using Flags = domain_element::Flags;
protected:
  using FlagsArray = SafeSTLArray<Flags, dim+1>;


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

  using topology_variant = typename FuncType::topology_variant;
  using eval_pts_variant = typename FuncType::eval_pts_variant;
#endif
  using topology_variant = TopologyVariants<dim_>;
  using eval_pts_variant = SubElemPtrVariants<Quadrature,dim_>;

private:
  /**
   * Default constructor. It does nothing but it is needed for the
   * <a href="http://www.boost.org/doc/libs/release/libs/serialization/">boost::serialization</a>
   * mechanism.
   */
  DomainHandler() = default;

public:
  DomainHandler(std::shared_ptr<DomainType> domain);


  ~DomainHandler();

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

  std::shared_ptr<DomainType> get_domain() const
  {
    return domain_;
  }


public:
//Is this really virtual?
  virtual void set_flags(const topology_variant &sdim,
                         const Flags &flag);
  template <int sdim>
  void set_flags(const Flags &flag)
  {
    this->set_flags(Topology<sdim>(), flag);
  }

  virtual void init_cache(ConstElementAccessor &elem,
                          const eval_pts_variant &quad) const;

  void init_cache(ElementConstIterator &elem,
                  const eval_pts_variant &quad) const
  {
    this->init_cache(*elem, quad);
  }

  virtual void fill_cache(const topology_variant &sdim,
                          ConstElementAccessor &elem,
                          const int s_id) const;

  void fill_cache(const topology_variant &sdim,
                  ElementConstIterator &elem,
                  const int s_id) const
  {
    this->fill_cache(sdim, *elem, s_id);
  }

  template <int sdim>
  void fill_cache(ElementConstIterator &elem,
                  const int s_id)
  {
    this->fill_cache(Topology<sdim>(), elem, s_id);
  }

protected:
  std::shared_ptr<typename ConstElementAccessor::CacheType>
  &get_element_cache(ConstElementAccessor &elem) const
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
    SetFlagsDispatcher(const Flags flag, FlagsArray &flags)
      :
      flag_(flag),
      flags_(flags)
    {}

    template<int sdim>
    void operator()(const Topology<sdim> &)
    {
      flags_[sdim] = flag_;
    }

    const Flags flag_;
    FlagsArray &flags_;
  };



  struct InitCacheDispatcher : boost::static_visitor<void>
  {
    InitCacheDispatcher(self_t const *domain_handler,
                        ConstElementAccessor &elem,
                        const FlagsArray &flags)
      :
      domain_handler_(domain_handler),
      elem_(elem),
      flags_(flags)
    {}


    template<int sdim>
    void operator()(const std::shared_ptr<const Quadrature<sdim>> &quad)
    {
      auto &cache = domain_handler_->get_element_cache(elem_);
      for (auto &s_id: UnitElement<dim_>::template elems_ids<sdim>())
      {
        auto &s_cache = cache->template get_sub_elem_cache<sdim>(s_id);
        const auto n_points = elem_.get_grid_element().template get_quad<sdim>()
                              ->get_num_points();
        s_cache.resize(flags_[sdim], n_points);
      }
    }

    self_t const *domain_handler_;
    ConstElementAccessor &elem_;
    const FlagsArray &flags_;
  };



  struct FillCacheDispatcher : boost::static_visitor<void>
  {
    FillCacheDispatcher(ConstElementAccessor &elem,
                        const int s_id)
      :
      elem_(elem),
      s_id_(s_id)
    {}


    template<int sdim>
    void operator()(const Topology<sdim> &)
    {
      using _Gradient = typename ElementAccessor::_Gradient;
      using _Measure = typename ElementAccessor::_Measure;

      const auto n_points = elem_.grid_elem_->template get_quad<sdim>()->get_num_points();

      auto &cache = elem_.local_cache_->template get_sub_elem_cache<sdim>(s_id_);

      if (cache.template status_fill<_Measure>())
      {
        auto &s_elem = UnitElement<dim_>::template get_elem<sdim>(s_id_);

        const auto &DF = cache.template get_data<_Gradient>();
        typename Domain<sdim, space_dim-sdim>::Gradient DF1;

        auto &measures = cache.template get_data<_Measure>();
        for (int pt = 0 ; pt < n_points; ++pt)
        {
          for (int l=0; l<sdim; ++l)
            DF1[l] = DF[pt][s_elem.active_directions[l]];

          measures[pt] = fabs(determinant<sdim,space_dim>(DF1));
        }
        cache.template set_status_filled<_Measure>(true);
      }


#if 0
      if (cache.template status_fill<_InvGradient>())
      {
        // TODO (pauletti, Nov 23, 2014): if also fill measure this could be done here
        const auto &DF = elem.template get_values<_Gradient, sdim>(s_id);
        auto &D_invF = cache.template get_data<_InvGradient>();
        Real det;
        for (int pt = 0 ; pt < n_points; ++pt)
          D_invF[pt] = inverse(DF[pt], det);

        cache.template set_status_filled<_InvGradient>(true);
      }

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

      if (cache.template status_fill<_BoundaryNormal>())
      {
        Assert(dim_ == sdim+1, ExcNotImplemented());
        const auto &D1_invF = cache.template get_data<_InvGradient>();
        const auto n_hat  = F_.get_grid()->template get_boundary_normals<sdim>(s_id)[0];
        auto &bndry_normal = cache.template get_data<_BoundaryNormal>();

        for (int pt = 0; pt < n_points; ++pt)
        {
          const auto D1_invF_t = co_tensor(transpose(D1_invF[pt]));
          bndry_normal[pt] = action(D1_invF_t, n_hat);
          bndry_normal[pt] /= bndry_normal[pt].norm();
        }

        cache.template set_status_filled<_BoundaryNormal>(true);
      }

      if (cache.template status_fill<_OuterNormal>())
      {
        Assert(sdim == dim_, ExcNotImplemented());
        Assert(codim_ == 1, ExcNotImplemented());

        const auto &DF = elem.template get_values<_Gradient, sdim>(s_id);
        auto &outer_normal = cache.template get_data<_OuterNormal>();

        for (int pt = 0; pt < n_points; ++pt)
        {
          outer_normal[pt] = cross_product<dim_, codim_>(DF[pt]);
          outer_normal[pt] /= outer_normal[pt].norm();
        }

        cache.template set_status_filled<_OuterNormal>(true);
      }


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

    ConstElementAccessor &elem_;
    const int s_id_;
  };



private:
  std::shared_ptr<DomainType> domain_;

  std::shared_ptr<GridHandler> grid_handler_;

  FlagsArray flags_;

  friend ElementAccessor;

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
    ar.template register_type<IgFunction<dim_,0,dim_+codim_,1> >();
    ar &boost::serialization::make_nvp("F_",F_);
    ar &boost::serialization::make_nvp("flags_",flags_);
  }
  ///@}
#endif
};

IGA_NAMESPACE_CLOSE

#endif

