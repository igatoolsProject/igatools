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
  using GridFuncHandler = typename GridFuncType::Handler;

  using ElementAccessor = DomainElement<dim_, codim_>;
  using ElementIterator = GridIterator<ElementAccessor>;

  using List = typename GridFuncType::List;
  using ListIt = typename GridFuncType::ListIt;
  using Flags = domain_element::Flags;

protected:
  using FlagsArray = SafeSTLArray<Flags, dim+1>;


  using topology_variant = TopologyVariants<dim_>;

  template<int k>
  using ConstQuad = const Quadrature<k>;
  using eval_pts_variant = SubElemPtrVariants<ConstQuad,dim_>;

private:

  DomainHandler() = delete;

public:
  DomainHandler(std::shared_ptr<DomainType> domain);


  ~DomainHandler() = default;


  std::shared_ptr<DomainType> get_domain() const;


public:
  void set_flags(const topology_variant &sdim,
                 const Flags &flag);

  template <int sdim>
  void set_flags(const Flags &flag);

  void set_element_flags(const Flags &flag);

  void init_cache(ElementAccessor &elem,
                  const eval_pts_variant &quad) const;

  void init_cache(ElementIterator &elem,
                  const eval_pts_variant &quad) const;

  void init_element_cache(ElementIterator &elem,
                          const std::shared_ptr<const Quadrature<dim_>> &quad) const;

  void init_element_cache(ElementAccessor &elem,
                          const std::shared_ptr<const Quadrature<dim_>> &quad) const;

  void fill_cache(const topology_variant &sdim,
                  ElementAccessor &elem,
                  const int s_id) const;

  template <int sdim>
  void fill_cache(ElementAccessor &elem,
                  const int s_id) const;

  template <int sdim>
  void fill_cache(ElementIterator &elem,
                  const int s_id) const;

  void fill_cache(const topology_variant &sdim,
                  ElementIterator &elem,
                  const int s_id) const;


  void fill_element_cache(ElementAccessor &elem) const;

  void fill_element_cache(ElementIterator &elem) const;


protected:
  typename ElementAccessor::CacheType
  &get_element_cache(ElementAccessor &elem) const;


private:
  /**
   * Alternative to
   * template <int sdim> set_flags()
   */
  struct SetFlagsDispatcher : boost::static_visitor<void>
  {
    SetFlagsDispatcher(const Flags cache_flag, FlagsArray &flags);

    template<int sdim>
    void operator()(const Topology<sdim> &);

    const Flags cache_flag_;
    FlagsArray &flags_;
  };



  struct InitCacheDispatcher : boost::static_visitor<void>
  {
    InitCacheDispatcher(const self_t &domain_handler,
                        ElementAccessor &elem,
                        const FlagsArray &flags);

    template<int sdim>
    void operator()(const std::shared_ptr<const Quadrature<sdim>> &quad);

    const self_t &domain_handler_;
    ElementAccessor &elem_;
    const FlagsArray &flags_;
  };



  struct FillCacheDispatcher : boost::static_visitor<void>
  {
    FillCacheDispatcher(ElementAccessor &elem,
                        const int s_id);


    template<int sdim>
    void operator()(const Topology<sdim> &);


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

