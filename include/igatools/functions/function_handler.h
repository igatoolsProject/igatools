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

#ifndef __FUNCTION_HANDLER_H_
#define __FUNCTION_HANDLER_H_

#include <igatools/base/config.h>
#include <igatools/geometry/domain.h>
#include <igatools/functions/function.h>


IGA_NAMESPACE_OPEN

template <int,int> class DomainHandler;

template <int, int, int, int> class FunctionElement;

/**
 * Function Class
 *
 */
template<int dim_, int codim_ = 0, int range_ = 1, int rank_ = 1>
class FunctionHandler :
  public std::enable_shared_from_this<FunctionHandler<dim_,codim_,range_,rank_> >
{
private:
  using self_t = FunctionHandler<dim_, codim_, range_, rank_>;

public:
  static const int space_dim = dim_ + codim_;
  static const int dim       = dim_;
  static const int codim     = codim_;
  static const int range     = range_;
  static const int rank      = rank_;

  using FuncType = const Function<dim_, codim_, range_, rank_>;
  using DomainType = const Domain<dim_, codim_>;
  using DomainHandlerType = typename DomainType::Handler;

  using ElementAccessor = typename FuncType::ElementAccessor;
  using ElementIterator = typename FuncType::ElementIterator;

  using List = typename DomainType::List;
  using ListIt = typename DomainType::ListIt;
  using Flags = function_element::Flags;

protected:
  using FlagsArray = SafeSTLArray<Flags, dim+1>;

#if 0
  /** Types for the input/output evaluation arguments */
  ///@{
  /**
   * Type for the input argument of the function.
   */
  using Point = Points<space_dim>;

  /**
   * Type for the return of the function.
   */
  using Value = Values<space_dim, range_, rank_>;

  /**
   * Type for the derivative of the function.
   */
  template <int order>
  using Derivative = Derivatives<space_dim, range_, rank_, order>;

  /**
   * Type for the gradient of the function.
   */
  using Gradient = Derivative<1>;

  /**
   * Type for the hessian of the function.
   */
  using Hessian = Derivative<2>;

  /**
   * Type for the divergence of function.
   */
  using Div = Values<space_dim, space_dim, rank_-1>;
  ///@}
#endif

  using topology_variant = TopologyVariants<dim_>;

  template<int k>
  using ConstQuad = const Quadrature<k>;
  using eval_pts_variant = SubElemPtrVariants<ConstQuad,dim_>;

  /** @name Constructors and destructor. */
  ///@{
protected:
  /**
   * Default constructor. Not allowed to be used.
   */
  FunctionHandler() = delete;

public:
  /** Constructor */
  FunctionHandler(std::shared_ptr<FuncType> func);

//  /**
//   * Copy constructor.
//   */
//  FunctionHandler(const self_t &func);

public:
  /** Destructor */
  virtual ~FunctionHandler() = default;
  ///@}


  static std::shared_ptr<self_t>
  create(std::shared_ptr<FuncType> func)
  {
    return std::shared_ptr<self_t>(new self_t(func));
  }


  static std::shared_ptr<const self_t>
  const_create(std::shared_ptr<FuncType> func)
  {
    return create(func);
  }


  std::shared_ptr<const FuncType> get_function() const
  {
    return func_;
  }

public:
  std::shared_ptr<const DomainHandlerType> get_domain_handler() const
  {
    return domain_handler_;
  }

protected:
  std::shared_ptr<DomainHandlerType> get_domain_handler()
  {
    return domain_handler_;
  }
public:
  virtual void set_flags(const topology_variant &sdim,
                         const Flags &flag);

  template <int sdim>
  void set_flags(const Flags &flag)
  {
    this->set_flags(Topology<sdim>(), flag);
  }

  void set_element_flags(const Flags &flag);

  virtual void init_cache(ElementAccessor &elem,
                          const eval_pts_variant &quad) const;

  void init_cache(ElementIterator &elem,
                  const eval_pts_variant &quad) const;


  virtual void fill_cache(const topology_variant &sdim,
                          ElementAccessor &elem,
                          const int s_id) const;

  void fill_cache(const topology_variant &sdim,
                  ElementIterator &elem,
                  const int s_id) const;

  template <int sdim>
  void fill_cache(ElementIterator &elem,
                  const int s_id)
  {
    this->fill_cache(Topology<sdim>(), elem, s_id);
  }

  template <int sdim>
  void fill_cache(ElementAccessor &elem,
                  const int s_id)
  {
    this->fill_cache(Topology<sdim>(), elem, s_id);
  }

  void fill_element_cache(ElementAccessor &elem);

  void fill_element_cache(ElementIterator &elem);

protected:
//  std::shared_ptr<typename ElementAccessor::CacheType>
  typename ElementAccessor::CacheType &
  get_element_cache(ElementAccessor &elem) const
  {
    return  elem.local_cache_;
  }

private:
  std::shared_ptr<FuncType> func_;

  std::shared_ptr<DomainHandlerType> domain_handler_;

  FlagsArray flags_;

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


  /**
   * Alternative to
   * template <int sdim> init_cache()
   */
  struct InitCacheDispatcher : boost::static_visitor<void>
  {
    InitCacheDispatcher(ElementAccessor &elem,
                        const FlagsArray &flags)
      :
      elem_(elem),
      flags_(flags)
    {}

    template<int sdim>
    void operator()(const std::shared_ptr<const Quadrature<sdim>> &quad)
    {
      auto &cache = elem_.local_cache_;

      const auto n_pts = elem_.
                         get_domain_element().
                         get_grid_function_element().
                         get_grid_element().
                         template get_quad<sdim>()->get_num_points();
      for (const auto s_id: UnitElement<dim_>::template elems_ids<sdim>())
      {
        auto &s_cache = cache.template get_sub_elem_cache<sdim>(s_id);
        s_cache.resize(flags_[sdim],n_pts);
      }
    }

    ElementAccessor &elem_;
    const FlagsArray &flags_;
  };


};

IGA_NAMESPACE_CLOSE

#endif
