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

#ifndef __FUNCTION_HANDLER_H_
#define __FUNCTION_HANDLER_H_

#include <igatools/base/config.h>
#include <igatools/base/tensor.h>
#include <igatools/geometry/unit_element.h>
#include <igatools/utils/value_vector.h>
#include <igatools/base/quadrature.h>
#include <igatools/geometry/grid_iterator.h>
#include <igatools/base/value_types.h>

IGA_NAMESPACE_OPEN

template <int,int> class Domain;
template <int,int> class DomainHandler;
template <int, int, int, int, class> class FunctionElementBase;
template <int, int, int, int> class FunctionElement;
template <int, int, int, int> class ConstFunctionElement;

/**
 * Function Class
 *
 * @ingroup serializable
 */
template<int dim_, int codim_ = 0, int range_ = 1, int rank_ = 1>
class FunctionElementHandler :
  public std::enable_shared_from_this<FunctionElementHandler<dim_,codim_,range_,rank_> >
{
private:
  using self_t = FunctionElementHandler<dim_, codim_, range_, rank_>;

public:
  using FuncType = const Function<dim_, codim_, range_, rank_>;

public:
  static const int space_dim = dim_ + codim_;
  static const int dim       = dim_;
  static const int codim     = codim_;
  static const int range     = range_;
  static const int rank      = rank_;

  using DomainType = Domain<dim_, codim_>;
  using DomainHandlerType = DomainHandler<dim_, codim_>;
  using ElementAccessor = typename FuncType::ElementAccessor;
  using ElementIterator = typename FuncType::ElementIterator;
  using ConstElementAccessor = typename FuncType::ConstElementAccessor;
  using ElementConstIterator = typename FuncType::ElementConstIterator;



  using List = typename DomainType::List;
  using ListIt = typename DomainType::ListIt;
  using Flags = function_element::Flags;
public:
  //using Flags = typename ElementAccessor::Flags;
protected:
  using FlagsArray = SafeSTLArray<Flags, dim+1>;

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

  using topology_variant = TopologyVariants<dim_>;

  template<int k>
  using ConstQuad = const Quadrature<k>;
  using eval_pts_variant = SubElemPtrVariants<ConstQuad,dim_>;

  /** @name Constructors and destructor. */
  ///@{
protected:
  /**
   * Default constructor. It does nothing but it is needed for the
   * <a href="http://www.boost.org/doc/libs/release/libs/serialization/">boost::serialization</a>
   * mechanism.
   */
  FunctionElementHandler() = default;

public:
  /** Constructor */
  FunctionElementHandler(std::shared_ptr<FuncType> func);

  /**
   * Copy constructor.
   */
  FunctionElementHandler(const self_t &func);

public:
  /** Destructor */
  virtual ~FunctionElementHandler() = default;
  ///@}



  //virtual std::shared_ptr<base_t> clone() const = 0;


  virtual void set_flags(const topology_variant &sdim,
                         const Flags &flag);

  virtual void init_cache(ConstElementAccessor &elem,
                          const eval_pts_variant &quad) const;

  void init_cache(ElementConstIterator &elem,
                  const eval_pts_variant &quad) const
  {
    this->init_cache(*elem, quad);
  }

//    void init_element_cache(ElementAccessor &elem) const;
//
//    void init_element_cache(ElementIterator &elem) const;

  virtual void fill_cache(const topology_variant &sdim,
                          ConstElementAccessor &elem,
                          const int s_id) const;

  void fill_cache(const topology_variant &sdim,
                  ElementConstIterator &elem,
                  const int s_id) const
  {
    this->fill_cache(sdim, *elem, s_id);
  }

//    void fill_element_cache(ElementAccessor &elem) const;
//    void fill_element_cache(ElementIterator &elem) const;
//    std::shared_ptr<typename ElementAccessor::CacheType>
//    &get_cache(ElementAccessor &elem);

  std::shared_ptr<const FuncType> get_function() const
  {
    return func_;
  }

protected:
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
    InitCacheDispatcher(const FlagsArray &flags,
                        ConstElementAccessor &elem)
      :
      flags_(flags),
      elem_(elem)
    {}

    template<int sdim>
    void operator()(const std::shared_ptr<const Quadrature<sdim>> &quad)
    {
      auto &cache = elem_.all_sub_elems_cache_;
      if (cache == nullptr)
      {
        using Cache = typename ConstElementAccessor::CacheType;
        cache = std::make_shared<Cache>();
      }

      const auto n_pts = quad->get_num_points();
      for (const auto s_id: UnitElement<dim_>::template elems_ids<sdim>())
      {
        auto &s_cache = cache->template get_sub_elem_cache<sdim>(s_id);
        s_cache.resize(flags_[sdim],n_pts);
      }
    }

    const FlagsArray &flags_;
    ConstElementAccessor &elem_;
  };


#ifdef MESH_REFINEMENT
private:
  std::shared_ptr<self_t> function_previous_refinement_;
public:
  const std::shared_ptr<self_t> &get_function_previous_refinement() const
  {
    return function_previous_refinement_;
  }
#endif // MESH_REFINEMENT

#ifdef SERIALIZATION
public:
  /**
   * Returns the unique identifier associated to each object instance.
   */
  Index get_object_id() const;

  /**
   * Get the name associated to the object instance.
   */
  const std::string &get_name() const;

  /**
   * Set the name associated to the object instance.
   */
  void set_name(const std::string &name);

private:
  /**
   * Unique identifier associated to each object instance.
   */
  Index object_id_;

  /**
   * Name associated to the object instance.
   */
  std::string name_;
public:
  /**
   * @name Functions needed for boost::serialization
   * @see <a href="http://www.boost.org/doc/libs/release/libs/serialization/">boost::serialization</a>
   */
  ///@{
  friend class boost::serialization::access;

  template<class Archive>
  void
  serialize(Archive &ar, const unsigned int version);
  /*
  {
      ar &boost::serialization::make_nvp("grid_elem_handler_",
                                         boost::serialization::base_object<Domain>(*this));

      ar &boost::serialization::make_nvp("flags_",flags_);
  }
  //*/
  ///@}
#endif // SERIALIZATION
};

IGA_NAMESPACE_CLOSE

#endif