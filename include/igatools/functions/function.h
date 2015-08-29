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

#ifndef __FUNCTION_H_
#define __FUNCTION_H_

#include <igatools/base/config.h>
#include <igatools/base/tensor.h>
#include <igatools/geometry/unit_element.h>
#include <igatools/utils/value_vector.h>
#include <igatools/base/quadrature.h>
#include <igatools/geometry/grid_iterator.h>

IGA_NAMESPACE_OPEN


namespace function_element
{
enum class Flags
{
  /** Fill nothing */
  none           =    0,

  /** Quadrature points on the element */
  value          =    1L << 1,

  /** Quadrature weigths on the element */
  gradient       =    1L << 2
};
}

template <int,int> class PhysicalDomain;
template <int, int, int, int, class> class FunctionElementBase;
template <int, int, int, int> class FunctionElement;
template <int, int, int, int> class ConstFunctionElement;

/**
 * Function Class
 *
 * @ingroup serializable
 */
template<int dim_, int codim_ = 0, int range_ = 1, int rank_ = 1>
class Function :
  public std::enable_shared_from_this<Function<dim_,codim_,range_,rank_> >
{
private:
  using self_t = Function<dim_, codim_, range_, rank_>;

public:
  static const int space_dim = dim_ + codim_;
  static const int dim       = dim_;
  static const int codim     = codim_;
  static const int range     = range_;
  static const int rank      = rank_;

  using DomainType = PhysicalDomain<dim_, codim_>;

  using ElementAccessor = FunctionElement<dim_, codim_, range_, rank_>;
  using ElementIterator = GridIterator<ElementAccessor>;
  using ConstElementAccessor = ConstFunctionElement<dim_, codim_, range_, rank_>;
  using ConstElementIterator = GridIterator<ConstElementAccessor>;

  using List = typename DomainType::List;
  using ListIt = typename DomainType::ListIt;
  using Flags = function_element::Flags;


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

protected:
  using FlagsArray = SafeSTLArray<Flags, dim+1>;
  /** @name Constructors and destructor. */
  ///@{
protected:
  /**
   * Default constructor. It does nothing but it is needed for the
   * <a href="http://www.boost.org/doc/libs/release/libs/serialization/">boost::serialization</a>
   * mechanism.
   */
  Function() = default;

  /** Constructor */
  Function(std::shared_ptr<const DomainType> phys_dom);

  /**
   * Copy constructor.
   */
  Function(const self_t &func);

public:
  /** Destructor */
  virtual ~Function() = default;
  ///@}

  std::shared_ptr<const DomainType> get_physical_domain() const
  {
    return phys_domain_;
  }

  //virtual std::shared_ptr<base_t> clone() const = 0;


  virtual void set_flags(const topology_variant &sdim,
                         const Flags &flag);

  virtual void init_cache(ElementAccessor &elem,
                          eval_pts_variant &quad) const;

  void init_cache(ElementIterator &elem,
                  eval_pts_variant &quad) const
  {
    this->init_cache(*elem, quad);
  }


//    void init_element_cache(ElementAccessor &elem) const;
//
//    void init_element_cache(ElementIterator &elem) const;

  virtual void fill_cache(const topology_variant &sdim,
                          ElementAccessor &elem,
                          const int s_id) const;

  void fill_cache(const topology_variant &sdim,
                  ElementIterator &elem,
                  const int s_id) const
  {
    this->fill_cache(sdim, *elem, s_id);
  }

//    void fill_element_cache(ElementAccessor &elem) const;
//
//    void fill_element_cache(ElementIterator &elem) const;


public:
  std::shared_ptr<ConstElementAccessor>
  create_element(const ListIt &index, const PropId &prop) const;



public:
  /** @name Functions involving the element iterator */
  ///@{
  /**
   * Returns a element iterator to the first element of the patch
   * with the property @p element_property.
   */
  ElementIterator begin(const PropId &element_property = ElementProperties::active);

  /**
   * Returns a element iterator to one-pass the end of patch
   * with the property @p element_property.
   */
  ElementIterator end(const PropId &element_property = ElementProperties::active);
  ///@}

  virtual void print_info(LogStream &out) const
  {
    Assert(false, ExcNotImplemented());
  }

//    std::shared_ptr<typename ElementAccessor::CacheType>
//    &get_cache(ElementAccessor &elem);

private:
  std::shared_ptr<const DomainType> phys_domain_;

protected:
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
                        ElementAccessor &elem)
      :
      flags_(flags),
      elem_(elem)
    {}

    template<int sdim>
    void operator()(std::shared_ptr<Quadrature<sdim>> &quad)
    {
      auto &cache = elem_.all_sub_elems_cache_;
      if (cache == nullptr)
      {
        using Cache = typename ElementAccessor::CacheType;
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
    ElementAccessor &elem_;
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
                                         boost::serialization::base_object<PhysicalDomain>(*this));

      ar &boost::serialization::make_nvp("flags_",flags_);
  }
  //*/
  ///@}
#endif // SERIALIZATION
};

IGA_NAMESPACE_CLOSE

#endif
