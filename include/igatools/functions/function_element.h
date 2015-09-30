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

#ifndef __FUNCTION_ELEMENT_H_
#define __FUNCTION_ELEMENT_H_


#include <igatools/functions/function.h>
#include <igatools/functions/function_handler.h>
#include <igatools/geometry/domain_element.h>

//#include <igatools/base/value_types.h>
//#include <igatools/basis_functions/values_cache.h>

IGA_NAMESPACE_OPEN

/**
 *
 * @ingroup elements
 */
template<int dim, int codim, int range, int rank, class ContainerType_>
class FunctionElementBase
{
private:
  using self_t = FunctionElementBase<dim, codim, range, rank, ContainerType_>;

public:

  using ContainerType = ContainerType_;
  using DomainElem = typename ContainerType_::DomainType::ConstElementAccessor;

  using ListIt = typename ContainerType_::ListIt;

  using Value = typename ContainerType_::Value;
//  using Gradient = typename ContainerType_::Gradient;
//  using Hessian  = typename ContainerType_::Hessian;
//  using Div      = typename ContainerType_::Div;
  template <int order>
  using Derivative = typename ContainerType_::template Derivative<order>;

  using Flags = function_element::Flags;
  using CacheFlags = function_element::CacheFlags;

protected:
#if 0
  /** @name Constructors */
  ///@{
  /**
   * Default constructor. It does nothing but it is needed for the
   * <a href="http://www.boost.org/doc/libs/release/libs/serialization/">boost::serialization</a>
   * mechanism.
   */
  FunctionElementBase() = default;
#endif

  FunctionElementBase() = delete;

public:
  /**
   * Construct an accessor pointing to the element with
   * flat index @p elem_index of the Grid @p grid.
   */
  FunctionElementBase(const std::shared_ptr<ContainerType_> func,
                      const ListIt &index,
                      const PropId &prop = ElementProperties::active);

  /**
   * Copy constructor. Not allowed to be used.
   */
  FunctionElementBase(const self_t &elem) = delete;

  /**
   * Move constructor.
   */
  FunctionElementBase(self_t &&elem) = default;

  /**
   * Destructor.
   */
  ~FunctionElementBase() = default;
  ///@}



  /** @name Assignment operators */
  ///@{
  /**
   * Copy assignment operator. Not allowed to be used.
   */
  self_t &operator=(const self_t &element) = delete;

  /**
   * Move assignment operator.
   */
  self_t &operator=(self_t &&elem) = default;
  ///@}

  /**
     * @name Comparison operators.
     *
     * @brief The comparison operators compares the <em>position</em> of the element in the grid.
     *
     * @warning To be comparable, two Function objects must be defined using the same Function
     * (and therefore on the same grid),
     * otherwise an assertion will be raised (in Debug mode).
     */
  ///@{
  /** Returns TRUE if the two elements have the same index on the grid. */
  bool operator==(const self_t &a) const;


  /** Returns TRUE if the two elements have different indices on the grid. */
  bool operator!=(const self_t &a) const;

  /**
   * Returns TRUE if the the index of the element on the left of the operator <tt> < </tt>
   * is smaller than the the index of the element on the right.
   * */
  bool operator<(const self_t &a) const;

  /**
   * Returns TRUE if the the index of the element on the left of the operator <tt> < </tt>
   * is bigger than the the index of the element on the right.
   * */
  bool operator>(const self_t &a) const;
  ///@}

  /*
  ListIt &operator++()
  {
    return (++(*domain_elem_));
  }
  //*/

  void operator++()
  {
    ++(*domain_elem_);
  }

  const DomainElem &get_domain_element() const;

  DomainElem &get_domain_element();

  void print_info(LogStream &out) const;

  void print_cache_info(LogStream &out) const;

  template<class ValueType, int sdim>
  auto
  get_values(const int s_id) const
  {
    Assert(local_cache_ != nullptr,ExcNullPtr());
    const auto &cache =
      local_cache_->template get_sub_elem_cache<sdim>(s_id);
    return cache.template get_data<ValueType>();
  }


#if 0
  /**
   * @name Methods for the for the evaluations of Functions's derivatives
   *  without the use of the cache.
   */
  ///@{
  /**
   * Returns a ValueTable with the values specified by the template parameter
   * <tt>ValueType</tt>
   * at each point (in the unit domain) specified by the input argument <tt>points</tt>.
   * @note This function does not use the cache and therefore can be called any time without
   * needing to pre-call init_cache()/fill_cache().
   * @warning The evaluation <tt>points</tt> must belong to the unit hypercube
   * \f$ [0,1]^{\text{dim}} \f$ otherwise, in Debug mode, an assertion will be raised.
   */
  template <class ValueType>
  decltype(auto) evaluate_at_points(const Quadrature<dim> &points)
  {
    func_->reset(ValueType::flag,points);
    const auto topology = Topology<dim>();
    func_->init_cache(*this,topology);
    func_->fill_cache(*this,topology,0);

    return this->template get_values<ValueType,dim>(0);
  }
  ///@}

#endif





  using _Value = function_element::_Value;
  using _Gradient = function_element::_Gradient;
  using _D2 = function_element::_D2;

private:


  using CType = boost::fusion::map<
                boost::fusion::pair<     _Value,DataWithFlagStatus<ValueVector<Value>>>,
                boost::fusion::pair<  _Gradient,DataWithFlagStatus<ValueVector<Derivative<1>>>>,
                boost::fusion::pair<  _D2,DataWithFlagStatus<ValueVector<Derivative<2>>>>
                >;



  using Cache = PointValuesCache<dim,CType>;

public:
  using CacheType = AllSubElementsCache<Cache>;

private:
//  std::shared_ptr<CacheType>
  CacheType &
  get_cache()
  {
//    Assert(local_cache_ != nullptr,ExcNullPtr());
    return local_cache_;
  }


private:
  std::shared_ptr<ContainerType_> func_;

  std::unique_ptr<DomainElem> domain_elem_;

  CacheType local_cache_;

  template <class Accessor> friend class GridIteratorBase;
  friend class FunctionHandler<dim, codim, range, rank>;






//
//    /** Returns the index of the element. */
//    IndexType get_index() const;



private:

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
  serialize(Archive &ar, const unsigned int version);
  ///@}
#endif // SERIALIZATION
#endif
};



template <int dim, int codim, int range, int rank>
class ConstFunctionElement
  : public FunctionElementBase<dim, codim, range, rank,
    const Function<dim,codim,range,rank> >
{
  using FunctionElementBase<dim, codim, range, rank,
        const Function<dim,codim,range,rank>>::FunctionElementBase;
};


template <int dim, int codim, int range, int rank>
class FunctionElement
  : public FunctionElementBase<dim, codim, range, rank, Function<dim,codim,range,rank> >
{

  using FunctionElementBase<dim, codim, range, rank,
        Function<dim,codim,range,rank>>::FunctionElementBase;
};

IGA_NAMESPACE_CLOSE

#endif
