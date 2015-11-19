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
template<int dim, int codim, int range, int rank>
class FunctionElement
{
private:
  using self_t = FunctionElement<dim, codim, range, rank>;

public:

  using ContainerType = const Function<dim,codim,range,rank>;
  using DomainElem = typename ContainerType::DomainType::ElementAccessor;

  using ListIt = typename ContainerType::ListIt;

  using IndexType = typename Grid<dim>::IndexType;

  using Value = typename ContainerType::Value;
//  using Gradient = typename ContainerType::Gradient;
//  using Hessian  = typename ContainerType::Hessian;
//  using Div      = typename ContainerType::Div;
  template <int order>
  using Derivative = typename ContainerType::template Derivative<order>;

  using Flags = function_element::Flags;

protected:

  FunctionElement() = delete;

public:
  /**
   * Construct an accessor pointing to the element with
   * flat index @p elem_index of the Grid @p grid.
   */
  FunctionElement(const std::shared_ptr<ContainerType> &func,
                  const ListIt &index,
                  const PropId &prop = ElementProperties::active);

  /**
   * Copy constructor. Not allowed to be used.
   */
  FunctionElement(const self_t &elem) = delete;

  /**
   * Move constructor.
   */
  FunctionElement(self_t &&elem) = default;

  /**
   * Destructor.
   */
  ~FunctionElement() = default;
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


  void operator++();

  void move_to(const IndexType &elem_id);

  const IndexType &get_index() const;

  const DomainElem &get_domain_element() const;

  DomainElem &get_domain_element();

  void print_info(LogStream &out) const;

  void print_cache_info(LogStream &out) const;

  template<class ValueType, int sdim>
  auto
  get_values_from_cache(const int s_id) const
  {
    const auto &cache =
      local_cache_.template get_sub_elem_cache<sdim>(s_id);
    return cache.template get_data<ValueType>();
  }


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
  decltype(auto) evaluate_at_points(const std::shared_ptr<const Quadrature<dim>> &points)
  {
    auto func_elem_handler = this->func_->create_cache_handler();
    func_elem_handler->template set_flags<dim>(ValueType::flag);
    func_elem_handler->init_cache(*this,points);
    func_elem_handler->template fill_cache<dim>(*this,0);

    return this->template get_values_from_cache<ValueType,dim>(0);
  }
  ///@}





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


protected:
  std::shared_ptr<ContainerType> func_;

private:
  std::unique_ptr<DomainElem> domain_elem_;

  CacheType local_cache_;

  template <class Accessor> friend class GridIteratorBase;
  friend class FunctionHandler<dim, codim, range, rank>;
};


IGA_NAMESPACE_CLOSE

#endif
