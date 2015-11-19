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

#ifndef __GRID_FUNCTION_ELEMENT_H_
#define __GRID_FUNCTION_ELEMENT_H_

#include <igatools/geometry/grid_function.h>
#include <igatools/geometry/grid_function_handler.h>
#include <igatools/geometry/grid_element.h>

IGA_NAMESPACE_OPEN

/**
 *
 * @ingroup elements
 */
template<int dim_, int space_dim_>
class GridFunctionElement
{
private:
  using self_t  = GridFunctionElement<dim_, space_dim_>;

public:
  using ContainerType = const GridFunction<dim_,space_dim_>;
  using GridElem = typename ContainerType::GridType::ElementAccessor;
  using ListIt = typename ContainerType::ListIt;

  using IndexType = typename Grid<dim_>::IndexType;

  using Value =  typename ContainerType::Value;
  template <int order>
  using Derivative = typename ContainerType::template Derivative<order>;


// using Gradient =  typename ContainerType_::Gradient;

  using Flags = grid_function_element::Flags;


  /** @name Constructors */
  ///@{
protected:
  /**
   * Default constructor. Not allowed to be used.
   */
  GridFunctionElement() = delete;

public:
  /**
   * Construct an accessor pointing to the element with
   * flat index @p elem_index of the Function @p func.
   */
  GridFunctionElement(const std::shared_ptr<ContainerType> &grid_function,
                      const ListIt &index,
                      const PropId &prop = ElementProperties::active);

  /**
   * Copy constructor. Not allowed to be used.
   */
  GridFunctionElement(const self_t &elem) = delete;

  /**
   * Move constructor.
   */
  GridFunctionElement(self_t &&elem) = default;

  /**
   * Destructor.
   */
  virtual ~GridFunctionElement() = default;
  ///@}


  /**
   * @name Comparison operators
   * @note In order to be meaningful, the comparison must be performed on
   * elements defined on
   * the <b>same</b> GridFunction
   * (in the sense that the pointer to the GridFunction held by the elements must
   * point to the same GridFunction object).
   */
  ///@{
  /**
   * True if the elements have the same index.
   *  @note In debug mode, it is also check they both refer to
   *  the same cartesian grid. No check is done on the cache.
   */
  virtual bool operator==(const self_t &elem) const;

  /**
   * True if the elements have different index.
   *  @note In debug mode, it is also check they both refer to
   *  the same cartesian grid. No check is done on the cache.
   */
  virtual bool operator!=(const self_t &elem) const;

  /**
   * True if the flat-index of the element on the left is smaller than
   * the flat-index of the element on the right.
   *  @note In debug mode, it is also check they both refer to
   *  the same cartesian grid. No check is done on the cache.
   */
  virtual bool operator<(const self_t &elem) const;

  /**
   * True if the flat-index of the element on the left is bigger than
   * the flat-index of the element on the right.
   *  @note In debug mode, it is also check they both refer to
   *  the same cartesian grid. No check is done on the cache.
   */
  virtual bool operator>(const self_t &elem) const;
  ///@}



  virtual void operator++();


  virtual void move_to(const IndexType &elem_id);

  const GridElem &get_grid_element() const;

  GridElem &get_grid_element();


  const IndexType &get_index() const;

  virtual void print_info(LogStream &out) const;

  void print_cache_info(LogStream &out) const;


public:

  template <class ValueType, int sdim>
  const auto &get_values_from_cache(const int s_id = 0) const
  {
    const auto &cache = local_cache_.template
                        get_sub_elem_cache<sdim>(s_id);
    return cache.template get_data<ValueType>();
  }

  const ValueVector<Value> &
  get_element_values_D0() const;

  const ValueVector<Derivative<1> > &
  get_element_values_D1() const;

  const ValueVector<Derivative<2> > &
  get_element_values_D2() const;


public:
  template <int order>
  using _D = grid_function_element::_D<order>;

private:

  template <class ValueType>
  struct IsInCache
  {
    const static bool value =
      std::is_same<ValueType,_D<0>>::value ||
      std::is_same<ValueType,_D<1>>::value ||
      std::is_same<ValueType,_D<2>>::value;
  };

  using CType = boost::fusion::map<
                boost::fusion::pair<_D<0>, DataWithFlagStatus<ValueVector<Value>>>,
                boost::fusion::pair<_D<1>,DataWithFlagStatus<ValueVector<Derivative<1>>>>,
                boost::fusion::pair<_D<2>,DataWithFlagStatus<ValueVector<Derivative<2>>>>
                >;

public:

  using Cache = PointValuesCache<dim_,CType>;

  using CacheType = AllSubElementsCache<Cache>;

protected:
  std::shared_ptr<ContainerType> grid_function_;

  std::unique_ptr<GridElem> grid_elem_;

private:
  CacheType local_cache_;

  template <class Accessor> friend class GridIteratorBase;
  friend class GridFunctionHandler<dim_, space_dim_>;


public:
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
  decltype(auto) evaluate_at_points(const std::shared_ptr<const Quadrature<dim_>> &quad,
                                    EnableIf< IsInCache<ValueType>::value > * = nullptr)
  {
    auto elem_handler = this->grid_function_->create_cache_handler();
    elem_handler->template set_flags<dim_>(ValueType::flag);
    elem_handler->init_cache(*this,quad);
    elem_handler->template fill_cache<dim_>(*this,0);

    return this->template get_values_from_cache<ValueType,dim_>(0);
  }

  template <class ValueType>
  decltype(auto) evaluate_at_points(const std::shared_ptr<const Quadrature<dim_>> &quad,
                                    EnableIf< !(IsInCache<ValueType>::value) > * = nullptr)
  {
    return grid_elem_->template evaluate_at_points<ValueType>(quad);
  }
  ///@}


  /**
   * Returns the quadrature weights corresponding to the <tt>dim</tt>
   * dimensional element (i.e. the element itself).
   *
   * @note The returned weights are the quadrature unit weights multiplied by the
   * <tt>dim</tt>-dimensional element measure.
   */
  const ValueVector<Real> &get_element_weights() const;

  /**
   * Returns the quadrature points corresponding to the <tt>dim</tt>
   * dimensional element (i.e. the element itself).
   */
  const ValueVector<Points<dim_>> &get_element_points() const;

public:

  bool same_grid_function_of(const self_t &elem) const;

};






IGA_NAMESPACE_CLOSE

#endif


