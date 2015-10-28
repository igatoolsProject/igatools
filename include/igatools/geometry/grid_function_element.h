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
template<int dim_, int space_dim_, class ContainerType_>
class GridFunctionElementBase
{
private:
  using self_t  = GridFunctionElementBase<dim_, space_dim_, ContainerType_>;

public:
  using ContainerType = ContainerType_;
  using GridElem = typename ContainerType_::GridType::ElementAccessor;
  using ListIt = typename ContainerType_::ListIt;

  using IndexType = typename Grid<dim_>::IndexType;

  using Value =  typename ContainerType_::Value;
  template <int order>
  using Derivative = typename ContainerType_::template Derivative<order>;


// using Gradient =  typename ContainerType_::Gradient;

  using Flags = grid_function_element::Flags;
  using CacheFlags = grid_function_element::CacheFlags;


  /** @name Constructors */
  ///@{
protected:
  /**
   * Default constructor. It does nothing but it is needed for the
   * <a href="http://www.boost.org/doc/libs/release/libs/serialization/">boost::serialization</a>
   * mechanism.
   */
  GridFunctionElementBase() = default;

public:
  /**
   * Construct an accessor pointing to the element with
   * flat index @p elem_index of the Function @p func.
   */
  GridFunctionElementBase(const std::shared_ptr<ContainerType_> grid_function,
                          const ListIt &index,
                          const PropId &prop = ElementProperties::active);

  /**
   * Copy constructor. Not allowed to be used.
   */
  GridFunctionElementBase(const self_t &elem) = delete;

  /**
   * Move constructor.
   */
  GridFunctionElementBase(self_t &&elem) = default;

  /**
   * Destructor.
   */
  ~GridFunctionElementBase() = default;
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
  bool operator==(const self_t &elem) const;

  /**
   * True if the elements have different index.
   *  @note In debug mode, it is also check they both refer to
   *  the same cartesian grid. No check is done on the cache.
   */
  bool operator!=(const self_t &elem) const;

  /**
   * True if the flat-index of the element on the left is smaller than
   * the flat-index of the element on the right.
   *  @note In debug mode, it is also check they both refer to
   *  the same cartesian grid. No check is done on the cache.
   */
  bool operator<(const self_t &elem) const;

  /**
   * True if the flat-index of the element on the left is bigger than
   * the flat-index of the element on the right.
   *  @note In debug mode, it is also check they both refer to
   *  the same cartesian grid. No check is done on the cache.
   */
  bool operator>(const self_t &elem) const;
  ///@}


public:

  void operator++()
  {
    ++(*grid_elem_);
  }


  void move_to(const IndexType &elem_id)
  {
    grid_elem_->move_to(elem_id);
  }

  const GridElem &get_grid_element() const
  {
    return *grid_elem_;
  }

  GridElem &get_grid_element()
  {
    return *grid_elem_;
  }


  void print_info(LogStream &out) const
  {
    using std::to_string;
    out.begin_item("GridElement<" + to_string(dim_) + "," + to_string(space_dim_) +">");
    grid_elem_->print_info(out);
    out.end_item();
  }

  void print_cache_info(LogStream &out) const
  {
    out.begin_item("GridElement's cache");
    grid_elem_->print_cache_info(out);
    out.end_item();

    local_cache_.print_info(out);
  }

public:

  template <class ValueType, int sdim>
  const auto &get_values_from_cache(const int s_id = 0) const
  {
    const auto &cache = local_cache_.template
                        get_sub_elem_cache<sdim>(s_id);
    return cache.template get_data<ValueType>();
  }


public:
  template <int order>
  using _D = grid_function_element::_D<order>;

private:
  using CType = boost::fusion::map<
                boost::fusion::pair<_D<0>, DataWithFlagStatus<ValueVector<Value>>>,
                boost::fusion::pair<_D<1>,DataWithFlagStatus<ValueVector<Derivative<1>>>>,
                boost::fusion::pair<_D<2>,DataWithFlagStatus<ValueVector<Derivative<2>>>>
                >;

public:

  using Cache = PointValuesCache<dim_,CType>;

  using CacheType = AllSubElementsCache<Cache>;

protected:
  std::shared_ptr<ContainerType_> grid_function_;

private:
  std::unique_ptr<GridElem> grid_elem_;

  CacheType local_cache_;

  template <class Accessor> friend class GridIteratorBase;
  friend class GridFunctionHandler<dim_, space_dim_>;
};



template <int dim, int space_dim>
class ConstGridFunctionElement
  : public GridFunctionElementBase<dim, space_dim,
    const GridFunction<dim,space_dim>>
{
  using GridFunctionElementBase<dim, space_dim,
        const GridFunction<dim,space_dim>>::GridFunctionElementBase;

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
  decltype(auto) evaluate_at_points(const std::shared_ptr<const Quadrature<dim>> &points)
  {
    auto grid_func_elem_handler = this->grid_function_->create_cache_handler();
    grid_func_elem_handler->template set_flags<dim>(ValueType::flag);
    grid_func_elem_handler->init_cache(*this,points);
    grid_func_elem_handler->template fill_cache<dim>(*this,0);

    return this->template get_values_from_cache<ValueType,dim>(0);
  }
  ///@}

};



template <int dim, int space_dim>
class GridFunctionElement
  : public GridFunctionElementBase<dim, space_dim,
    GridFunction<dim,space_dim>>
{
  using GridFunctionElementBase<dim, space_dim,
        GridFunction<dim,space_dim>>::GridFunctionElementBase;
};

IGA_NAMESPACE_CLOSE

#endif


