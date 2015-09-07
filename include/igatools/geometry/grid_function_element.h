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

  using Value =  typename ContainerType_::Value;
  using Gradient =  typename ContainerType_::Gradient;

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
   * the <b>same grid</b>
   * (in the sense that the pointer to the grid held by the element must
   * point to the same grid object).
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
  ListIt &operator++()
  {
    return (++(*grid_elem_));
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
    AssertThrow(false,ExcNotImplemented());
  }

  void print_cache_info(LogStream &out) const
  {
    AssertThrow(false,ExcNotImplemented());
  }

public:
  template<int sdim>
  ValueVector<Real> const &get_measures(const int s_id) const
  {
    return get_values_from_cache<_Measure,sdim>(s_id);
  }

  template<int sdim>
  auto const &get_points(const int s_id) const
  {
    return get_values_from_cache<_Point,sdim>(s_id);
  }

  template<int sdim>
  ValueVector<Real> get_w_measures(const int s_id) const;

#if 0
  ValueVector<SafeSTLArray<Value, space_dim_> >
  get_exterior_normals() const;
#endif

private:
  template <class ValueType, int sdim>
  auto &get_values_from_cache(const int s_id = 0) const
  {
    Assert(local_cache_ != nullptr,ExcNullPtr());
    const auto &cache = local_cache_->template
                        get_sub_elem_cache<sdim>(s_id);
    return cache.template get_data<ValueType>();
  }

public:
  using _Point     = grid_function_element::_Point;
  using _Measure   = grid_function_element::_Measure;
  using _Gradient  = grid_function_element::_Gradient;

private:
  using CType = boost::fusion::map<
                boost::fusion::pair< _Point,    DataWithFlagStatus<ValueVector<Value>> >,
                boost::fusion::pair< _Measure,  DataWithFlagStatus<ValueVector<Real>> >,
                boost::fusion::pair< _Gradient, DataWithFlagStatus<ValueVector<Gradient>> >
                >;
//                ,
//                  boost::fusion::pair<   _InvGradient,DataWithFlagStatus<ValueVector<InvDerivative<1>>>>,
//                  boost::fusion::pair<    _InvHessian,DataWithFlagStatus<ValueVector<InvDerivative<2>>>>,
//                  boost::fusion::pair<_BoundaryNormal,DataWithFlagStatus<ValueVector<Points<space_dim>>>>,
//                  boost::fusion::pair<   _OuterNormal,DataWithFlagStatus<ValueVector<Points<space_dim>>>>,
//                  boost::fusion::pair<     _Curvature,DataWithFlagStatus<ValueVector<SafeSTLVector<Real>>>>
//                  >;

  using Cache = PointValuesCache<dim_,CType>;

public:
  using CacheType = AllSubElementsCache<Cache>;

private:
  std::shared_ptr<ContainerType_> grid_function_;

  std::unique_ptr<GridElem> grid_elem_;

  std::shared_ptr<CacheType> local_cache_;

  template <class Accessor> friend class GridIteratorBase;
  friend class GridFunctionHandler<dim_, space_dim_>;
};



template <int dim, int codim>
class ConstGridFunctionElement
  : public GridFunctionElementBase<dim, codim,
    const GridFunction<dim,codim>>
{
  using GridFunctionElementBase<dim, codim,
        const GridFunction<dim,codim>>::GridFunctionElementBase;
};



template <int dim, int codim>
class GridFunctionElement
  : public GridFunctionElementBase<dim, codim,
    GridFunction<dim,codim>>
{
  using GridFunctionElementBase<dim, codim,
        GridFunction<dim,codim>>::GridFunctionElementBase;
};

IGA_NAMESPACE_CLOSE

#endif


