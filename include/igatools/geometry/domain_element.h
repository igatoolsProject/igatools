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

#ifndef __DOMAIN_ELEMENT_H_
#define __DOMAIN_ELEMENT_H_

#include <igatools/geometry/domain.h>
#include <igatools/geometry/domain_handler.h>
#include <igatools/geometry/grid_function_element.h>

IGA_NAMESPACE_OPEN

/**
 *
 * @ingroup elements
 */
template<int dim_, int codim_, class ContainerType_>
class DomainElementBase
{
private:
  using self_t  = DomainElementBase<dim_, codim_, ContainerType_>;

public:
  using ContainerType = ContainerType_;
  using GridFuncElem = typename ContainerType_::GridFuncType::ConstElementAccessor;
  using ListIt = typename ContainerType_::ListIt;

  using Point =  typename ContainerType_::Point;
  using Gradient =  typename ContainerType_::Gradient;

  using Flags = domain_element::Flags;
  using CacheFlags = domain_element::CacheFlags;


  /** @name Constructors */
  ///@{
protected:
  /**
   * Default constructor. It does nothing but it is needed for the
   * <a href="http://www.boost.org/doc/libs/release/libs/serialization/">boost::serialization</a>
   * mechanism.
   */
  DomainElementBase() = default;

public:
  /**
   * Construct an accessor pointing to the element with
   * flat index @p elem_index of the Function @p func.
   */
  DomainElementBase(const std::shared_ptr<ContainerType_> domain,
                    const ListIt &index,
                    const PropId &prop = ElementProperties::active);

  /**
   * Copy constructor. Not allowed to be used.
   */
  DomainElementBase(const self_t &elem) = delete;

  /**
   * Move constructor.
   */
  DomainElementBase(self_t &&elem) = default;

  /**
   * Destructor.
   */
  ~DomainElementBase() = default;
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

#if 0

  template <int order>
  using Derivative = typename FuncType::template Derivative<order>;
#endif

public:
  ListIt &operator++()
  {
    return (++(*grid_func_elem_));
  }


  const GridFuncElem &get_grid_function_element() const
  {
    return *grid_func_elem_;
  }

  GridFuncElem &get_grid_function_element()
  {
    return *grid_func_elem_;
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
  auto const get_points(const int s_id) const
  {
    return grid_func_elem_->template
           get_values<grid_function_element::_D<0>, sdim>(s_id);
  }

  template<int sdim>
  ValueVector<Real> get_w_measures(const int s_id) const;

  ValueVector<SafeSTLArray<Point, codim_> >
  get_exterior_normals() const;

#if 0
  using MetricTensor =
    Tensor<dim, 1, tensor::covariant, Tensor<dim, 1, tensor::contravariant, Tdouble> >;

  ValueVector<MetricTensor> compute_inv_first_fundamental_form() const;

  ValueVector<MetricTensor> compute_second_fundamental_form() const;

  ValueVector< Derivative<1> > get_D_external_normals() const;

  const ValueVector<SafeSTLVector<Real> > &get_principal_curvatures() const;


  template<int sub_dim>
  const ValueVector<Points<space_dim> > &
  get_boundary_normals(const int s_id) const
  {
#if 0
    Assert(dim==sub_dim+1, ExcNotImplemented());
    ValueVector<Points<space_dim>> res;
    const auto &DF_inv = get_values_from_cache<_InvGradient, sub_dim>(s_id);
    const auto n_hat  = this->get_grid()->template get_boundary_normals<sub_dim>(s_id)[0];

    const auto n_points = DF_inv.get_num_points();
    res.resize(n_points);
    for (int pt = 0; pt < n_points; ++pt)
    {
      const auto DF_inv_t = co_tensor(transpose(DF_inv[pt]));
      res[pt] = action(DF_inv_t, n_hat);
      res[pt] /= res[pt].norm();
    }
    return res;
#endif
    return get_values_from_cache<_BoundaryNormal,sub_dim>(s_id);
  }
#endif

public:
  template <class ValueType, int sdim>
  auto &get_values_from_cache(const int s_id = 0) const
  {
    const auto &cache = local_cache_.template
                        get_sub_elem_cache<sdim>(s_id);
    return cache.template get_data<ValueType>();
  }


  using _Measure = domain_element::_Measure;

  using _InvJacobian = domain_element::_InvJacobian;

  using _BoundaryNormal = domain_element::_BoundaryNormal;

  using _ExtNormal = domain_element::_ExtNormal;

private:
  template<int order>
  using InvDerivative = Derivatives<dim_+codim_,dim_,1,order>;


  using CType = boost::fusion::map<
                boost::fusion::pair<_Measure       ,DataWithFlagStatus<ValueVector<Real>> >,
                boost::fusion::pair<_InvJacobian   ,DataWithFlagStatus<ValueVector<InvDerivative<1>>>>,
                boost::fusion::pair<_BoundaryNormal,DataWithFlagStatus<ValueVector<Points<dim_+codim_>>>>,
                boost::fusion::pair<_ExtNormal     ,DataWithFlagStatus<ValueVector<SafeSTLArray<Point,codim_>>>>
                >;
//                ,
//                  boost::fusion::pair<    _InvHessian,DataWithFlagStatus<ValueVector<InvDerivative<2>>>>,
//                  boost::fusion::pair<     _Curvature,DataWithFlagStatus<ValueVector<SafeSTLVector<Real>>>>
//                  >;

  using Cache = PointValuesCache<dim_,CType>;



public:
  using CacheType = AllSubElementsCache<Cache>;

private:
  std::shared_ptr<ContainerType_> domain_;

  std::unique_ptr<GridFuncElem> grid_func_elem_;

  CacheType local_cache_;

  template <class Accessor> friend class GridIteratorBase;
  friend class DomainHandler<dim_, codim_>;
};



template <int dim, int codim>
class ConstDomainElement
  : public DomainElementBase<dim, codim,
    const Domain<dim,codim>>
{
  using DomainElementBase<dim, codim,
        const Domain<dim,codim>>::DomainElementBase;
};



template <int dim, int codim>
class DomainElement
  : public DomainElementBase<dim, codim,
    Domain<dim,codim>>
{
  using DomainElementBase<dim, codim,
        Domain<dim,codim>>::DomainElementBase;
};

IGA_NAMESPACE_CLOSE

#endif // PHYSICAL_DOMAIN_ELEMENT_H_


