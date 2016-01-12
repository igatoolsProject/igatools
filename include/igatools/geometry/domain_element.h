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
template<int dim_, int codim_>
class DomainElement
{
private:
  using self_t  = DomainElement<dim_, codim_>;

public:
  using ContainerType = const Domain<dim_,codim_>;
  using GridFunc = typename ContainerType::GridFuncType;
  using GridFuncElem = typename GridFunc::ElementAccessor;
  using ListIt = typename ContainerType::ListIt;

  using IndexType = typename Grid<dim_>::IndexType;

  using Point =  typename ContainerType::Point;
  using Jacobian = typename ContainerType::Gradient;
  using Hessian = typename ContainerType::template Derivative<2>;

  using Flags = domain_element::Flags;


  /** @name Constructors */
  ///@{
protected:
  /**
   * Default constructor. It does nothing but it is needed for the
   * <a href="http://www.boost.org/doc/libs/release/libs/serialization/">boost::serialization</a>
   * mechanism.
   */
  DomainElement() = default;

public:
  DomainElement(const std::shared_ptr<ContainerType> &domain,
                std::unique_ptr<GridFuncElem> &&grid_func_elem);

  /**
   * Copy constructor. Not allowed to be used.
   */
  DomainElement(const self_t &elem) = delete;

  /**
   * Move constructor.
   */
  DomainElement(self_t &&elem) = default;

  /**
   * Destructor.
   */
  ~DomainElement() = default;
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
  ///@}

#if 0

  template <int order>
  using Derivative = typename FuncType::template Derivative<order>;
#endif

public:

  void operator++();

  void move_to(const IndexType &elem_id);

  const IndexType &get_index() const;

  const GridFuncElem &get_grid_function_element() const;

  GridFuncElem &get_grid_function_element();


  void print_info(LogStream &out) const;

  void print_cache_info(LogStream &out) const;

public:


  template<int sdim>
  const ValueVector<Point> &get_points(const int s_id) const;

  template<int sdim>
  const ValueVector<Jacobian> &get_jacobians(const int s_id) const;

  template<int sdim>
  const ValueVector<Hessian> &get_hessians(const int s_id) const;

  template<int sdim>
  const ValueVector<Real> &get_measures(const int s_id) const;

  template<int sdim>
  const ValueVector<Real> &get_w_measures(const int s_id) const;


  const ValueVector<SafeSTLArray<Point, codim_> > &
  get_exterior_normals() const;

  const ValueVector<typename GridFunc::template Derivative<1> > &
  get_exterior_normals_D1() const;

  template <int sdim>
  const ValueVector<Points<dim_+codim_> > &
  get_boundary_normals(const int s_id, EnableIf<(sdim >= 0)> * = nullptr) const;


  const ValueVector<Point> &get_element_points() const;

  const ValueVector<Jacobian> &get_element_jacobians() const;

  const ValueVector<Hessian> &get_element_hessians() const;

  const ValueVector<Real> &get_element_measures() const;

  const ValueVector<Real> &get_element_w_measures() const;


  const ValueVector<SafeSTLVector<Real> > &get_principal_curvatures() const;


#if 0
  using MetricTensor =
    Tensor<dim, 1, tensor::covariant, Tensor<dim, 1, tensor::contravariant, Tdouble> >;

  ValueVector<MetricTensor> compute_inv_first_fundamental_form() const;

  ValueVector<MetricTensor> compute_second_fundamental_form() const;

  ValueVector< Derivative<1> > get_D_external_normals() const;



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

  using _W_Measure = domain_element::_W_Measure;

  using _InvJacobian = domain_element::_InvJacobian;

  using _InvHessian = domain_element::_InvHessian;

  using _BoundaryNormal = domain_element::_BoundaryNormal;

  using _ExtNormal = domain_element::_ExtNormal;
  using _ExtNormalD1 = domain_element::_ExtNormalD1;

  using _Curvature = domain_element::_Curvature;

  using _FirstFundamentalForm = domain_element::_FirstFundamentalForm;

  using _SecondFundamentalForm = domain_element::_SecondFundamentalForm;

private:
  template<int order>
  using Derivative = typename GridFunc::template Derivative<order>;

  template<int order>
  using InvDerivative = Derivatives<dim_+codim_,dim_,1,order>;

  using MetricTensor =
    Tensor<dim_,1,tensor::covariant,Tensor<dim_,1,tensor::contravariant,Tdouble> >;

  template <class ValueType>
  struct IsInCache
  {
    const static bool value =
      std::is_same<ValueType,_Measure>::value ||
      std::is_same<ValueType,_W_Measure>::value ||
      std::is_same<ValueType,_InvJacobian>::value ||
      std::is_same<ValueType,_InvHessian>::value ||
      std::is_same<ValueType,_BoundaryNormal>::value ||
      std::is_same<ValueType,_ExtNormal>::value ||
      std::is_same<ValueType,_Curvature>::value ||
      std::is_same<ValueType,_FirstFundamentalForm>::value ||
      std::is_same<ValueType,_SecondFundamentalForm>::value ||
      std::is_same<ValueType,_ExtNormalD1>::value;
  };

  using CType = boost::fusion::map<
                boost::fusion::pair<_Measure       ,DataWithFlagStatus<ValueVector<Real>> >,
                boost::fusion::pair<_W_Measure     ,DataWithFlagStatus<ValueVector<Real>> >,
                boost::fusion::pair<_InvJacobian   ,DataWithFlagStatus<ValueVector<InvDerivative<1>>>>,
                boost::fusion::pair<_InvHessian    ,DataWithFlagStatus<ValueVector<InvDerivative<2>>>>,
                boost::fusion::pair<_BoundaryNormal,DataWithFlagStatus<ValueVector<Points<dim_+codim_>>>>,
                boost::fusion::pair<_ExtNormal     ,DataWithFlagStatus<ValueVector<SafeSTLArray<Point,codim_>>>>,
                boost::fusion::pair<_Curvature     ,DataWithFlagStatus<ValueVector<SafeSTLVector<Real>>>>,
                boost::fusion::pair<_FirstFundamentalForm,DataWithFlagStatus<ValueVector<MetricTensor>>>,
                boost::fusion::pair<_SecondFundamentalForm,DataWithFlagStatus<ValueVector<MetricTensor>>>,
                boost::fusion::pair<_ExtNormalD1   ,DataWithFlagStatus<ValueVector<Derivative<1>>>>
                >;

  using Cache = PointValuesCache<dim_,CType>;



public:
  using CacheType = AllSubElementsCache<Cache>;

private:
  std::shared_ptr<ContainerType> domain_;

  std::unique_ptr<GridFuncElem> grid_func_elem_;

  CacheType local_cache_;

  template <class Accessor> friend class GridIteratorBase;
  friend class DomainHandler<dim_, codim_>;


public:
  template <class ValueType>
  decltype(auto) evaluate_at_points(const std::shared_ptr<const Quadrature<dim_>> &quad,
                                    EnableIf< IsInCache<ValueType>::value > * = nullptr)
  {
    auto elem_handler = this->domain_->create_cache_handler();
    elem_handler->set_element_flags(ValueType::flag);
    elem_handler->init_cache(*this,quad);
    elem_handler->fill_element_cache(*this);

    return this->template get_values_from_cache<ValueType,dim_>(0);
  }

  template <class ValueType>
  decltype(auto) evaluate_at_points(const std::shared_ptr<const Quadrature<dim_>> &quad,
                                    EnableIf< !(IsInCache<ValueType>::value) > * = nullptr)
  {
    return grid_func_elem_->template
           evaluate_at_points<typename ValueType::GridFuncElemType>(quad);
  }


  /**
   * Return TRUE if the element index is referring to a valid element in the Grid.
   *
   * @note An element with non valido position can happens when we use the ++ operator
   * on an element that is the last in the list.
   */
  bool has_valid_position() const
  {
    return grid_func_elem_->has_valid_position();
  }


};



IGA_NAMESPACE_CLOSE

#endif // PHYSICAL_DOMAIN_ELEMENT_H_


