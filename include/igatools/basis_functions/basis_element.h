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


#ifndef __BASIS_ELEMENT_H_
#define __BASIS_ELEMENT_H_

#include <igatools/base/config.h>
#include <igatools/basis_functions/values_cache.h>

#include <igatools/base/quadrature.h>

#include <igatools/utils/static_multi_array.h>
#include <igatools/utils/cartesian_product_indexer.h>

#include <igatools/basis_functions/spline_space.h>

#include <igatools/basis_functions/basis.h>
#include <igatools/geometry/grid_element.h>

#include <igatools/linear_algebra/dense_vector.h>


IGA_NAMESPACE_OPEN



/**
 *
 * @ingroup elements
 */
template<int dim_,int codim_,int range_,int rank_>
class BasisElement
{
protected:

  using Bs = const Basis<dim_,codim_,range_,rank_>;

private:
  using self_t = BasisElement<dim_,codim_,range_,rank_>;

public:

  using GridType = Grid<dim_>;
  using IndexType = typename GridType::IndexType;
  using List = typename GridType::List;
  using ListIt = typename GridType::ListIt;

  using GridElem = GridElement<dim_>;



  using Point = typename Bs::Point;
  using Value = typename Bs::Value;
  template <int order>
  using Derivative = typename Bs::template Derivative<order>;
  using Div = typename Bs::Div;

  static const int dim = dim_;
  static const int space_dim = dim_+codim_;




  using ContainerType = Bs;

  /**
   * For each component gives a product array of the dimension
   */
  template<class T>
  using ComponentContainer = typename SplineSpace<dim_,range_,rank_>::template ComponentContainer<T>;
  using TensorSizeTable = typename SplineSpace<dim_,range_,rank_>::TensorSizeTable;
  ///@}


  using Flags = basis_element::Flags;


private:
  /** @name Constructors */
  ///@{
  /**
   * \brief Default constructor. Not allowed to be used.
   */
  BasisElement() = delete;

  /**
   * \brief Copy constructor. Not allowed to be used.
   */
  BasisElement(const self_t &elem) = delete;

  /**
   * \brief Move constructor. Not allowed to be used.
   */
  BasisElement(self_t &&elem) = delete;
  ///@}


  /** @name Assignment operators */
  ///@{
  /**
   * \brief Copy assignment operator. Not allowed to be used.
   */
  self_t &operator=(const self_t &element) = delete;

  /**
   * \brief Move assignment operator.
   */
  self_t &operator=(self_t &&elem) = default;
  ///@}

public:

  BasisElement(const std::shared_ptr<Bs> &basis);



  /**
   * \brief Destructor.
   */
  virtual ~BasisElement() = default;




  /**
   * @name Comparison operators.
   *
   * @brief The comparison operators compares the <em>position</em> of the element in the grid.
   *
   * @warning To be comparable, two BasisElement objects must be defined on the same Basis
   * (and therefore on the same grid),
   * otherwise an assertion will be raised (in Debug mode).
   */
  ///@{
  /** \brief Returns TRUE if the two elements have the same index on the grid. */
  virtual bool operator==(const self_t &a) const;


  /** \brief Returns TRUE if the two elements have different indices on the grid. */
  virtual bool operator!=(const self_t &a) const;

  /**
   * \brief Returns true if two elements belongs from the same Basis.
   */
  bool has_same_basis_of(const self_t &elem) const;
  ///@}


  /**
   * \brief Move the element to the next one with the same property.
   */
  virtual void operator++() = 0;

  /**
   * \brief Move the element to the one specified by <tt>elem_id</tt>.
   *
   * @warning Use this function only if you know what you are doing
   */
  virtual void move_to(const IndexType &elem_id) = 0;


  /**
   * @name Functions for getting information about the element connectivity.
   */
  ///@{
  /**
   * \brief Returns the global dofs of the local (non zero) basis functions
   * on the element.
   *
   * @note The dofs can be filtered invoking the function with the argument @p dof_property.
   * If @p dof_property is equal to DofProperties::active, then no filter is applied.
   *
   * For example:
   * \code
     auto loc_to_glob_all = elem->get_local_to_global(DofProperties::active);
     // loc_to_glob_all[0] is the global id of the first basis function on the element
     // loc_to_glob_all[1] is the global id of the second basis function on the element
     // ...
     auto loc_to_glob_active = elem->get_local_to_global(DofProperties::active);
     // loc_to_glob_active[0] is the global id of the first active basis function on the element
     // loc_to_glob_active[1] is the global id of the second active basis function on the element
     // ...
    \endcode
   *
   */
  SafeSTLVector<Index>
  get_local_to_global(const std::string &dofs_property = DofProperties::active) const;

  /**
   * \brief Returns the patch dofs of the local (non zero) basis functions
   * on the element.
   *
   * @note The dofs can be filtered invoking the function with the argument @p dof_property.
   * If @p dof_property is equal to DofProperties::active, then no filter is applied.
   *
   */
  SafeSTLVector<Index>
  get_local_to_patch(const std::string &dofs_property = DofProperties::active) const;


  SafeSTLVector<Index>
  get_local_dofs(const std::string &dofs_property = DofProperties::active) const;


  /**
   * \brief Returns the number of non zero basis functions with the given
   * @p dofs_property,
   * that are non-zero over the current element.
   */
  Size get_num_basis(const std::string &dofs_property = DofProperties::active) const;
  ///@}


  /**
   * \brief Return a reference to the underlying GridElement.
   *
   * \warning Use this function only if you know what you are doing.
   */
  virtual GridElem &get_grid_element() = 0;

  /**
   * \brief Return a const-reference to the underlying GridElement.
   */
  virtual const GridElem &get_grid_element() const = 0;

  /**
   * \brief Returns the index of the element.
   */
  const IndexType &get_index() const;


  /**
   * \name Functions for retrieving data from the cache.
   */
  ///@{
  /**
   * \brief Returns the basis function <tt>ValueType</tt> data from the cache
   * relative to the <tt>sdim</tt>-dimensional  <tt>s_id</tt>-th sub-element.
   *
   * @note The returned data is filtered by the <tt>dofs_property</tt>.
   */
  template <class ValueType, int sdim = dim_>
  auto
  get_basis_data(const int s_id, const std::string &dofs_property = DofProperties::active) const
  {
    const auto &values_all_elem_dofs = all_sub_elems_cache_.
                                       template get_sub_elem_cache<sdim>(s_id).
    template get_data<ValueType>();

    //--------------------------------------------------------------------------------------
    // filtering the values that correspond to the dofs with the given property --- begin
    SafeSTLVector<Index> dofs_global;
    SafeSTLVector<Index> dofs_local_to_patch;
    SafeSTLVector<Index> dofs_local_to_elem;

    this->basis_->get_spline_space()->get_element_dofs(
      this->get_index(),
      dofs_global,
      dofs_local_to_patch,
      dofs_local_to_elem,
      dofs_property);

    const auto n_filtered_dofs = dofs_local_to_elem.size();
    const auto n_pts = values_all_elem_dofs.get_num_points();

    using VType = typename std::remove_reference<decltype(values_all_elem_dofs)>::type;
    VType values_filtered_elem_dofs(n_filtered_dofs,n_pts);

    int fn = 0;
    for (const auto loc_dof : dofs_local_to_elem)
    {
      const auto values_all_elem_dofs_fn = values_all_elem_dofs.get_function_view(loc_dof);

      const auto values_filtered_elem_dofs_fn = values_filtered_elem_dofs.get_function_view(fn);

      std::copy(values_all_elem_dofs_fn.begin(),
                values_all_elem_dofs_fn.end(),
                values_filtered_elem_dofs_fn.begin());

      ++fn;
    }
    // filtering the values that correspond to the dofs with the given property --- end
    //--------------------------------------------------------------------------------------

    return values_filtered_elem_dofs;
  }

  /**
   * \brief Returns the basis function <tt>ValueType</tt> data from the cache
   * relative to the <tt>dim</tt>-dimensional element (i.e. the element itself).
   *
   * @note The returned data is filtered by the <tt>dofs_property</tt>.
   */
  template <class ValueType>
  auto
  get_basis_data_element(const std::string &dofs_property = DofProperties::active) const
  {
    return this->template get_basis_data<ValueType,dim_>(0,dofs_property);
  }

  template <class ValueType, int sdim = dim_>
  auto
  linear_combination(const SafeSTLVector<Real> &loc_coefs,
                     const int s_id,
                     const std::string &dofs_property) const
  {
    const auto &basis_values =
      this->template get_basis_data<ValueType,sdim>(s_id,dofs_property);
    return basis_values.evaluate_linear_combination(loc_coefs) ;
  }

  /**
   * \brief Returns the <b>basis values</b> evaluated at the evaluation points on the
   * <dim>-dimensional element (i.e. the element itself).
   *
   * @note Before call this function, the BasisElement cache must be filled with
   * the appropriate data.
   */
  ValueTable<Value>
  get_element_values(const std::string &dofs_property = DofProperties::active) const;

  /**
   * \brief Returns the <b>basis gradients</b> evaluated at the evaluation points on the
   * <dim>-dimensional element (i.e. the element itself).
   *
   * @note Before call this function, the BasisElement cache must be filled with
   * the appropriate data.
   */
  ValueTable<Derivative<1> >
  get_element_gradients(const std::string &dofs_property = DofProperties::active) const;

  /**
   * \brief Returns the <b>basis hessians</b> evaluated at the evaluation points on the
   * <dim>-dimensional element (i.e. the element itself).
   *
   * @note Before call this function, the BasisElement cache must be filled with
   * the appropriate data.
   */
  ValueTable<Derivative<2> >
  get_element_hessians(const std::string &dofs_property = DofProperties::active) const;

  /**
   * \brief Returns the <b>basis divergences</b> evaluated at the evaluation points on the
   * <dim>-dimensional element (i.e. the element itself).
   *
   * \note Before call this function, the BasisElement cache must be filled with
   * the appropriate data.
   *
   * \warning The divergence is defined only if ...
   *
   * \todo Complete the documentation (03 Nov 2015)
   */
  ValueTable<Div>
  get_element_divergences(const std::string &dofs_property = DofProperties::active) const;

  /**
   * \brief Returns the <tt>sdim</tt>-dimensional <tt>s_id</tt>-th sub-element measure
   * multiplied by the weights of the quadrature.
   */
  template <int sdim>
  ValueVector<Real> get_w_measures(const int s_id) const;

  /**
   * \brief Returns the <tt>dim</tt>-dimensional element measure
   * multiplied by the weights of the quadrature.
   */
  ValueVector<Real> get_element_w_measures() const;
  ///@}

  /**
   * @name Functions for the basis evaluations without the use of the cache.
   */
  ///@{
  /**
   * \brief Returns a ValueTable with the quantity specified by the template parameter @p ValueType,
   * computed for all local basis function,
   * at each point (in the unit domain) specified by the input argument <tt>points</tt>.
   * @note This function does not use the cache and therefore can be called any time without
   * needing to pre-call init_cache()/fill_cache().
   * @warning The evaluation <tt>points</tt> must belong to the unit hypercube
   * \f$ [0,1]^{\text{dim}} \f$ otherwise, in Debug mode, an assertion will be raised.
   */
  template <class ValueType>
  auto
  evaluate_basis_at_points(
    const std::shared_ptr<const Quadrature<dim_>> &points,
    const std::string &dofs_property)
  {
    auto elem_handler = this->basis_->create_cache_handler();
    elem_handler->template set_flags<dim_>(ValueType::flag);
    elem_handler->init_element_cache(*this,points);
    elem_handler->fill_element_cache(*this);

    return this->template get_basis_data<ValueType,dim_>(0,dofs_property);
  }

  ///@}


  /**
   * \brief Computes and returns the <b>local mass-matrix</b> for the basis function
   * with the given <tt>dofs_property</tt>, i.e. the matrix \f$M\f$ in which
   * its (i,j) entry is:
   * \f$ m_{ij} = \int_{\omega} \varphi_i \cdot \varphi_j \; d\omega \f$
   *
   * The integration domain \f$ \omega \f$ depends on the value of the template parameter <tt>sdim</tt>
   * and on the input argument <tt>s\_id</tt>.
   * If <tt>sdim == dim </tt> the unique valid value for <tt>s_id</tt> is <tt>0</tt> and
   * the integration is made on the element \f$ \Omega^e \f$,
   * otherwise the integration is perfomed in the sub-elements of topological
   * dimension <tt>sdim</tt>, identified by <tt>s_id</tt>
   *
   * \note If the <tt>dofs_property</tt> is omitted, then the mass matrix refers to
   * all basis function that have support on the element.
   *
   * \pre In order to call this function the following quantities must be available
   * in the element's cache:
   * - basis values
   * - weights multiplied by element's measure (w_measures)
   *
   */
  template <int sdim>
  DenseMatrix
  integrate_u_v(const int s_id,
                const PropId &dofs_property = DofProperties::active);

  virtual
  DenseMatrix
  integrate_u_v_sum_factorization_impl(const TopologyVariants<dim_> &topology,
                                       const int s_id,
                                       const PropId &dofs_property = DofProperties::active)
  {
    AssertThrow(false,ExcMessage("This function must be implemented in a derived class."));
    DenseMatrix mat;
    return mat;
  }

  virtual
  DenseMatrix
  integrate_gradu_gradv_sum_factorization_impl(const TopologyVariants<dim_> &topology,
                                               const int s_id,
                                               const PropId &dofs_property = DofProperties::active)
  {
    AssertThrow(false,ExcMessage("This function must be implemented in a derived class."));
    DenseMatrix mat;
    return mat;
  }


  /**
   * \brief Computes and returns the <b>local stiffness-matrix</b> (of the Poisson problem)
   * for the basis function
   * with the given <tt>dofs_property</tt>, i.e. the matrix \f$M\f$ in which
   * its (i,j) entry is:
   * \f$ m_{ij} = \int_{\omega} \nabla \varphi_i \cdot \nabla \varphi_j \; d\omega \f$
   *
   * The integration domain \f$ \omega \f$ depends on the value of the template parameter <tt>sdim</tt>
   * and on the input argument <tt>s\_id</tt>.
   * If <tt>sdim == dim </tt> the unique valid value for <tt>s_id</tt> is <tt>0</tt> and
   * the integration is made on the element \f$ \Omega^e \f$,
   * otherwise the integration is perfomed in the sub-elements of topological
   * dimension <tt>sdim</tt>, identified by <tt>s_id</tt>
   *
   * \note If the <tt>dofs_property</tt> is omitted, then the stiffness matrix refers to
   * all basis function that have support on the element.
   *
   * \pre In order to call this function the following quantities must be available
   * in the element's cache:
   * - basis gradients
   * - weights multiplied by element's measure (w_measures)
   *
   */
  template <int sdim>
  DenseMatrix
  integrate_gradu_gradv(const int s_id,
                        const PropId &dofs_property = DofProperties::active);

  /**
   * \brief Computes and returns the vector \f$R\f$ in which
   * its (i) entry is:
   * \f$ R_{i} = \int_{\omega} \nabla \varphi_i \cdot f \; d\omega \f$
   *
   * The integration domain \f$ \omega \f$ depends on the value of the template parameter <tt>sdim</tt>
   * and on the input argument <tt>s\_id</tt>.
   * If <tt>sdim == dim </tt> the unique valid value for <tt>s_id</tt> is <tt>0</tt> and
   * the integration is made on the element \f$ \Omega^e \f$,
   * otherwise the integration is perfomed in the sub-elements of topological
   * dimension <tt>sdim</tt>, identified by <tt>s_id</tt>
   *
   *
   * \note If the <tt>dofs_property</tt> is omitted, then the stiffness matrix refers to
   * all basis function that have support on the element.
   *
   * \pre In order to call this function the following quantities must be available
   * in the element's cache:
   * - basis values
   * - weights multiplied by element's measure (w_measures)
   *
   */
  template <int sdim>
  DenseVector
  integrate_u_func(const ValueVector<Value> &func_at_points,
                   const int s_id,
                   const PropId &dofs_property = DofProperties::active);


  /**
   * \name Function for printing internal information.
   */
  ///@{
  /**
   * \brief Prints some information about the object.
   * @note Mostly used for testing and debugging purposes.
   */
  virtual void print_info(LogStream &out) const;

  /**
   * \brief Prints some information about the local cache inside the object.
   * @note Mostly used for testing and debugging purposes.
   */
  virtual void print_cache_info(LogStream &out) const;
  ///@}

private:

  /**
   * \brief Basis upon which the element refers from.
   */
  std::shared_ptr<Bs> basis_;

public:
  using _Value =  basis_element::_Value;
  using _Gradient = basis_element::_Gradient;
  using _Hessian = basis_element::_Hessian;
  using _Divergence = basis_element::_Divergence;

protected:

  using CType = boost::fusion::map<
                boost::fusion::pair<     _Value,DataWithFlagStatus<ValueTable<Value>>>,
                boost::fusion::pair<  _Gradient,DataWithFlagStatus<ValueTable<Derivative<1>>>>,
                boost::fusion::pair<   _Hessian,DataWithFlagStatus<ValueTable<Derivative<2>>>>,
                boost::fusion::pair<_Divergence,DataWithFlagStatus<ValueTable<Div>>>
                >;

public:
  using Cache = BasisValuesCache<dim_,CType>;

  using CacheType = AllSubElementsCache<Cache>;

protected:
  /** The local (element and sub-element(s)) cache. */
  CacheType all_sub_elems_cache_;

public:

  /**
   * \brief Returns the Basis upon which the element refers from.
   */
  std::shared_ptr<Bs> get_basis() const;




  friend class BasisHandler<dim_,codim_,range_,rank_>;


  /**
   * Return TRUE if the element index is referring to a valid element in the Grid.
   *
   * @note An element with non valido position can happens when we use the ++ operator
   * on an element that is the last in the list.
   */
  bool has_valid_position() const
  {
    return this->get_grid_element().has_valid_position();
  }


  /**
   * Returns the number of 1D spline functions in the space component <tt>comp</tt>,
   * along each coordinate direction.
   *
   * @note The returned numbers refer to all 1D splines in the element, without distinction
   * between their property (if any).
   */
  TensorSize<dim_> get_num_splines_1D(const int comp) const;
};


IGA_NAMESPACE_CLOSE



#endif // #ifndef __BASIS_ELEMENT_H_

