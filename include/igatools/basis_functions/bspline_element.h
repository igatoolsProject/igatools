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

#ifndef BSPLINE_ELEMENT_H_
#define BSPLINE_ELEMENT_H_

#include <igatools/base/config.h>
#include <igatools/basis_functions/reference_element.h>

#include <igatools/linear_algebra/dense_matrix.h>
#include <igatools/basis_functions/bernstein_basis.h>
#include <igatools/basis_functions/values1d_const_view.h>

IGA_NAMESPACE_OPEN

template <int dim, int range, int rank> class BSpline;
template <int dim, int range, int rank> class BSplineElementHandler;
template <class Accessor> class GridIterator;

/**
 * See module on \ref accessors_iterators for a general overview.
 * @ingroup elements
 */
template <int dim, int range, int rank>
class BSplineElement :
  public ReferenceElement<dim,range,rank>
{
private:
  using self_t = BSplineElement<dim,range,rank>;
  using parent_t = ReferenceElement<dim,range,rank>;

public:

  /** Type required by the GridIterator templated iterator */
  using ContainerType = const BSpline<dim, range, rank> ;

  /** Type required for the generic algorithm on the spaces (plots??) */
  using Basis = BSpline<dim, range, rank> ;


  using GridType = Grid<dim>;
  using IndexType = typename GridType::IndexType;
  using List = typename GridType::List;
  using ListIt = typename GridType::ListIt;

  using GridElem = GridElement<dim>;

public:
  template <int order>
  using Derivative = typename parent_t::template Derivative<order>;
  using typename parent_t::Point;
  using typename parent_t::Value;

  /** @name Constructors */
  ///@{
public:

  /**
   * Default constructor.
   */
  BSplineElement() = delete;


  /**
   * Constructs an accessor to element number index of a
   * BSpline space.
   */
  BSplineElement(const std::shared_ptr<ContainerType> space,
                 const ListIt &index,
                 const PropId &prop = ElementProperties::active);

  /**
   * Copy constructor. Not allowed to be used.
   */
  BSplineElement(const self_t &elem) = delete;

  /**
   * Move constructor.
   */
  BSplineElement(self_t &&elem) = default;

  /**
   * Destructor.
   */
  virtual ~BSplineElement() = default;
  ///@}

  /** @name Assignment operators */
  ///@{
  /**
   * Copy assignment operator. Not allowed to be used.
   */
  self_t &operator=(const self_t &elem) = delete;

  /**
   * Move assignment operator.
   */
  self_t &operator=(self_t &&elem) = default;
  ///@}





  /**
   * Returns the BSpline upon which the element is defined.
   */
  std::shared_ptr<const Basis> get_bspline_basis() const;

private:

  template <class Accessor> friend class GridIterator;
  friend class BSplineElementHandler<dim, range, rank>;


  template<class T>
  using ComponentContainer =
    typename Basis::template ComponentContainer<T>;


  using Splines1D = BasisValues1d;
  using Splines1DTable = ComponentContainer<SafeSTLArray<Splines1D,dim>>;
  using AllSplines1DTable = SafeSTLArray<SafeSTLVector<Splines1DTable>,dim+1>;
  AllSplines1DTable all_splines_1D_table_;


public:
#if 0
  ComponentContainer<SafeSTLArray<ValueTable<Real>,dim> >
  evaluate_univariate_derivatives_at_points(const int deriv_order,
                                            const Quadrature<dim> &pts) const
  {
    ComponentContainer<SafeSTLArray<ValueTable<Real>,dim> > values;

    Assert(false,ExcNotImplemented());

    return values;
  }
#endif

  virtual void print_cache_info(LogStream &out) const override final;



private:
  std::unique_ptr<GridElem> grid_elem_;


public:
  /**
   * Return a reference to the GridElement.
   */
  GridElem &get_grid_element() override final;

  /**
   * Return a const-reference to the GridElement.
   */
  const GridElem &get_grid_element() const override final;

  virtual void operator++() override;

  /**
   * Move the element to the one specified by <tt>elem_id</tt>.
   *
   * @warning Use this function only if you know what you are doing
   */
  virtual void move_to(const IndexType &elem_id) override final;


  virtual void print_info(LogStream &out) const override final;

};

IGA_NAMESPACE_CLOSE

#endif
