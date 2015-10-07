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


#ifndef NURBS_ELEMENT_H_
#define NURBS_ELEMENT_H_

#include <igatools/base/config.h>

#ifdef NURBS

#include <igatools/basis_functions/reference_element.h>


//#include <igatools/linear_algebra/dense_matrix.h>
//#include <igatools/basis_functions/bernstein_basis.h>
//#include <igatools/basis_functions/bspline_element_scalar_evaluator.h>



IGA_NAMESPACE_OPEN




template <int dim, int range, int rank> class NURBSSpace;
template <int,int,int> class NURBSElementHandler;
template <class Accessor> class GridIterator;

/**
 * @brief NURBS element
 *
 * See module on \ref accessors_iterators for a general overview.
 * @ingroup elements
 * @ingroup serialization
 */
template <int dim, int range, int rank>
class NURBSElement :
  public ReferenceElement<dim,range,rank>
{
private:
  using self_t = NURBSElement<dim,range,rank>;
  using parent_t = ReferenceElement<dim,range,rank>;

public:

  /** Type required by the GridIterator templated iterator */
  using ContainerType = const NURBSSpace<dim, range, rank> ;

  /** Type required for the generic algorithm on the spaces (plots??) */
  using Space = NURBSSpace<dim, range, rank> ;


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
protected:
  /**
   * Default constructor. It does nothing but it is needed for the
   * <a href="http://www.boost.org/doc/libs/release/libs/serialization/">boost::serialization</a>
   * mechanism.
   */
  NURBSElement() = default;

public:
  /**
   * Constructs an accessor to element number index of a
   * BsplineSpace space.
   */
  NURBSElement(const std::shared_ptr<ContainerType> space,
               const ListIt &index,
               const PropId &prop = ElementProperties::active);

#if 0
  /**
   * Copy constructor.
   * It can be used with different copy policies (i.e. deep copy or shallow copy).
   * The default behaviour (i.e. using the proper interface of a classic copy constructor)
   * uses the deep copy.
   */
  NURBSElement(const self_t &elem,
               const CopyPolicy &copy_policy = CopyPolicy::deep);
//*/
#endif

  /**
   * Move constructor.
   */
  NURBSElement(self_t &&elem) = default;

  /**
   * Destructor.
   */
  virtual ~NURBSElement() = default;
  ///@}

  /** @name Assignment operators */
  ///@{
  /**
   * Copy assignment operator.
   * @note Creates a new element cache, but it shares
   * the one dimensional cache with the copied element.
   */
  self_t &operator=(const self_t &elem) = default;

  /**
   * Move assignment operator.
   */
  self_t &operator=(self_t &&elem) = default;
  ///@}

#if 0
  /** @name Functions/operators for moving the element in the NURBSSpace.*/
  ///@{
  /**
   * Sets the index of the element using the flatten representation.
   * @note This function also updates the index for the tensor representation.
   * @warning This may be a dangerous function, be careful when using it
   * as it is easy to use incorrectly. Only use it if you know what you
   * are doing.
   */
  virtual void move_to(const Index flat_index) override final;
  ///@}
#endif

  /**
   * Returns the NURBSSpace in which the NURBSElement is defined.
   */
  std::shared_ptr<const Space> get_nurbs_space() const;

private:

  using SpSpace = typename Space::SpSpace;

  typename SpSpace::ElementAccessor bspline_elem_;

  using WeightFunction = typename Space::WeightFunction;
  using WeightElem = typename WeightFunction::ConstElementAccessor;
  std::unique_ptr<WeightElem> weight_elem_;

public:

  friend class NURBSElementHandler<dim, range, rank>;


  /**
   * Return a reference to the GridElement.
   */
  GridElem &get_grid_element() override final
  {
    return bspline_elem_.get_grid_element();
  }

  /**
   * Return a const-reference to the GridElement.
   */
  const GridElem &get_grid_element() const override final
  {
    return bspline_elem_.get_grid_element();
  }


  virtual void operator++() override final
  {
    ++(*weight_elem_);
    ++bspline_elem_;
    Assert(false,ExcNotImplemented());
  }

  /**
   * Move the element to the one specified by <tt>elem_id</tt>.
   *
   * @warning Use this function only if you know what you are doing
   */
  virtual void move_to(const IndexType &elem_id) override final
  {
    AssertThrow(false,ExcNotImplemented());
  }

};

IGA_NAMESPACE_CLOSE

#endif // #ifdef NURBS


#endif // end of #ifndef NURBS_ELEMENT_H_



