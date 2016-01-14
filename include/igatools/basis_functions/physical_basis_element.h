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


#ifndef PHYSICAL_SPACE_ELEMENT_H
#define PHYSICAL_SPACE_ELEMENT_H

#include <igatools/base/config.h>

#include <igatools/base/quadrature.h>
#include <igatools/basis_functions/bspline_element.h>
#include <igatools/basis_functions/nurbs_element.h>
#include <igatools/basis_functions/physical_basis.h>
#include <igatools/geometry/push_forward.h>



IGA_NAMESPACE_OPEN

template <class Accessor> class GridIterator;

/**
 *
 * @ingroup elements
 */
template<int dim_,int range_,int rank_,int codim_>
class PhysicalBasisElement
  :
  public BasisElement<dim_,codim_,range_,rank_>
{
  using self_t = PhysicalBasisElement<dim_,range_,rank_,codim_>;
  using parent_t = BasisElement<dim_,codim_,range_,rank_>;

public :
  using PhysSpace = PhysicalBasis<dim_,range_,rank_,codim_>;
  /** Type required by the GridIterator templated iterator */
  using ContainerType = const PhysSpace;

  using Basis = PhysSpace;
  using RefBasis = typename PhysSpace::RefBasis;
  using PushFwd = typename PhysSpace::PushFwd;
//    using RefElemAccessor = BasisElement<RefBasis::dim,0,RefBasis::range,RefBasis::rank,Transformation::h_grad>;
  using RefElemAccessor = ReferenceBasisElement<RefBasis::dim,RefBasis::range,RefBasis::rank>;

  using PhysDomain = Domain<dim_, codim_>;
  using PhysDomainElem = DomainElement<dim_, codim_>;

  static const auto dim = PushFwd::dim;
  static const auto space_dim = PushFwd::space_dim;
  static const auto codim = PushFwd::codim;

  using PhysPoint = typename Basis::Point;

  using GridType = Grid<dim>;
  using IndexType = typename GridType::IndexType;
  using List = typename GridType::List;
  using ListIt = typename GridType::ListIt;


  /**
   * @name Constructors
   */
  ///@{
public:

  /**
   * Default constructor. Not allowed to be used.
   */
  PhysicalBasisElement() = delete;

  PhysicalBasisElement(const std::shared_ptr<ContainerType> &space,
                       GridIterator<RefElemAccessor> &&ref_space_element,
                       GridIterator<PhysDomainElem> &&phys_domain_element);


  /**
   * Copy constructor. Not allowed to be used.
   */
  PhysicalBasisElement(const self_t &in) = delete;


  /**
   * Move constructor.
   */
  PhysicalBasisElement(self_t &&in) = default;

  /**
   * Destructor.
   */
  virtual ~PhysicalBasisElement() = default;

  ///@}

  /**
   * @name Assignment operators
   */
  ///@{
  /**
   * Copy assignment operator. Not allowed to be used.
   */
  self_t &
  operator=(const self_t &in) = delete;

  /**
   * Move assignment operator.
   */
  self_t &
  operator=(self_t &&in) = default;

  ///@}



  /**
   * @name Getting quantities that are geometry-related
   */
  ///@{
  /**
   * Returns the <tt>k</tt> dimensional j-th sub-element measure
   * multiplied by the weights of the quadrature.
   */
  template <int k>
  const ValueVector<Real> get_w_measures(const int j) const;

  /**
   * Returns the gradient determinant of the map at the dilated quadrature points.
   */
  template <int k>
  const ValueVector<Real> &get_measures(const int j) const;

  const ValueVector<Real> get_element_w_measures() const;

  template <int k = dim>
  const ValueVector<PhysPoint> get_points(const int j = 0) const;

  const ValueVector<PhysPoint> get_element_points() const;

  template<int k>
  const ValueVector<Points<space_dim> > &
  get_boundary_normals(const int s_id) const;


  /**
   * Prints internal information about the BSplineElementAccessor.
   * Its main use is for testing and debugging.
   */
  void print_info(LogStream &out) const override final;

  virtual void print_cache_info(LogStream &out) const override final;
  ///@}

public:

  /**
   * Return a reference to the GridElement.
   */
  virtual GridElement<dim_> &get_grid_element() override final;

  /**
   * Return a const-reference to the GridElement.
   */
  virtual const GridElement<dim_> &get_grid_element() const override final;

  /**
   * Return a const reference of the reference space element.
   */
  const RefElemAccessor &get_ref_space_element() const;

  /**
   * Return a non-const reference of the reference space element.
   */
  RefElemAccessor &get_ref_space_element();


  /**
   * Return a const reference of the DomainElement.
   */
  const PhysDomainElem &get_domain_element() const;

  /**
   * Return a non-const reference of the DomainElement.
   */
  PhysDomainElem &get_domain_element();

public:
  using parent_t::get_num_basis;



  /** Returns the index of the element. */
  IndexType get_index() const;


  /** @name Functions/operators for moving the element in the Grid.*/
  ///@{
  virtual void operator++() override final;

  /**
   * Move the element to the one specified by <tt>elem_id</tt>.
   *
   * @warning Use this function only if you know what you are doing
   */
  virtual void move_to(const IndexType &elem_id) override final;
  ///@}


  /**
   * @name Comparison operators.
   *
   * @brief The comparison operators compares the <em>position</em> of the element in the grid.
   *
   * @warning To be comparable, two BasisElement objects must be defined on the same space
   * (and therefore on the same grid),
   * otherwise an assertion will be raised (in Debug mode).
   */
  ///@{
  /** Returns TRUE if the two elements have the same index on the grid. */
  virtual bool operator==(const parent_t &a) const override final;


  /** Returns TRUE if the two elements have different indices on the grid. */
  virtual bool operator!=(const parent_t &a) const override final;
  ///@}


  std::shared_ptr<const PhysSpace> get_physical_space() const;


private:
  template <class Accessor> friend class GridIteratorBase;
  template <int,int,int,int> friend class PhysicalBasisHandler;

  GridIterator<RefElemAccessor> ref_space_element_;

  GridIterator<PhysDomainElem> phys_domain_element_;

  std::shared_ptr<const PhysSpace> phys_space_;

  /**
   * Returns true if two elements belongs from the same Basis.
   */
  bool is_comparable_with(const self_t &elem) const;

};


IGA_NAMESPACE_CLOSE

#endif
