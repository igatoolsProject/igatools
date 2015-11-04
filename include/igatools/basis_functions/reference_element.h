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

#ifndef REFERENCE_ELEMENT_H_
#define REFERENCE_ELEMENT_H_

#include <igatools/base/config.h>
#include <igatools/basis_functions/space_element.h>
#include <igatools/basis_functions/reference_element_handler.h>

IGA_NAMESPACE_OPEN

template <int, int, int> class ReferenceSpaceBasis;

/**
 *
 * @ingroup elements
 */
template <int dim, int range, int rank>
class ReferenceElement : public SpaceElement<dim,0,range,rank>
{
public:

  /** Type required by the GridForwardIterator templated iterator */
  using ContainerType = const ReferenceSpaceBasis<dim,range,rank> ;

  using Basis = ReferenceSpaceBasis<dim,range,rank>;
  using ConstSpace = const ReferenceSpaceBasis<dim,range,rank>;

  using parent_t = SpaceElement<dim,0,range,rank>;

  using RefPoint = typename Basis::RefPoint;
  using Point = typename Basis::Point;
  using Value = typename Basis::Value;

  template <int order>
  using Derivative = typename Basis::template Derivative<order>;

  using Div = typename Basis::Div;


  using GridType = Grid<dim>;
  using IndexType = typename GridType::IndexType;
  using List = typename GridType::List;
  using ListIt = typename GridType::ListIt;

  using GridElem = GridElement<dim>;

public:
  ReferenceElement() = delete;

  /**
   * Copy constructor. Not allowed to be used.
   */
  ReferenceElement(const ReferenceElement<dim,range,rank> &elem) = delete;

  /**
   * Constructs an accessor to element number index of a
   * ReferenceSpaceBasis basis.
   */
  ReferenceElement(const std::shared_ptr<ConstSpace> &basis,
                   const ListIt &index,
                   const PropId &prop = ElementProperties::active);


  virtual ~ReferenceElement() = default;

  /**
   * Returns the <tt>k</tt> dimensional j-th sub-element measure
   * multiplied by the weights of the quadrature on the unit element.
   */
  template <int sdim>
  ValueVector<Real> get_w_measures(const int s_id) const
  {
    return this->get_grid_element().template get_weights<sdim>(s_id);
  }

  /**
   * Returns the gradient determinant of the identity map at the dilated quadrature points.
   */
  ValueVector<Real> get_element_w_measures() const;


//    using OffsetTable = typename Basis::template ComponentContainer<int>;
  using OffsetTable = SafeSTLArray<int,Basis::n_components+1>;

  using TensorSizeTable = typename Basis::TensorSizeTable;

protected:

  /** Number of scalar basis functions along each direction, for all space components. */
  TensorSizeTable n_basis_direction_;

  /**
   * Offset of the scalar basis functions across the different components.
   *
   * @note The last entry of the array contains the total number of scalar basis functions.
   */
  OffsetTable comp_offset_;

  using Indexer = CartesianProductIndexer<dim>;
  using IndexerPtr = std::shared_ptr<Indexer>;
  using IndexerPtrTable = typename Basis::template ComponentContainer<IndexerPtr>;

  /** Hash table for fast conversion between flat-to-tensor basis function ids. */
  IndexerPtrTable basis_functions_indexer_;

  std::shared_ptr<const Basis> basis_;

public:
  using parent_t::get_num_basis;


  /**
   * Returns the basis function ID offset between the different components.
   */
  OffsetTable get_basis_offset() const;

  /**
   * Number of non-zero scalar basis functions associated
   * with the i-th space component on the element.
   * This makes sense as a reference B-spline space
   * is only allowed to be of the cartesian product type
   * V = V1 x V2 x ... X Vn.
   */
  int get_num_basis_comp(const int i) const;

  virtual void print_info(LogStream &out) const override;


#if 0
  std::shared_ptr<const Basis> get_ig_space() const;
#endif
};

IGA_NAMESPACE_CLOSE

#endif // #ifndef REFERENCE_ELEMENT_H_
