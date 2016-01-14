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

#ifndef REFERENCE_SPACE_H_
#define REFERENCE_SPACE_H_

#include <igatools/base/config.h>
#include <igatools/base/array_utils.h>
//#include <igatools/functions/function_element.h>
#include <igatools/utils/cartesian_product_array.h>
#include <igatools/basis_functions/basis.h>
#include <igatools/geometry/grid.h>
#include <igatools/basis_functions/basis_element.h>
#include <igatools/basis_functions/spline_space.h>
#include <igatools/basis_functions/dof_distribution.h>
#include <igatools/utils/multi_array_utils.h>

IGA_NAMESPACE_OPEN

template <int, int, int ,int> class PhysicalBasis;

template <int, int, int> class ReferenceBasisElement;
template <int,int,int> class ReferenceBasisHandler;

template <int, int, int> class BSpline;
template <int, int, int> class NURBS;


template <int,int,int> class DofDistribution;


/**
 * @brief Base abstract class for reference spaces (i.e BSpline and NURBS).
 *
 * @ingroup containers
 * @ingroup serializable
 */
template<int dim_, int range_ = 1, int rank_ = 1>
class ReferenceBasis :
  public Basis<dim_,0,range_,rank_>
{
public:
  using base_t = Basis<dim_,0,range_,rank_>;
  using self_t = ReferenceBasis<dim_,range_,rank_>;

  static const int dim       = dim_;
  static const int codim     = 0;
  static const int space_dim = dim_;
  static const int range     = range_;
  static const int rank      = rank_;
  static const bool is_physical_space = false;

  /**
   * See documentation in \ref Basis
   *
   * @see Basis
   */

  using RefBasis = ReferenceBasis<dim_,range_,rank_>;

  template <int order>
  using Derivative = typename base_t::template Derivative<order>;
  using Point = typename base_t::Point;
  using Value = typename base_t::Value;
  using Div   = typename base_t::Div;
  using RefPoint = Point;


  using GridType = Grid<dim>;



  /** Type for the element accessor. */
  using ElementAccessor = ReferenceBasisElement<dim,range,rank>;

  /** Type for iterator over the elements.  */
  using ElementIterator = GridIterator<ElementAccessor>;

  using Handler = ReferenceBasisHandler<dim_, range_, rank_>;

  using SpSpace = SplineSpace<dim_,range_,rank_>;

  using EndBehaviour = typename SpSpace::EndBehaviour;
  using EndBehaviourTable = typename SpSpace::EndBehaviourTable;

  template <class T>
  using ComponentContainer = typename SpSpace::template ComponentContainer<T>;

  using ComponentMap = typename SpSpace::template ComponentContainer<int>::ComponentMap;

  static const auto n_components = SpSpace::n_components;

protected:
  /**
   * Default constructor. It does nothing but it is needed for the serialization
   * mechanism.
   */
  ReferenceBasis() = default;


public:
  virtual ~ReferenceBasis() = default;


  template <int sdim>
  using SubGridMap = typename GridType::template SubGridMap<sdim>;

  template <int k>
  using InterBasisMap = SafeSTLVector<Index>;

  template <int k>
  using SubRefSpace = ReferenceBasis<k, range, rank>;

  template <int k>
  using SubSpace = PhysicalBasis<k,range,rank, dim-k>;

  virtual bool is_bspline() const = 0;





  /**
   * Returns a const reference to the end behaviour table of the BSpline space.
   */
  virtual const EndBehaviourTable &get_end_behaviour_table() const = 0;







  template<int sdim>
  std::shared_ptr<const SubRefSpace<sdim> >
  get_ref_sub_space(const int s_id,
                    InterBasisMap<sdim> &dof_map,
                    const std::shared_ptr<const Grid<sdim>> &sub_grid = nullptr) const;

  template<int sdim>
  std::shared_ptr<const SubSpace<sdim> >
  get_sub_space(const int s_id, InterBasisMap<sdim> &dof_map,
                SubGridMap<sdim> &elem_map) const;


  virtual std::unique_ptr<ReferenceBasisElement<dim_,range_,rank_> >
  create_ref_element_begin(const PropId &property) const = 0;

  virtual std::unique_ptr<ReferenceBasisElement<dim_,range_,rank_> >
  create_ref_element_end(const PropId &property) const = 0;

protected:

  std::shared_ptr<const RefBasis> ref_basis_previous_refinement_ = nullptr;



#ifdef MESH_REFINEMENT

public:
  std::shared_ptr<const self_t> get_basis_previous_refinement() const;

  void create_connection_for_insert_knots(const std::shared_ptr<self_t> &space);


#endif // MESH_REFINEMENT

private:

#ifdef SERIALIZATION
  /**
   * @name Functions needed for serialization
   */
  ///@{
  friend class cereal::access;

  template<class Archive>
  void
  serialize(Archive &ar);
  ///@}
#endif // SERIALIZATION
};

IGA_NAMESPACE_CLOSE

#endif // #ifndef REFERENCE_SPACE_H_
