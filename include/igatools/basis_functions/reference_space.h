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
#include <igatools/functions/function_element.h>
#include <igatools/utils/cartesian_product_array.h>
#include <igatools/basis_functions/space.h>
#include <igatools/geometry/grid.h>
#include <igatools/basis_functions/space_element.h>
#include <igatools/basis_functions/spline_space.h>
#include <igatools/basis_functions/dof_distribution.h>
#include <igatools/utils/multi_array_utils.h>

IGA_NAMESPACE_OPEN

template<Transformation,int, int> class PushForward;

template <int, int, int ,int,Transformation> class PhysicalSpace;

template <int, int, int> class ReferenceElement;
template <int,int,int> class ReferenceElementHandler;

template <int, int, int> class BSplineSpace;
template <int, int, int> class NURBSSpace;


template <int,int,int> class DofDistribution;


/**
 * @brief Base abstract class for reference spaces (i.e BSplineSpace and NURBSSpace).
 *
 * @ingroup containers
 * @ingroup serializable
 */
template<int dim_, int range_ = 1, int rank_ = 1>
class ReferenceSpace :
  public Space<dim_,0,range_,rank_,Transformation::h_grad>
{
public:
  using base_t = Space<dim_,0,range_,rank_,Transformation::h_grad>;
  using self_t = ReferenceSpace<dim_,range_,rank_>;

  static const int dim       = dim_;
  static const int codim     = 0;
  static const int space_dim = dim_;
  static const int range     = range_;
  static const int rank      = rank_;
  static const bool is_physical_space = false;

  /**
   * See documentation in \ref Space
   *
   * @see Space
   */
  using PushForwardType = PushForward<Transformation::h_grad, dim, codim>;

  using RefSpace = ReferenceSpace<dim_,range_,rank_>;

  using Func = Function<dim, 0, range, rank>;

  template <int order>
  using Derivative = typename Func::template Derivative<order>;
  using Point = typename Func::Point;
  using Value = typename Func::Value;
  using Div   = typename Func::Div;
  using RefPoint = Point;


  using GridType = Grid<dim>;



  /** Type for the element accessor. */
  using ElementAccessor = ReferenceElement<dim,range,rank>;

  /** Type for iterator over the elements.  */
  using ElementIterator = GridIterator<ElementAccessor>;

  using ElementHandler = ReferenceElementHandler<dim_, range_, rank_>;

  using SpaceData = SplineSpace<dim_,range_,rank_>;

  using Degrees = typename SpaceData::Degrees;
  using Multiplicity = typename SpaceData::Multiplicity;
  using EndBehaviour = typename SpaceData::EndBehaviour;
  using Periodicity = typename SpaceData::Periodicity;

  using KnotsTable = typename SpaceData::KnotsTable;
  using DegreeTable = typename SpaceData::DegreeTable;
  using MultiplicityTable = typename SpaceData::MultiplicityTable;
  using TensorSizeTable = typename SpaceData::TensorSizeTable;
  using PeriodicityTable = typename SpaceData::PeriodicityTable;
  using EndBehaviourTable = typename SpaceData::EndBehaviourTable;

  template <class T>
  using ComponentContainer = typename SpaceData::template ComponentContainer<T>;

  using ComponentMap = typename SpaceData::template ComponentContainer<int>::ComponentMap;

  static const auto n_components = SpaceData::n_components;

protected:
  /**
   * Default constructor. It does nothing but it is needed for the serialization
   * mechanism.
   */
  ReferenceSpace() = default;


public:
  virtual ~ReferenceSpace() = default;


  template <int sdim>
  using SubGridMap = typename GridType::template SubGridMap<sdim>;

  template <int k>
  using InterSpaceMap = SafeSTLVector<Index>;

  template <int k>
  using SubRefSpace = ReferenceSpace<k, range, rank>;

  template <int k>
  using SubSpace = PhysicalSpace<k,range,rank, dim-k, Transformation::h_grad>;

  virtual bool is_bspline() const = 0;

  /**
   * Returns the degree of the BSpline space for each component and for each coordinate direction.
   * \return The degree of the BSpline space for each component and for each coordinate direction.
   * The first index of the returned object is the component id, the second index is the direction id.
   */
  virtual const DegreeTable &get_degree_table() const = 0;


  /**
   * Return the maximum value of the polynomial degree, for each component, for each direction;
   */
  virtual int get_max_degree() const override final;





  /**
   * Returns a const reference to the end behaviour table of the BSpline space.
   */
  virtual const EndBehaviourTable &get_end_behaviour_table() const = 0;


  virtual const PeriodicityTable &get_periodicity() const = 0;







  template<int k>
  std::shared_ptr< SubRefSpace<k> >
  get_ref_sub_space(const int s_id,
                    InterSpaceMap<k> &dof_map,
                    std::shared_ptr<Grid<k>> sub_grid = nullptr) const;

  template<int k>
  std::shared_ptr<SubSpace<k> >
  get_sub_space(const int s_id, InterSpaceMap<k> &dof_map,
                std::shared_ptr<Grid<k>> sub_grid,
                SubGridMap<k> &elem_map) const;


  virtual std::unique_ptr<ReferenceElement<dim_,range_,rank_> >
  create_ref_element(const typename base_t::ListIt &index, const PropId &property) const = 0;

protected:

  std::shared_ptr<const RefSpace> ref_space_previous_refinement_ = nullptr;



#ifdef MESH_REFINEMENT

public:
  std::shared_ptr<const base_t> get_space_previous_refinement() const override final
  {
    return ref_space_previous_refinement_;
  }
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
  serialize(Archive &ar)
  {
    using std::to_string;
    const std::string base_name = "Space_" +
                                  to_string(dim_) + "_" +
                                  to_string(0) + "_" +
                                  to_string(range_) + "_" +
                                  to_string(rank_) + "_hgrad";

    ar &make_nvp(base_name,base_class<base_t>(this));

#ifdef MESH_REFINEMENT
    auto tmp = std::const_pointer_cast<RefSpace>(ref_space_previous_refinement_);
    ar &make_nvp("ref_space_previous_refinement_",tmp);
    ref_space_previous_refinement_ = std::const_pointer_cast<const RefSpace>(tmp);
#endif
  }

  ///@}
#endif // SERIALIZATION
};

IGA_NAMESPACE_CLOSE

#endif // #ifndef REFERENCE_SPACE_H_
