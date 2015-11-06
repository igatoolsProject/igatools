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

#ifndef __BSPLINE_SPACE_H_
#define __BSPLINE_SPACE_H_

#include <igatools/base/config.h>
#include <igatools/base/logstream.h>

#include <igatools/basis_functions/reference_space_basis.h>
#include <igatools/basis_functions/bernstein_extraction.h>
#include <igatools/geometry/domain.h>
#include <igatools/basis_functions/physical_space_basis.h>

IGA_NAMESPACE_OPEN


template <int, int, int> class BSplineElement;
template <int, int, int> class BSplineElementHandler;
/**
 * Multivariate (tensor product) scalar, vector or k-tensor
 * valued B-spline space.
 * This object can be thought of as providing the
 * B-spline basis functions of a spline space.
 * The space is determined by:
 * - the knot vectors (implemented in the class Grid)
 * - the multiplicity vectors
 * - and the degree
 *
 * BSpline allows the use of different multiplicity vectors
 * and degrees for each direction and for each component.
 *
 * \section const Constructing a BSpline
 * Similar to the mechanism use in Grid we use
 * the create technique to create a smartpointer for each constructor
 * the class provides.
 * @todo enter a glossary for create idiom technique and refer from here
 * \code
 * auto space = BSpline<dim>::create();
 * \endcode
 *
 * \section eval Evaluating basis function
 * Basis function are evaluated iterating on the elements.
 * Similarly we can obtain other information such as the local to global.
 *
 * \section dofs Degrees of freedom (dofs)
 * Each basis function is assigned a degree of freedom (an integer),
 * the order of this assignment is done following some prescribed
 * recipe or policy whose knowledge is not Really needed by most users.
 * They are internally stored in a grid-like multiarray container,
 * called the index space.
 * It works together with the index_space_offset which for each element
 * provides a view of the index
 * space to know which
 * are the non-zero basis function on each element, what we generally
 * refer to as the local to global mapping.
 *
 * \section bezier Storage of the basis functions (Bezier Extraction)
 * The basis functions on each element are stored in the Bspline space
 * as the 1D Bezier extraction operator.
 * When they need to be evaluated the operator applied to the
 * Berenstein polynomials,
 * allows to compute the values at the quadrature points.
 *
 * @todo write a module about cache optimization and handling.
 *
 * \section hom_range Optimizing homogeneous range type vector spaces
 *
 * \author martinelli, 2012, 2013, 2014
 * \author pauletti, 2012, 2013, 2014
 *
 *
 * @ingroup containers
 * @ingroup serializable
 */
template<int dim_, int range_ = 1, int rank_ = 1>
class BSpline :
  public ReferenceSpaceBasis<dim_, range_, rank_>
{
private:
  using BaseSpace = ReferenceSpaceBasis<dim_, range_, rank_>;

  /** Type for current class. */
  using self_t = BSpline<dim_,range_,rank_>;

public:
  /** see documentation in \ref Basis */

  using GridType = Grid<dim_>;
  using ElementHandler = BSplineElementHandler<dim_, range_, rank_>;


  static const int dim       = dim_;
  static const int codim     = 0;
  static const int space_dim = dim_;
  static const int range     = range_;
  static const int rank      = rank_;
  static const bool is_physical_space = false;

public:
  using typename BaseSpace::Func;
  using typename BaseSpace::Point;
  using typename BaseSpace::Value;

  template <int order>
  using Derivative = typename BaseSpace::template Derivative<order>;

  using typename BaseSpace::Div;

  /**
   * See documentation in \ref Basis
   *
   * @see Basis
   */

//  using RefBasis = typename BaseSpace::RefBasis;

  using RefPoint = Point;

public:
  /** Type for the element accessor. */
  using ElementAccessor = BSplineElement<dim,range,rank>;

  /** Type for iterator over the elements.  */
  using ElementIterator = GridIterator<SpaceElement<dim,0,range,rank>>;


  using SpSpace = SplineSpace<dim_,range_,rank_>;
  using EndBehaviour = typename SpSpace::EndBehaviour;
  using EndBehaviourTable = typename SpSpace::EndBehaviourTable;

  using BaseSpace::ComponentContainer;

  using IndexType = TensorIndex<dim_>;
  using PropertyList = PropertiesIdContainer<IndexType>;
  using List = typename PropertyList::List;
  using ListIt = typename PropertyList::List::iterator;

public:
  /**
   * @name Creators.
   */
  ///@{
#if 0
  /**
   * Builds and returns a maximum regularity (non-const) BSpline space
   * over a (non-const) Grid
   * @p grid for the given @p degree in all directions and homogeneous
   * in all components.
   */
  static std::shared_ptr<self_t>
  create(const int degree,
         const std::shared_ptr<GridType> &grid,
         const InteriorReg interior_reg = InteriorReg::maximum,
         const bool periodic = false,
         const BasisEndBehaviour end_b = BasisEndBehaviour::interpolatory);

  /**
   * Builds and returns a maximum regularity (const) BSpline space
   * over a (const) Grid
   * @p grid for the given @p degree in all directions and homogeneous
   * in all components.
   */
  static std::shared_ptr<const self_t>
  const_create(const int degree,
               const std::shared_ptr<const GridType> &grid,
               const InteriorReg interior_reg = InteriorReg::maximum,
               const bool periodic = false,
               const BasisEndBehaviour end_b = BasisEndBehaviour::interpolatory);

  /**
   * Builds and returns a maximum regularity (non-const) BSpline space
   * over a (non-const) Grid
   * @p grid for the given @p degree[i] in the i-th direction and homogeneous
   * in all components.
   */
  static std::shared_ptr<self_t>
  create(const Degrees &degree,
         const std::shared_ptr<GridType> &grid,
         const InteriorReg interior_reg = InteriorReg::maximum,
         const Periodicity &periodic = Periodicity(false),
         const EndBehaviour &end_b = EndBehaviour(BasisEndBehaviour::interpolatory));

  /**
   * Builds and returns a maximum regularity (const) BSpline space
   * over a (const) Grid
   * @p grid for the given @p degree[i] in the i-th direction and homogeneous
   * in all components.
   */
  static std::shared_ptr<const self_t>
  const_create(const Degrees &degree,
               const std::shared_ptr<const GridType> &grid,
               const InteriorReg interior_reg = InteriorReg::maximum,
               const Periodicity &periodic = Periodicity(false),
               const EndBehaviour &end_b = EndBehaviour(BasisEndBehaviour::interpolatory));

  /**
   * Builds and returns a (non-const) BSpline space
   * over a (non-const) Grid
   * @p grid with the given multiplicity vector @p mult_vectors
   * for each component
   * and the given @p degree for each direction and for each
   * component.
   */
  static std::shared_ptr<self_t>
  create(const DegreeTable &deg,
         const std::shared_ptr<GridType> &grid,
         const MultiplicityTable &interior_mult,
         const PeriodicityTable &periodic,
         const EndBehaviourTable &end_b);

  /**
   * Builds and returns a (const) BSpline space
   * over a (const) Grid
   * @p grid with the given multiplicity vector @p mult_vectors
   * for each component
   * and the given @p degree for each direction and for each
   * component.
   */
  static std::shared_ptr<const self_t>
  const_create(const DegreeTable &deg,
               const std::shared_ptr<const GridType> &grid,
               const MultiplicityTable &interior_mult,
               const PeriodicityTable &periodic,
               const EndBehaviourTable &end_b);
#endif

  static std::shared_ptr<self_t>
  create(const std::shared_ptr<SpSpace> &spline_space,
         const EndBehaviourTable &end_b);

  static std::shared_ptr<const self_t>
  const_create(const std::shared_ptr<const SpSpace> &spline_space,
               const EndBehaviourTable &end_b);

  static std::shared_ptr<self_t>
  create(const std::shared_ptr<SpSpace> &spline_space,
         const BasisEndBehaviour &end_b = BasisEndBehaviour::interpolatory);

  static std::shared_ptr<const self_t>
  const_create(const std::shared_ptr<const SpSpace> &spline_space,
               const BasisEndBehaviour &end_b = BasisEndBehaviour::interpolatory);

  ///@}

  /**
   * Create an element (defined on this space) with a given @p index.
   */
  virtual std::unique_ptr<SpaceElement<dim_,0,range_,rank_> >
  create_element(const ListIt &index, const PropId &property) const override final;

  /**
   * Create an element (defined on this space) with a given @p index.
   */
  virtual std::unique_ptr<ReferenceElement<dim_,range_,rank_> >
  create_ref_element(const ListIt &index, const PropId &property) const override final;


  /** Destructor. */
  virtual ~BSpline() = default;



protected:
  /** @name Constructors */
  ///@{
  /**
   * Default constructor. It does nothing but it is needed for the
   * <a href="http://www.boost.org/doc/libs/release/libs/serialization/">boost::serialization</a>
   * mechanism.
   */
  BSpline() = default;
#if 0
  /**
   * Constructs a maximum regularity BSpline space
   * over a (const or non-const) Grid
   * @p grid for the given @p degree in all directions and homogeneous
   * in all components.
   */
  explicit BSpline(const int degree,
                   const SharedPtrConstnessHandler<GridType> &grid,
                   const InteriorReg interior_reg,
                   const bool periodic,
                   const BasisEndBehaviour end_b);


  /**
   * Constructs a maximum regularity BSpline space over
   * over a (const or non-const) Grid
   * @p grid for the given @p degree[i] in the i-th direction and homogeneous
   * in all components.
   */
  explicit BSpline(const Degrees &degree,
                   const SharedPtrConstnessHandler<GridType> &grid,
                   const InteriorReg interior_reg,
                   const Periodicity &periodic,
                   const EndBehaviour &end_b);


  /**
   * Constructs a BSpline space over
   * over a (const or non-const) Grid
   * @p grid with the given multiplicity vector @p mult_vectors
   * for each component
   * and the given @p degree for each direction and for each
   * component.
   */
  explicit BSpline(const DegreeTable &deg,
                   const SharedPtrConstnessHandler<GridType> &grid,
                   const MultiplicityTable &interior_mult,
                   const PeriodicityTable &periodic,
                   const EndBehaviourTable &end_b);
#endif


  explicit BSpline(const SharedPtrConstnessHandler<SpSpace> &spline_space,
                   const EndBehaviourTable &end_b);


  /**
   * Copy constructor. Not allowed to be used.
   */
  BSpline(const self_t &space) = delete;
  ///@}

  /** @name Assignment operators */
  ///@{
  /** Copy assignment. Not allowed to be used. */
  self_t &
  operator=(const self_t &space) = delete;
  ///@}

public:
  virtual std::shared_ptr<const Grid<dim_>> get_grid() const override final;


  virtual void get_element_dofs(
    const IndexType element_id,
    SafeSTLVector<Index> &dofs_global,
    SafeSTLVector<Index> &dofs_local_to_patch,
    SafeSTLVector<Index> &dofs_local_to_elem,
    const std::string &dofs_property = DofProperties::active) const override final;



  template <int k>
  using InterSpaceMap = typename BaseSpace::template InterSpaceMap<k>;



  /**
   * Construct a sub space of dimension k conforming to
   * the subspace sub element sub_elem_id and a map from the elements of
   * the sub_element grid to the corresponding element of the current
   * grid.
   */
  template<int sdim>
  std::shared_ptr<const BSpline<sdim,range_,rank_> >
  get_sub_bspline_space(const int sub_elem_id,
                        InterSpaceMap<sdim> &dof_map,
                        const std::shared_ptr<const Grid<sdim>> &sub_grid) const;

#if 0
  template<int k>
  std::shared_ptr<SubSpace<k> >
  get_sub_space(const int s_id, InterSpaceMap<k> &dof_map,
                SubGridMap<k> &elem_map) const;
#endif




  /**
   * /brief Returns the SplineSpace used to build the BSpline basis.
   */
  std::shared_ptr<const SplineSpace<dim_,range_,rank_> >
  get_spline_space() const override final;



public:


  /**
   * Returns a reference to the end behaviour table of the BSpline space.
   */
  virtual const EndBehaviourTable &get_end_behaviour_table() const override final;

  /**
   * Prints internal information about the space.
   * @note Mostly used for debugging and testing.
   */
  virtual void print_info(LogStream &out) const override final;




private:

  SharedPtrConstnessHandler<SpSpace > spline_space_;


  EndBehaviourTable end_b_;


  /** Bezier extraction operator. */
  BernsteinExtraction<dim, range, rank> operators_;


  /** If end knots are not in the repeated knot vector */
  using EndIntervalTable = typename BaseSpace::template
                           ComponentContainer<SafeSTLArray<std::pair<Real,Real>, dim>>;
  EndIntervalTable end_interval_;



  friend class BSplineElement<dim, range, rank>;
  friend class BSplineElementHandler<dim, range, rank>;

#ifdef MESH_REFINEMENT
  /**
   * Rebuild the internal state of the object after an insert_knots() function is invoked.
   *
   * @pre Before invoking this function, must be invoked the function grid_->insert_knots().
   * @note This function is connected to the Grid's signal for the refinement, and
   * it is necessary in order to avoid infinite loops in the insert_knots() function calls.
   *
   * @ingroup h_refinement
   */
  void rebuild_after_insert_knots(
    const SafeSTLArray<SafeSTLVector<Real>,dim> &knots_to_insert,
    const Grid<dim> &old_grid) override final;

public:
  virtual void refine_h(const Size n_subdivisions = 2) override final;

#endif // MESH_REFINEMENT

  /**
   * Returns the current object wrapped by a std::shared_ptr.
   *
   * @note Internally uses the shared_from_this() function.
   */
  std::shared_ptr<const self_t > get_this_basis() const;


public:
  DeclException1(ExcScalarRange, int,
                 << "Range " << arg1 << "should be 0 for a scalar valued"
                 << " space.");


  virtual bool is_bspline() const override final;

  virtual std::unique_ptr<SpaceElementHandler<dim_,0,range_,rank_>>
      create_cache_handler() const override final;


private:
#if 0
  /**
   * Lookup table for the local dof id in each element component
   */
  typename SpSpace::template ComponentContainer<SafeSTLVector<TensorIndex<dim> > >
  dofs_tensor_id_elem_table_;
#endif

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



#ifdef SERIALIZATION

#include <igatools/basis_functions/bspline.serial>

#endif // SERIALIZATION



#endif /* __BSPLINE_SPACE_H_ */

