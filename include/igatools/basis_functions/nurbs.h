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

#ifndef NURBS_SPACE_H_
#define NURBS_SPACE_H_

#include <igatools/base/config.h>

#include <igatools/basis_functions/bspline.h>
#include <igatools/basis_functions/bspline_element.h>
#include <igatools/functions/ig_grid_function.h>

#ifdef USE_NURBS

IGA_NAMESPACE_OPEN

template <int, int, int> class NURBSElement;
template <int, int, int> class NURBSElementHandler;

/**
 * Multivariate (tensor product) scalar, vector or k-tensor
 * valued NURBS space.
 *
 * @ingroup containers
 * @ingroup serializable
 */
template <int dim_, int range_ = 1, int rank_ = 1>
class NURBS :
  public ReferenceSpaceBasis<dim_,range_,rank_>
{
private:
  using BaseSpace = ReferenceSpaceBasis<dim_,range_,rank_>;
  using self_t = NURBS<dim_, range_, rank_>;

public:
  using BSpBasis = BSpline<dim_, range_, rank_>;


  /** see documentation in \ref Basis */

  using GridType = Grid<dim_>;
  static const int dim       = dim_;
  static const int codim     = 0;
  static const int space_dim = dim_;
  static const int range     = range_;
  static const int rank      = rank_;
  static const bool is_physical_space = false;

  static const auto n_components = SplineSpace<dim_, range_, rank_>::n_components;


  /**
   * See documentation in \ref Basis
   *
   * @see Basis
   */

//  using RefBasis = typename BaseSpace::RefBasis;


  using IndexType = typename GridType::IndexType;
  using PropertyList = PropertiesIdContainer<IndexType>;
  using List = typename PropertyList::List;
  using ListIt = typename PropertyList::List::iterator;


public:
  template <int order>
  using Derivative = typename BSpBasis::template Derivative<order>;
  using Point = typename BSpBasis::Point;
  using Value = typename BSpBasis::Value;
  using Div   = typename BSpBasis::Div;

  using RefPoint = typename BSpBasis::RefPoint;

public:

  /** Type for the element accessor. */
  using ElementAccessor = NURBSElement<dim, range, rank> ;

  /** Type for iterator over the elements.  */
  using ElementIterator = GridIterator<ReferenceElement<dim,range,rank> >;

  using ElementHandler = NURBSElementHandler<dim_, range_, rank_>;




  template <int k>
  using InterSpaceMap = typename BaseSpace::template InterSpaceMap<k>;


  /**
   * Construct a sub space of dimension <tt>sdim</tt> conforming to
   * the subspace sub element sub_elem_id and a map from the elements of
   * the sub_element grid to the corresponding element of the current
   * grid.
   */
  template<int sdim>
  std::shared_ptr<const NURBS<sdim,range_,rank_> >
  get_sub_nurbs_space(const int sub_elem_id,
                      InterSpaceMap<sdim> &dof_map,
                      const std::shared_ptr<const Grid<sdim>> &sub_grid) const;
#if 0
  template<int k>
  std::shared_ptr<SubSpace<k> >
  get_sub_space(const int s_id, InterSpaceMap<k> &dof_map,
                SubGridMap<k> &elem_map) const;
#endif

public:
//    /** Container indexed by the components of the space */
  template< class T>
  using ComponentContainer = typename BSpBasis::template ComponentContainer<T>;

  using EndBehaviour = typename BSpBasis::EndBehaviour;
  using EndBehaviourTable = typename BSpBasis::EndBehaviourTable;


  using WeightSpace = BSpline<dim_,1,1>;
  using WeightFunction = IgGridFunction<dim_,1>;
  using WeightFunctionPtr = std::shared_ptr<WeightFunction>;
  using Weights = DynamicMultiArray<Real,dim>;
  using WeightsTable = ComponentContainer<Weights>;

public:
  /**
   * Returns a shared_ptr wrapping a (non-const) NURBS from a
   * (non-const) BSpline and a scalar weight function.
   */
  static std::shared_ptr<self_t>
  create(const std::shared_ptr<BSpBasis> &bs_space,
         const std::shared_ptr<WeightFunction> &weight_func);

  /**
   * Returns a shared_ptr wrapping a (const) NURBS from a
   * (const) BSpline and a scalar weight function.
   */
  static std::shared_ptr<const self_t>
  const_create(const std::shared_ptr<const BSpBasis> &bs_space,
               const std::shared_ptr<const WeightFunction> &weight_func);

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

  /** Destructor */
  virtual ~NURBS() = default;

protected:
  /** @name Constructor */
  ///@{
  /**
   * Default constructor. It does nothing but it is needed for the
   * serialization mechanism.
   */
  NURBS() = default;

  /**
   * Construct a NURBS from a BSpline and a scalar weight function.
   */
  explicit NURBS(const SharedPtrConstnessHandler<BSpBasis> &bsp_basis,
                 const SharedPtrConstnessHandler<WeightFunction> &weight_func);

  /**
   * Copy constructor. Not allowed to be used.
   */
  NURBS(const self_t &space) = delete;

  ///@}

public:
  /** @name Getting information about the space */
  ///@{

  virtual bool is_bspline() const override final;


  virtual void get_element_dofs(
    const IndexType &element_id,
    SafeSTLVector<Index> &dofs_global,
    SafeSTLVector<Index> &dofs_local_to_patch,
    SafeSTLVector<Index> &dofs_local_to_elem,
    const std::string &dofs_property = DofProperties::active) const override final;

  ///@}

  const std::shared_ptr<const BSpBasis> get_bspline_basis() const;

  virtual std::shared_ptr<const Grid<dim_>> get_grid() const override final;


  /**
   * Get the weights function of the NURBS.
   */
  std::shared_ptr<const WeightFunction> get_weight_func() const;


  /**
   * Returns a const reference to the end behaviour table of the BSpline space.
   */
  virtual const EndBehaviourTable &get_end_behaviour_table() const override final;


  /**
   * Prints internal information about the space.
   * @note Mostly used for debugging and testing.
   */
  virtual void print_info(LogStream &out) const override final;



  /**
   * /brief Returns the SplineSpace used to build the BSpline basis.
   */
  std::shared_ptr<const SplineSpace<dim_,range_,rank_>>
                                                     get_spline_space() const override final;


private:
  /**
   * B-spline space
   */
  SharedPtrConstnessHandler<BSpBasis> bsp_basis_;

  /**
   * Weight function.
   */
  SharedPtrConstnessHandler<WeightFunction> weight_func_;



  friend class NURBSElement<dim, range, rank>;
  friend class NURBSElementHandler<dim, range, rank>;


  /**
   * Returns the current object wrapped by a std::shared_ptr.
   *
   * @note Internally uses the shared_from_this() function.
   */
  std::shared_ptr<const self_t > get_this_basis() const;


public:
  virtual std::unique_ptr<SpaceElementHandler<dim_,0,range_,rank_> >
  create_cache_handler() const override final;


private:

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
    const SafeSTLArray<SafeSTLVector<Real>,dim_> &knots_to_insert,
    const Grid<dim_> &old_grid) override final;

public:
  virtual void refine_h(const Size n_subdivisions = 2) override final;

//  void create_connection_for_insert_knots(const std::shared_ptr<self_t> &space);
#endif // MESH_REFINEMENT



#ifdef SERIALIZATION
private:
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
    const std::string base_name = "ReferenceSpaceBasis_" +
                                  to_string(dim_) + "_" +
                                  to_string(0) + "_" +
                                  to_string(range_) + "_" +
                                  to_string(rank_);

    ar &make_nvp(base_name,base_class<BaseSpace>(this));
    ar &make_nvp("bsp_basis_",bsp_basis_);

    ar &make_nvp("weight_func_",weight_func_);
  }
  ///@}
#endif // SERIALIZATION
};

IGA_NAMESPACE_CLOSE




#ifdef SERIALIZATION

#include <igatools/basis_functions/nurbs.serial>

#endif // SERIALIZATION


#endif /* #ifdef USE_NURBS */


#endif /* NURBS_SPACE_H_ */


