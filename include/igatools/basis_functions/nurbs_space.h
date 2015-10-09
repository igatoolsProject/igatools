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

#include <igatools/basis_functions/bspline_space.h>
#include <igatools/basis_functions/bspline_element.h>
#include <igatools/functions/ig_grid_function.h>

#ifdef NURBS

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
class NURBSSpace :
  public ReferenceSpace<dim_,range_,rank_>
{
private:
  using BaseSpace = ReferenceSpace<dim_,range_,rank_>;
  using self_t = NURBSSpace<dim_, range_, rank_>;

public:
  using BSpSpace = BSplineSpace<dim_, range_, rank_>;


  /** see documentation in \ref Space */

  using GridType = Grid<dim_>;
  static const int dim       = dim_;
  static const int codim     = 0;
  static const int space_dim = dim_;
  static const int range     = range_;
  static const int rank      = rank_;
  static const bool is_physical_space = false;

  static const auto n_components = SplineSpace<dim_, range_, rank_>::n_components;


  /**
   * See documentation in \ref Space
   *
   * @see Space
   */

  using RefSpace = typename BaseSpace::RefSpace;


  using IndexType = TensorIndex<dim_>;
  using PropertyList = PropertiesIdContainer<IndexType>;
  using List = typename PropertyList::List;
  using ListIt = typename PropertyList::List::iterator;


public:
  using Func = typename BSpSpace::Func;
  template <int order>
  using Derivative = typename BSpSpace::template Derivative<order>;
  using Point = typename BSpSpace::Point;
  using Value = typename BSpSpace::Value;
  using Div   = typename BSpSpace::Div;

  using RefPoint = typename BSpSpace::RefPoint;

public:

  /** Type for the element accessor. */
  using ElementAccessor = NURBSElement<dim, range, rank> ;

  /** Type for iterator over the elements.  */
  using ElementIterator = GridIterator<ReferenceElement<dim,range,rank> >;

  using ElementHandler = NURBSElementHandler<dim_, range_, rank_>;



  template <int sdim>
  using SubGridMap = typename GridType::template SubGridMap<sdim>;

//    using InterGridMap = std::map<Index,Index>;

  template <int k>
  using InterSpaceMap = SafeSTLVector<Index>;

  template <int k>
  using SubRefSpace = NURBSSpace<k, range, rank>;

  template <int k>
  using SubSpace = PhysicalSpace<k,range,rank, dim-k>;

  /**
   * Construct a sub space of dimension k conforming to
   * the subspace sub element sub_elem_id and a map from the elements of
   * the sub_element grid to the corresponding element of the current
   * grid.
   */
  template<int k>
  std::shared_ptr<SubRefSpace<k> >
  get_ref_sub_space(const int sub_elem_id,
                    InterSpaceMap<k> &dof_map,
                    std::shared_ptr<Grid<k>> sub_grid = nullptr) const;

  template<int k>
  std::shared_ptr<SubSpace<k> >
  get_sub_space(const int s_id, InterSpaceMap<k> &dof_map,
                std::shared_ptr<Grid<k>> sub_grid,
                SubGridMap<k> &elem_map) const;

public:
//    /** Container indexed by the components of the space */
  template< class T>
  using ComponentContainer = typename BSpSpace::template ComponentContainer<T>;

  using Degrees = typename BSpSpace::Degrees;
  using Multiplicity = typename BSpSpace::Multiplicity;
  using EndBehaviour = typename BSpSpace::EndBehaviour;
  using Periodicity = typename BSpSpace::Periodicity;

  using KnotsTable = typename BSpSpace::KnotsTable;
  using DegreeTable = typename BSpSpace::DegreeTable;
  using MultiplicityTable = typename BSpSpace::MultiplicityTable;
  using TensorSizeTable = typename BSpSpace::TensorSizeTable;
  using PeriodicityTable = typename BSpSpace::PeriodicityTable;
  using EndBehaviourTable = typename BSpSpace::EndBehaviourTable;


  using WeightSpace = BSplineSpace<dim_,1,1>;
  using WeightFunction = GridFunction<dim_,1>;
  using WeightFunctionPtr = std::shared_ptr<WeightFunction>;
  using Weights = DynamicMultiArray<Real,dim>;
  using WeightsTable = ComponentContainer<Weights>;

public:
  /**
   * Returns a shared_ptr wrapping a (non-const) NURBSSpace from a
   * (non-const) BSplineSpace and a scalar weight function.
   */
  static std::shared_ptr<self_t>
  create(const std::shared_ptr<BSpSpace> &bs_space,
         const std::shared_ptr<WeightFunction> &weight_func);

  /**
   * Returns a shared_ptr wrapping a (const) NURBSSpace from a
   * (const) BSplineSpace and a scalar weight function.
   */
  static std::shared_ptr<const self_t>
  const_create(const std::shared_ptr<const BSpSpace> &bs_space,
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
  virtual ~NURBSSpace() = default;

protected:
  /** @name Constructor */
  ///@{
  /**
   * Default constructor. It does nothing but it is needed for the
   * serialization mechanism.
   */
  NURBSSpace() = default;

  /**
   * Construct a NURBSSpace from a BSplineSpace and a scalar weight function.
   */
  explicit NURBSSpace(const SharedPtrConstnessHandler<BSpSpace> &bsp_space,
                      const SharedPtrConstnessHandler<WeightFunction> &weight_func);

  /**
   * Copy constructor. Not allowed to be used.
   */
  NURBSSpace(const self_t &space) = delete;

  ///@}

public:
  /** @name Getting information about the space */
  ///@{

  virtual bool is_bspline() const override final;

  virtual const DegreeTable &get_degree_table() const override final;

  virtual void get_element_dofs(
    const IndexType element_id,
    SafeSTLVector<Index> &dofs_global,
    SafeSTLVector<Index> &dofs_local_to_patch,
    SafeSTLVector<Index> &dofs_local_to_elem,
    const std::string &dofs_property = DofProperties::active) const override final;

  ///@}

  const std::shared_ptr<const BSpSpace> get_spline_space() const;

  virtual std::shared_ptr<Grid<dim_>> get_ptr_grid() const override final
  {
    return bsp_space_.get_ptr_data()->get_ptr_grid();
  }

  virtual std::shared_ptr<const Grid<dim_>> get_ptr_const_grid() const override final
  {
    return bsp_space_.get_ptr_const_data()->get_ptr_const_grid();
  }


  /**
   * Get the weights function of the NURBSSpace.
   */
  std::shared_ptr<const WeightFunction> get_weight_func() const;

  const PeriodicityTable &get_periodicity() const override final;


  /**
   * Returns a const reference to the end behaviour table of the BSpline space.
   */
  virtual const EndBehaviourTable &get_end_behaviour_table() const override final;


  /**
   * Prints internal information about the space.
   * @note Mostly used for debugging and testing.
   */
  virtual void print_info(LogStream &out) const override final;

  std::shared_ptr<const DofDistribution<dim, range, rank> >
  get_ptr_const_dof_distribution() const override final;

  std::shared_ptr<DofDistribution<dim, range, rank> >
  get_ptr_dof_distribution() override final;



private:
  /**
   * B-spline space
   */
  SharedPtrConstnessHandler<BSpSpace> bsp_space_;

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
  std::shared_ptr<const self_t > get_this_space() const;


public:
  virtual std::unique_ptr<SpaceElementHandler<dim_,0,range_,rank_>>
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
    const Grid<dim_> &old_grid);

  void create_connection_for_insert_knots(std::shared_ptr<self_t> space);
#endif // MESH_REFINEMENT



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
    const std::string base_name = "ReferenceSpace_" +
                                  to_string(dim_) + "_" +
                                  to_string(0) + "_" +
                                  to_string(range_) + "_" +
                                  to_string(rank_) + "_hgrad";

    ar &make_nvp(base_name,base_class<BaseSpace>(this));
    ar &make_nvp("bsp_space_",bsp_space_);

    ar &make_nvp("weight_func_",weight_func_);
  }
  ///@}
#endif // SERIALIZATION
};

IGA_NAMESPACE_CLOSE




#ifdef SERIALIZATION

#include <igatools/basis_functions/nurbs_space.serial>

#endif // SERIALIZATION


#endif /* #ifdef NURBS */


#endif /* NURBS_SPACE_H_ */


