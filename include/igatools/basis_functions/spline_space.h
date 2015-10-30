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

#ifndef SPLINE_SPACE_H_
#define SPLINE_SPACE_H_

#include <igatools/base/config.h>
#include <igatools/base/array_utils.h>
#include <igatools/functions/function.h>
#include <igatools/functions/function_element.h>
#include <igatools/utils/cartesian_product_array.h>
#include <igatools/utils/static_multi_array.h>
#include <igatools/utils/dynamic_multi_array.h>
#include <igatools/utils/shared_ptr_constness_handler.h>
//#include <igatools/basis_functions/space.h>
#include <igatools/geometry/grid.h>

IGA_NAMESPACE_OPEN




enum class BasisEndBehaviour
{
  /**
   * Interpolatory basis functions at knots boundary (i.e. open knot vector).
   */
  interpolatory,

  end_knots,

  periodic
};

// For the interior multiplicities
// maximum regularity
// minimul regularity discontinous
enum class InteriorReg
{
  maximum,
  minimun
};

/**
 * @brief Tensor product spline space
 *
 * A one dimensional polynomial spline space is determined by:
 * - a domain, the interval [a,b]
 * - the polynomial order
 * - a partition of [a,b], the knots
 * - the interior knots smoothness (given by the interior multiplicity)
 *
 * This class provides the realization of a scalar, vector or k-tensor
 * value tensor product spline space.
 *
 * @note This is independent of the basis functions one may wish to use
 * for the given space.
 *
 * @author pauletti, 2014
 * @author M. Martinelli, 2014, 2015
 *
 * @ingroup serializable
 */
template<int dim, int range = 1, int rank = 1>
class SplineSpace
{

private:
  using GridType = Grid<dim>;

public:

  using Func = Function<dim, 0, range, rank>;

public:
  template <int order>
  using Derivative = typename Func::template Derivative<order>;
  using Point = typename Func::Point;
  using Value = typename Func::Value;
  using Div   = typename Func::Div;

public:
  template<class> class ComponentContainer;
  static constexpr int n_components = ComponentContainer<Size>::n_entries;
  static const SafeSTLArray<Size, n_components> components;

public:
  using KnotCoordinates = typename GridType::KnotCoordinates;
  using BoundaryKnots = SafeSTLArray<CartesianProductArray<Real,2>, dim>;
  using Degrees  = TensorIndex<dim>;
  using Multiplicity = CartesianProductArray<Size, dim>;
  using Periodicity = SafeSTLArray<bool, dim>;
  using EndBehaviour = SafeSTLArray<BasisEndBehaviour, dim>;

  using DegreeTable = ComponentContainer<Degrees>;
  using MultiplicityTable = ComponentContainer<Multiplicity>;
  using BoundaryKnotsTable = ComponentContainer<BoundaryKnots>;
  using KnotsTable = ComponentContainer<KnotCoordinates>;
  using PeriodicityTable = ComponentContainer<Periodicity>;
  using EndBehaviourTable = ComponentContainer<EndBehaviour>;

  /**
   * Component container holding tensor size
   */
  class TensorSizeTable : public ComponentContainer<TensorSize<dim> >
  {
  public:
    using base_t = ComponentContainer<TensorSize<dim>>;


    TensorSizeTable() = default;

    TensorSizeTable(const base_t &n_basis);


    TensorSizeTable(const TensorSizeTable &in) = default;
    TensorSizeTable(TensorSizeTable &&in) = default;
    TensorSizeTable &operator=(const TensorSizeTable &in) = default;
    TensorSizeTable &operator=(TensorSizeTable &&in) = default;

    ~TensorSizeTable() = default;

    //*/



    Size get_component_size(const int comp) const;

    Size total_dimension() const;

    ComponentContainer<Size> get_offset() const;

    void print_info(LogStream &out) const;

  };




public:
  static std::shared_ptr<SplineSpace<dim,range,rank> > create(
    const DegreeTable &deg,
    const std::shared_ptr<GridType> &grid,
    const MultiplicityTable &interior_mult,
    const PeriodicityTable &periodic = PeriodicityTable(SafeSTLArray<bool,dim>(false)));

  static std::shared_ptr<const SplineSpace<dim,range,rank> > const_create(
    const DegreeTable &deg,
    const std::shared_ptr<const GridType> &grid,
    const MultiplicityTable &interior_mult,
    const PeriodicityTable &periodic = PeriodicityTable(SafeSTLArray<bool,dim>(false)));

private:
  /**
   * Construct a spline space with the knots (inferred from the (cons or non-const) @p grid
   * and @p interior_multiplicity), degree and as periodicity conditions
   */
  explicit SplineSpace(const DegreeTable &deg,
                       const SharedPtrConstnessHandler<GridType> &grid,
                       const MultiplicityTable &interior_multiplicity,
                       const PeriodicityTable &periodic =
                         PeriodicityTable(SafeSTLArray<bool,dim>(false)));


public:
  const DegreeTable &get_degree_table() const;

  const PeriodicityTable &get_periodicity() const;

  const SafeSTLArray<Index,n_components> &get_components_map() const;

  const SafeSTLVector<Index> &get_active_components_id() const;

  /** @name Getting information about the space */
  ///@{
  /**
   * Total number of basis functions. This is the dimensionality
   * of the space.
   */
  Size get_num_basis() const;

  /**
   * Total number of basis functions
   * for the comp space component.
   */
  Size get_num_basis(const int comp) const;

  /**
   *  Total number of basis functions for the comp space component
   *  and the dir direction.
   */
  Size get_num_basis(const int comp, const int dir) const;

  /**
   * Component-direction indexed table with the number of basis functions
   * in each direction and component
   */
  const TensorSizeTable &get_num_basis_table() const;
  ///@}




  template<int k>
  using SubSpace = SplineSpace<k, range, rank>;

  template<int k>
  typename SubSpace<k>::MultiplicityTable
  get_sub_space_mult(const Index s_id) const;

  template<int k>
  typename SubSpace<k>::DegreeTable
  get_sub_space_degree(const Index s_id) const;

  template<int k>
  typename SubSpace<k>::PeriodicityTable
  get_sub_space_periodicity(const Index s_id) const;


public:

  KnotsTable compute_knots_with_repetition(const EndBehaviourTable &ends,
                                           const BoundaryKnotsTable &boundary_knots = BoundaryKnotsTable()) const;

  // KnotsTable compute_knots_with_repetition(const EndBehaviourTable &ends) const;

  /**
   * For each element and for each component there is an initial
   * tensor index in the Index space from where all non-zero basis
   * function can be determined.
   */
  MultiplicityTable accumulated_interior_multiplicities() const;

  /**
   * Fill the multiplicy for the <p>regularity</p> requested
   *  of the given number of knots
   */
  static
  MultiplicityTable
  get_multiplicity_from_regularity(const InteriorReg regularity,
                                   const DegreeTable &deg,
                                   const TensorSize<dim> &n_elem);


  void print_info(LogStream &out) const;

  std::shared_ptr<Grid<dim> > get_ptr_grid();
  std::shared_ptr<const Grid<dim> > get_grid() const;


  const ComponentContainer<SafeSTLVector<TensorIndex<dim> > > &
  get_dofs_tensor_id_elem_table() const;

private:

  SharedPtrConstnessHandler<Grid<dim> > grid_;

  MultiplicityTable interior_mult_;

  DegreeTable deg_;

  /** Table with the dimensionality of the space in each component and direction */
  TensorSizeTable space_dim_;

  //EndBehaviourTable end_behaviour_;
  PeriodicityTable periodic_;


  /**
   * Lookup table for the local dofs id in each element component
   */
  ComponentContainer<SafeSTLVector<TensorIndex<dim> > > dofs_tensor_id_elem_table_;



public:

  /** Returns the multiplicity of the internal knots that defines the space. */
  const MultiplicityTable &get_interior_mult() const;

  const PeriodicityTable &get_periodic_table() const;






  /**
   *  Class to manage the component quantities with the knowledge of
   * uniform range spaces
   *
   * @ingroup serializable
   */
  template<class T>
  class ComponentContainer : public StaticMultiArray<T,range,rank>
  {
    using base_t = StaticMultiArray<T,range,rank>;
    using self_t = ComponentContainer<T>;
  public:
    /** Type of the iterator. */
    using iterator =  MultiArrayIterator<ComponentContainer<T>>;

    /** Type of the const iterator. */
    using const_iterator =  MultiArrayConstIterator<ComponentContainer<T>>;
  public:
    using base_t::n_entries;

    using  ComponentMap = SafeSTLArray<Index, n_entries>;

    ComponentContainer(const ComponentMap &comp_map =
                         sequence<n_entries>());

    ComponentContainer(const ComponentMap &comp_map, const T &val);

    ComponentContainer(bool uniform, const T &val)
      : ComponentContainer(ComponentMap(0), val)
    {}

    /**
     * Construct a homogenous range table with val value
     */
    ComponentContainer(const T &val);

    ComponentContainer(std::initializer_list<T> list);

    ComponentContainer(const self_t &in) = default;
    ComponentContainer(self_t &&in) = default;

    self_t &operator=(const self_t &in) = default;
    self_t &operator=(self_t &&in) = default;

    const_iterator
    cbegin() const
    {
      return const_iterator(*this,0);
    }

    const_iterator
    cend() const
    {
      return const_iterator(*this,IteratorState::pass_the_end);
    }


    const_iterator
    begin() const
    {
      return cbegin();
    }

    const_iterator
    end() const
    {
      return cend();
    }

    iterator
    begin()
    {
      return iterator(*this,0);
    }

    iterator
    end()
    {
      return iterator(*this,IteratorState::pass_the_end);
    }

    /**
     *  Flat index access operator (non-const version).
     */
    T &operator[](const Index i);

    /**
     *  Flat index access operator (const version).
     */
    const T &operator[](const Index i) const;

    const Index active(const Index i) const
    {
      return comp_map_[i];
    }

    const SafeSTLVector<Index> &get_active_components_id() const
    {
      return active_components_id_;
    }

    const SafeSTLVector<Index> &get_inactive_components_id() const
    {
      return inactive_components_id_;
    }

    void
    print_info(LogStream &out) const
    {
      out.begin_item("Raw components: ");
      base_t::print_info(out);
      out.end_item();

      out.begin_item("Components map: ");
      comp_map_.print_info(out);
      out.end_item();

      out.begin_item("Active components ids: ");
      active_components_id_.print_info(out);
      out.end_item();

      out.begin_item("Inactive components ids: ");
      inactive_components_id_.print_info(out);
      out.end_item();
    }

    const SafeSTLArray <Index, n_entries> &get_comp_map() const
    {
      return comp_map_;
    }

    /**
     * Equality comparison operator.
     */
    bool operator==(const self_t &table) const;


  private:
    /** For each component return the index of the active component */
    SafeSTLArray<Index, n_entries> comp_map_;

    /** list of the active components */
    SafeSTLVector<Index> active_components_id_;

    /** list of the inactive components */
    SafeSTLVector<Index> inactive_components_id_;

#ifdef SERIALIZATION
    /**
     * @name Functions needed for boost::serialization
     * @see <a href="http://uscilab.github.io/cereal/serialization_functions.html">Cereal serialization</a>
     */
    ///@{
    friend class cereal::access;

    template<class Archive>
    void
    serialize(Archive &ar)
    {
      ar &make_nvp("ComponentContainer_base_t",base_class<base_t>(this));

      ar &make_nvp("comp_map_",comp_map_);

      ar &make_nvp("active_components_id_", active_components_id_);

      ar &make_nvp("inactive_components_id_", inactive_components_id_);
    }
    ///@}
#endif // SERIALIZATION
  };


  std::shared_ptr<const SplineSpace<dim,range,rank> > spline_space_previous_refinement_;

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
    const Grid<dim> &old_grid);

  void create_connection_for_insert_knots(std::shared_ptr<SplineSpace<dim,range,rank>> space);

  void refine_h(const Size n_subdivisions);

public:
  std::shared_ptr<const SplineSpace<dim,range,rank> >
  get_spline_space_previous_refinement() const
  {
    return spline_space_previous_refinement_;
  }
#endif // MESH_REFINEMENT

protected:

  /** This function initialize the member variables from the constructor
   * arguments or after an h-refinement. */
  void init();


private:
  /**
   * Default constructor. It does nothing but it is needed for the
   * <a href="http://www.boost.org/doc/libs/release/libs/serialization/">boost::serialization</a>
   * mechanism.
   */
  SplineSpace() = default;

#ifdef SERIALIZATION
  /**
   * @name Functions needed for the serialization
   * @see <a href="http://uscilab.github.io/cereal/serialization_functions.html">Cereal serialization</a>
   */
  ///@{
  friend class cereal::access;

  template<class Archive>
  void
  serialize(Archive &ar);
  ///@}
#endif // SERIALIZATION
};


template <class T, int dim>
inline
SafeSTLVector<T>
unique_container(SafeSTLArray <T, dim> a)
{
  auto it = std::unique(a.begin(), a.end());
  return SafeSTLVector<T>(a.begin(), it);
}



template<int dim, int range, int rank>
template<class T>
SplineSpace<dim, range, rank>::
ComponentContainer<T>::
ComponentContainer(const ComponentMap &comp_map)
  :
  base_t(),
  comp_map_(comp_map),
  active_components_id_(unique_container<Index, n_entries>(comp_map)),
  inactive_components_id_(n_entries)
{
  auto all = sequence<n_entries>();
  auto it=std::set_difference(all.begin(), all.end(),
                              active_components_id_.begin(),active_components_id_.end(),
                              inactive_components_id_.begin());

  inactive_components_id_.resize(it-inactive_components_id_.begin());
}



template<int dim, int range, int rank>
template<class T>
SplineSpace<dim, range, rank>::
ComponentContainer<T>::
ComponentContainer(const ComponentMap &comp_map, const T &val)
  :
  base_t(),
  comp_map_(comp_map),
  active_components_id_(unique_container<Index, n_entries>(comp_map)),
  inactive_components_id_(n_entries)
{
  auto all = sequence<n_entries>();
  auto it=std::set_difference(all.begin(), all.end(),
                              active_components_id_.begin(),active_components_id_.end(),
                              inactive_components_id_.begin());

  inactive_components_id_.resize(it-inactive_components_id_.begin());

  for (auto i : active_components_id_)
    base_t::operator[](i) = val;
}



template<int dim, int range, int rank>
template<class T>
SplineSpace<dim, range, rank>::
ComponentContainer<T>::
ComponentContainer(std::initializer_list<T> list)
  :
  base_t(list),
  comp_map_(sequence<n_entries>()),
  active_components_id_(unique_container<Index, n_entries>(comp_map_))
{}



template<int dim, int range, int rank>
template<class T>
SplineSpace<dim, range, rank>::
ComponentContainer<T>::
ComponentContainer(const T &val)
  :
  comp_map_(0),
  active_components_id_(1,0),
  inactive_components_id_(n_entries-1)
{
  for (int i=1; i<n_entries; ++i)
    inactive_components_id_[i-1] = i;

  base_t::operator[](0) = val;
}



template<int dim, int range, int rank>
template<class T>
bool
SplineSpace<dim, range, rank>::
ComponentContainer<T>::
operator==(const self_t &table) const
{
  const bool same_comp_map = (comp_map_ == table.comp_map_);

  const bool same_active_components_id = (active_components_id_ == table.active_components_id_);

  const bool same_inactive_components_id = (inactive_components_id_ == table.inactive_components_id_);

  bool same_data = false;
  if (same_comp_map && same_active_components_id && same_inactive_components_id)
  {
    same_data = true;
    for (const auto comp : active_components_id_)
      same_data = same_data && (base_t::operator[](comp) == table[comp]);
  }

  return same_data;
}


template<int dim, int range, int rank>
template<class T>
T &
SplineSpace<dim, range, rank>::
ComponentContainer<T>::
operator[](const Index i)
{
  return base_t::operator[](comp_map_[i]);
}


template<int dim, int range, int rank>
template<class T>
const T &
SplineSpace<dim, range, rank>::
ComponentContainer<T>::
operator[](const Index i) const
{
  return base_t::operator[](comp_map_[i]);
}

IGA_NAMESPACE_CLOSE

#ifdef SERIALIZATION

#include <igatools/basis_functions/spline_space.serial>

#endif // SERIALIZATION

#endif
