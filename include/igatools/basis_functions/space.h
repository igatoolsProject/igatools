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

#ifndef __SPACE_H_
#define __SPACE_H_

#include <igatools/base/config.h>
#include <igatools/utils/shared_ptr_constness_handler.h>
#include <igatools/geometry/grid.h>
//#include <igatools/base/function.h>

//#include <igatools/basis_functions/space_element.h>



IGA_NAMESPACE_OPEN

template <int,int> class Domain;

template <int,int,int,int> class Function;
//template <int> class NonConstGridElement;
template <int,int,int,int> class SpaceElement;
template <int,int,int,int> class SpaceElementHandler;




template <int,int,int> class DofDistribution;


/**
 * * @brief This is an auxiliary class used represent the "concept" of isogeometric function space, defined
 * over <tt>dim</tt>-dimensional Grid.

 * It is used as base class of ReferenceSpace and PhysicalSpace.
 *
 * @author martinelli, 2015.
 *
 * @ingroup serializable
 */
template <int dim_,int codim_,int range_,int rank_>
class Space
  :
  public std::enable_shared_from_this<Space<dim_,codim_,range_,rank_> >
{
private:
  using self_t = Space<dim_,codim_,range_,rank_>;



  /** @name Constructor and destructor. */
  ///@{
protected:
  /**
   * Default constructor. It does nothing, apart assigning a unique id to the object.
   */
  Space();



  /** Copy constructor. */
  Space(const self_t &) = delete;

  /** Move constructor. */
  Space(self_t &&) = default;

public:
  /** Destructor. */
  virtual ~Space() = default;
  ///@}

  /** @name Assignment operator. */
  ///@{
  /** Copy assignment operator. Not allowed to be used. */
  self_t &operator=(const self_t &) = delete;

  /** Move assignment operator. Not allowed to be used. */
  self_t &operator=(self_t &&) = delete;
  ///@}


public:

  using IndexType = TensorIndex<dim_>;
  using PropertyList = PropertiesIdContainer<IndexType>;
  using List = typename PropertyList::List;
  using ListIt = typename PropertyList::List::iterator;


  static const int dim = dim_;
  static const int codim = codim_;
  static const int range = range_;
  static const int rank = rank_;

  /**
   * Returns the unique identifier associated to each object instance.
   */
  Index get_object_id() const;


  virtual std::shared_ptr<Grid<dim_>> get_ptr_grid() const = 0;

  virtual std::shared_ptr<const Grid<dim_>> get_ptr_const_grid() const = 0;

  /**
   * Get the name associated to the object instance.
   */
  const std::string &get_name() const;

  /**
   * Set the name associated to the object instance.
   */
  void set_name(const std::string &name);


  virtual std::shared_ptr<const DofDistribution<dim_,range_,rank_> >
  get_ptr_const_dof_distribution() const = 0;


  virtual std::shared_ptr<DofDistribution<dim_,range_,rank_> >
  get_ptr_dof_distribution() = 0;


  /**
   * Returns the dofs that are on the interior of the <tt>dim</tt>-dimensional hypercube
   * (i.e. the dofs that are not on the boundary).
   */
  std::set<Index> get_interior_dofs() const;

  using topology_variant = TopologyVariants<dim_>;

  std::set<Index> get_boundary_dofs(const int s_id, const topology_variant &topology) const;

  template<int k>
  std::set<Index> get_boundary_dofs(const int s_id) const
  {
    return this->get_boundary_dofs(s_id,Topology<k>());
  }

  /** @name Functions for retrieving information about the number of basis function. */
  ///@{
  Size get_num_basis() const;

  Size get_num_basis(const int comp) const;

  Size get_num_basis(const int comp, const int dir) const;

  Size get_elem_num_basis() const;

  /**
   * This function returns the global dof id corresponding to the basis function
   * with tensor index <p>tensor_index</p> on the @p comp component of the space.
   */
  Index
  get_global_dof_id(const TensorIndex<dim> &tensor_index,
                    const Index comp) const;
  ///@}



  virtual void get_element_dofs(
    const IndexType element_id,
    SafeSTLVector<Index> &dofs_global,
    SafeSTLVector<Index> &dofs_local_to_patch,
    SafeSTLVector<Index> &dofs_local_to_elem,
    const std::string &dofs_property = DofProperties::active) const  = 0;
//  {
//    Assert(false,ExcMessage("This function must be pure abstract!"));
//  }


  /**
   * Return the maximum value of the polynomial degree, for each component, for each direction;
   */
  virtual int get_max_degree() const = 0;


  /**
   * Create and element (defined on this space) with a given flat_index
   */
  virtual std::unique_ptr<SpaceElement<dim_,codim_,range_,rank_>>
      create_element(const ListIt &index, const PropId &property) const = 0;



  virtual std::unique_ptr<SpaceElementHandler<dim_,codim_,range_,rank_> >
  create_cache_handler() const = 0;


  virtual void print_info(LogStream &out) const = 0;



  using ElementAccessor = SpaceElement<dim_,codim_,range_,rank_>;
  using ElementIterator = GridIterator<ElementAccessor>;

  /** @name Functions involving the element iterator */
  ///@{
  /**
   * Returns an element iterator pointing to the first element of the patch
   * with the property @p element_property.
   */
  ElementIterator begin(const PropId &element_property = ElementProperties::active) const ;


  /**
   * Returns a element iterator pointing to one-pass the end of patch
   * with the property @p element_property.
   */
  ElementIterator end(const PropId &element_property = ElementProperties::active) const;

  /**
   * Returns an element iterator pointing to the first element of the patch
   * with the property @p element_property.
   */
  ElementIterator cbegin(const PropId &element_property = ElementProperties::active) const ;


  /**
   * Returns a element iterator pointing to one-pass the end of patch
   * with the property @p element_property.
   */
  ElementIterator cend(const PropId &element_property = ElementProperties::active) const;
///@}


#ifdef MESH_REFINEMENT

  /**
   * Perform the h-refinement of the space in all the directions.
   *
   * Each interval in the unrefined grid is uniformly divided in @p n_subdivisions
   * sub-intervals.
   *
   * @ingroup h_refinement
   */
  void refine_h(const Size n_subdivisions = 2);


  /**
   * Rebuild the internal state of the object after an insert_knots() function is invoked.
   *
   * @pre Before invoking this function, must be invoked the function grid_->insert_knots().
   * @note This function is connected to the Grid's signal for the refinement, and
   * it is necessary in order to avoid infinite loops in the insert_knots() function calls.
   *
   * @ingroup h_refinement
   */
  virtual void rebuild_after_insert_knots(
    const SafeSTLArray<SafeSTLVector<Real>,dim> &knots_to_insert,
    const Grid<dim> &old_grid) = 0;

  void create_connection_for_insert_knots(const std::shared_ptr<self_t> &space);


#if 0
  virtual std::shared_ptr<const self_t> get_space_previous_refinement() const = 0;
#endif

#endif


private:

  /**
   * Unique identifier associated to each object instance.
   */
  Index object_id_ = 0;

  /**
   * Name associated to the object instance.
   */
  std::string name_;



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


#endif // __SPACE_H_
