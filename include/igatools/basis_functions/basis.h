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



IGA_NAMESPACE_OPEN

template <int,int> class Domain;

template <int,int,int,int> class Function;
//template <int> class NonConstGridElement;
template <int,int,int,int> class SpaceElement;
template <int,int,int,int> class SpaceElementHandler;


template <int,int,int> class SplineSpace;


template <int,int,int> class DofDistribution;


/**
 * @brief This is an auxiliary class used represent the "concept" of isogeometric basis function
 * in which its space is defined over <tt>dim</tt>-dimensional Grid.
 *
 * It is used as base class of ReferenceSpaceBasis and PhysicalSpaceBasis.
 *
 * @author martinelli, 2015.
 *
 * @ingroup containers
 * @ingroup serializable
 */
template <int dim_,int codim_,int range_,int rank_>
class Basis
  :
  public std::enable_shared_from_this<Basis<dim_,codim_,range_,rank_> >
{
private:
  using self_t = Basis<dim_,codim_,range_,rank_>;


public:
  static const int space_dim = dim_+codim_;

  /** Types for the input/output evaluation arguments */
  ///@{
  /**
   * Type for the input argument of the function.
   */
  using Point = Points<space_dim>;

  /**
   * Type for the return of the function.
   */
  using Value = Values<space_dim, range_, rank_>;

  /**
   * Type for the derivative of the function.
   */
  template <int order>
  using Derivative = Derivatives<space_dim, range_, rank_, order>;

  /**
   * Type for the gradient of the function.
   */
  using Gradient = Derivative<1>;

  /**
   * Type for the hessian of the function.
   */
  using Hessian = Derivative<2>;

  /**
   * Type for the divergence of function.
   */
  using Div = Values<space_dim, space_dim, rank_-1>;
  ///@}


  /** @name Constructor and destructor. */
  ///@{
protected:
  /**
   * Default constructor. It does nothing, apart assigning a unique id to the object.
   */
  Basis();



  /** Copy constructor. */
  Basis(const self_t &) = delete;

  /** Move constructor. */
  Basis(self_t &&) = default;

public:
  /** Destructor. */
  virtual ~Basis() = default;
  ///@}

  /** @name Assignment operator. */
  ///@{
  /** Copy assignment operator. Not allowed to be used. */
  self_t &operator=(const self_t &) = delete;

  /** Move assignment operator. Not allowed to be used. */
  self_t &operator=(self_t &&) = delete;
  ///@}


public:

  using GridType = Grid<dim_>;
  using IndexType = typename GridType::IndexType;
  using PropertyList = PropertiesIdContainer<IndexType>;
  using List = typename PropertyList::List;
  using ListIt = typename PropertyList::List::iterator;



  /**
   * \brief Returns the unique identifier associated to each object instance.
   */
  Index get_object_id() const;


  /**
   * \brief Get the name associated to the object instance.
   */
  const std::string &get_name() const;

  /**
   * \brief Set the name associated to the object instance.
   */
  void set_name(const std::string &name);


  /**
   * \brief Returns the Grid upon which the space is built.
   */
  virtual std::shared_ptr<const Grid<dim_>> get_grid() const = 0;


  using topology_variant = TopologyVariants<dim_>;

  /**
   * \name Functions for getting information about the dofs
   */
  ///@{
  std::shared_ptr<const DofDistribution<dim_,range_,rank_> >
  get_dof_distribution() const;




  virtual void get_element_dofs(
    const IndexType &element_id,
    SafeSTLVector<Index> &dofs_global,
    SafeSTLVector<Index> &dofs_local_to_patch,
    SafeSTLVector<Index> &dofs_local_to_elem,
    const std::string &dofs_property = DofProperties::active) const = 0;
  ///@}

  /** @name Functions for retrieving information about the number of basis function. */
  ///@{
  Size get_num_basis() const;


//  const std::set<Index> &get_global_dofs(const std::string &dof_prop = DofProperties::active) const ;

  ///@}


  /**
   * Create and element (defined on this space) with a given flat_index
   */
  virtual std::unique_ptr<SpaceElement<dim_,codim_,range_,rank_>>
      create_element(const ListIt &index, const PropId &property) const = 0;



  virtual std::unique_ptr<SpaceElementHandler<dim_,codim_,range_,rank_> >
  create_cache_handler() const = 0;


  virtual void print_info(LogStream &out) const = 0;


  virtual
  std::shared_ptr<const SplineSpace<dim_,range_,rank_> >
  get_spline_space() const = 0;


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
  virtual void refine_h(const Size n_subdivisions = 2) = 0;


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
    const SafeSTLArray<SafeSTLVector<Real>,dim_> &knots_to_insert,
    const Grid<dim_> &old_grid) = 0;

//  void create_connection_for_insert_knots(const std::shared_ptr<self_t> &space);


#if 0
  virtual std::shared_ptr<const self_t> get_space_previous_refinement() const = 0;
#endif

#endif


public:
  static const int dim = dim_;
  static const int codim = codim_;
  static const int range = range_;
  static const int rank = rank_;

private:

  /**
   * \brief Unique identifier associated to each object instance.
   */
  Index object_id_ = 0;

  /**
   * \brief Name associated to the object instance.
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
