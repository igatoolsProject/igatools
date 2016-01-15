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

#ifndef __GRID_FUNCTION_H_
#define __GRID_FUNCTION_H_

#include <igatools/base/config.h>
#include <igatools/geometry/grid.h>
#include <igatools/geometry/grid_handler.h>
#include <igatools/utils/shared_ptr_constness_handler.h>

IGA_NAMESPACE_OPEN

template <int, int> class GridFunctionElement;
template <int, int> class GridFunctionHandler;

/**
 * @brief Base class for function mapping the reference (or parametric) domain
 * \f$ \hat{\Omega} \subset \mathbb{R}^{dim} \f$ (represented by Grid) to \f$\mathbb{R}^{range}\f$.
 *
 * @ingroup serializable
 */
template<int dim_, int range_>
class GridFunction :
  public std::enable_shared_from_this<GridFunction<dim_,range_> >
{
private:
  using self_t = GridFunction<dim_, range_>;

public:
  static const int dim = dim_;
  static const int range = range_;

  using GridType = Grid<dim_>;

  using ElementAccessor = GridFunctionElement<dim_,range_>;
  using ElementIterator = GridIterator<ElementAccessor>;

  using Handler = GridFunctionHandler<dim_,range_>;

  using List = typename GridType::List;
  using ListIt = typename GridType::ListIt;

public:
  using GridPoint = typename GridType::Point;
  using Value = Values<dim,range, 1>;
  template <int order>
  using Derivative = Derivatives<dim,range,1,order>;

  using Gradient = Derivative<1>;
  ///@}

  /**
   * Default constructor. It sets the unique value for the object ID.
   *
   * @note  This constructor is needed for the serialization mechanism.
   */
  GridFunction();

  GridFunction(const SharedPtrConstnessHandler<GridType> &grid);


  virtual ~GridFunction() = default;




  std::shared_ptr<const GridType> get_grid() const;


  virtual std::unique_ptr<Handler>
  create_cache_handler() const = 0;

  virtual std::unique_ptr<ElementAccessor>
  create_element_begin(const PropId &prop) const;

  virtual std::unique_ptr<ElementAccessor>
  create_element_end(const PropId &prop) const;


  ///@name Iterating of grid elements
  ///@{
#if 0
  /**
   * This function returns a element iterator to the first element of the patch.
   */
  ElementIterator begin(const PropId &prop = ElementProperties::active);

  /**
   * This function returns a element iterator to one-pass the end of patch.
   */
  ElementIterator end(const PropId &prop = ElementProperties::active);
#endif

  /**
   * This function returns a element (const) iterator to the first element of the patch.
   */
  ElementIterator begin(const PropId &prop = ElementProperties::active) const;

  /**
   * This function returns a element (const) iterator to one-pass the end of patch.
   */
  ElementIterator end(const PropId &prop = ElementProperties::active) const;

  /**
   * This function returns a element (const) iterator to the first element of the patch.
   */
  virtual ElementIterator cbegin(const PropId &prop = ElementProperties::active) const;

  /**
   * This function returns a element (const) iterator to one-pass the end of patch.
   */
  virtual ElementIterator cend(const PropId &prop = ElementProperties::active) const;
  ///@}

  /**
   * Get the name associated to the object instance.
   */
  const std::string &get_name() const;

  /**
   * Set the name associated to the object instance.
   */
  void set_name(const std::string &name);

  /**
   * Returns the unique identifier associated to each object instance.
   */
  Index get_object_id() const;


  virtual void print_info(LogStream &out) const = 0;

#ifdef MESH_REFINEMENT
  std::shared_ptr<const self_t>
  get_grid_function_previous_refinement() const;

private:
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

public:

  /**
   *  Connect a slot (i.e. a function pointer) to the refinement signals
   *  which will be
   *  emitted whenever a insert_knots() function is called by the underlying
   *  a Grid member.
   */
  boost::signals2::connection
  connect_insert_knots(const typename Grid<dim_>::SignalInsertKnotsSlot &subscriber);

  void create_connection_for_insert_knots(const std::shared_ptr<self_t> &grid_function);
#endif // MESH_REFINEMENT

private:
  SharedPtrConstnessHandler<Grid<dim_>> grid_;

protected:
  /// Name.
  std::string name_;

private:
  /**
   * Unique identifier associated to each object instance.
   */
  const Index object_id_;

  friend class GridFunctionElement<dim_, range_>;

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

#ifdef MESH_REFINEMENT

protected:
  std::shared_ptr<const self_t> grid_function_previous_refinement_;
#endif // MESH_REFINEMENT


public:
  virtual const SafeSTLVector<typename GridType::IndexType> &
  get_elements_with_property(const PropId &elems_property) const
  {
    return grid_->get_elements_with_property(elems_property);
  }

};




IGA_NAMESPACE_CLOSE

#endif

