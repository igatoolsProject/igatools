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
 *
 */
template<int dim_, int space_dim_>
class GridFunction :
  public std::enable_shared_from_this<GridFunction<dim_,space_dim_> >
{
private:
  using self_t = GridFunction<dim_, space_dim_>;

public:
  static const int space_dim = space_dim_;
  static const int dim = dim_;

  using GridType = Grid<dim_>;

  using ElementAccessor = GridFunctionElement<dim_, space_dim_>;
  using ElementIterator = GridIterator<ElementAccessor>;

  using ElementHandler = GridFunctionHandler<dim_, space_dim_>;

  using List = typename GridType::List;
  using ListIt = typename GridType::ListIt;

public:
  using GridPoint = typename GridType::Point;
  using Value = Values<dim, space_dim, 1>;
  template <int order>
  using Derivative = Derivatives<dim, space_dim, 1, order>;

  using Gradient = Derivative<1>;
  ///@}

  /**
   * Default constructor. It does nothing but it is needed for the
   * serialization mechanism.
   */
  GridFunction() = default;

  GridFunction(const SharedPtrConstnessHandler<GridType> &grid);


  virtual ~GridFunction() = default;




  std::shared_ptr<const GridType> get_grid() const;


  virtual std::unique_ptr<ElementHandler>
  create_cache_handler() const = 0;

  virtual std::unique_ptr<ElementAccessor>
  create_element(const ListIt &index, const PropId &prop) const;

#if 0
  std::unique_ptr<ElementAccessor>
  create_element(const ListIt &index, const PropId &prop);
#endif

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


  virtual void print_info(LogStream &out) const = 0;

#ifdef MESH_REFINEMENT
  std::shared_ptr<const self_t>
  get_grid_function_previous_refinement() const
  {
    return grid_function_previous_refinement_;
  }

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

  friend class GridFunctionElement<dim_, space_dim_>;

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
    ar &make_nvp("grid_",grid_);
  }
  ///@}
#endif // SERIALIZATION

#ifdef MESH_REFINEMENT

protected:
  std::shared_ptr<const self_t> grid_function_previous_refinement_;
#endif // MESH_REFINEMENT


public:
  virtual const SafeSTLSet<typename GridType::IndexType> &
  get_elements_with_property(const PropId &elems_property) const
  {
    return grid_->get_elements_with_property(elems_property);
  }

};



template <int,int,int> class SubGridFunctionElement;
template <int,int,int> class SubGridFunctionHandler;

template<int sdim,int dim,int space_dim>
class SubGridFunction :
  public GridFunction<sdim,space_dim>
{
private:
  using self_t = SubGridFunction<sdim,dim,space_dim>;


public:
  using base_t  = GridFunction<sdim,space_dim>;
  using SupFunc = GridFunction<dim,space_dim>;

  using GridType = Grid<sdim>;
  using SuperGrid = Grid<dim>;

  using ElementAccessor = SubGridFunctionElement<sdim,dim,space_dim>;
  using ElementIterator = GridIterator<ElementAccessor>;
  using ElementHandler = SubGridFunctionHandler<sdim,dim,space_dim>;

  using List = typename GridType::List;
  using ListIt = typename GridType::ListIt;


  //    template <int j>
  using SubGridMap = typename SuperGrid::template SubGridMap<sdim>;


public:

  SubGridFunction(const SharedPtrConstnessHandler<SupFunc> &sup_func,
                  const int s_id,
                  const SubGridMap &sub_grid_elem_map,
                  const SharedPtrConstnessHandler<GridType> &grid)
    :
    base_t(grid),
    sup_func_(sup_func),
    s_id_(s_id),
    elems_property_("boundary"),
    constant_directions_(UnitElement<dim>::template get_elem<sdim>(s_id).constant_directions),
    sub_grid_elem_map_(sub_grid_elem_map)
  {
    /*
    LogStream out;
    out.begin_item("Grid:");
    this->get_grid()->print_info(out);
    out.end_item();


    out.begin_item("Sup. Grid:");
    auto sup_grid = std::const_pointer_cast<SuperGrid>(func->get_grid());
    sup_grid->print_info(out);
    out.end_item();
    //*/
    /*
         sup_grid->add_property(elems_property_);
         for (const auto &elem_id : sub_grid_elem_map)
           sup_grid->set_property_status_elem(elems_property_,elem_id.second,true);
    //*/

//     sup_grid->print_info(out);

    /*
         auto func_elem = this->begin();
         auto func_elem_end = this->end();
         int elem_id = 0;
         for (; func_elem != func_elem_end ; ++func_elem, ++elem_id)
         {
           out.begin_item("Element " + std::to_string(elem_id));

           out << "Element ID: " << func_elem->get_index() << std::endl;

           out.end_item();
         }
    //*/


    const auto &constant_values = UnitElement<dim>::template get_elem<sdim>(s_id).constant_values;
    int i = 0;
    for (const auto &v : constant_values)
    {
      const auto dir = constant_directions_[i];
      const auto &knots = grid->get_knot_coordinates(dir);

      constant_coordinates_[i] = (v == 0)? knots.front() : knots.back();
      ++i;
    }
    LogStream out;
    for (const auto &elems_id : sub_grid_elem_map_)
    {
      id_elems_sub_grid_.insert(elems_id.first);
      id_elems_sup_grid_.insert(elems_id.second);
      out << "Sub elem ID: " << elems_id.first << "    Sup elem ID: " << elems_id.second << std::endl;
    }
    /*
         out.begin_item("Constant coordinates:");
         constant_coordinates_.print_info(out);
         out.end_item();
    //*/
//     AssertThrow(false,ExcNotImplemented());
  }

  virtual ~SubGridFunction() = default;

  virtual GridIterator<GridFunctionElement<sdim,space_dim>>
                                                         cbegin(const PropId &prop) const override
  {
    LogStream out;
    out << "cbegin" <<std::endl;

    auto elem = std::make_unique<SubGridFunctionElement<sdim,dim,space_dim>>
                (
                  std::dynamic_pointer_cast<const self_t>(this->shared_from_this()),
                  id_elems_sub_grid_.begin(),
                  prop);

    elem->print_info(out);


    //TODO: (martinelli, Nov 16,2015): the iterator is not using the property!
    return GridIterator<GridFunctionElement<sdim,space_dim>>(
             std::move(elem)
           );
  }

  virtual GridIterator<GridFunctionElement<sdim,space_dim>>
                                                         cend(const PropId &prop) const override
  {
    //TODO: (martinelli, Nov 16,2015): the iterator is not using the property!
    return GridIterator<GridFunctionElement<sdim,space_dim>>(
             std::move(std::make_unique<SubGridFunctionElement<sdim,dim,space_dim>>
                       (
                         std::dynamic_pointer_cast<const self_t>(this->shared_from_this()),
                         id_elems_sub_grid_.end(),
                         prop))
           );
  }
#if 0
  auto begin(const PropId &prop) const
  {
    return this->cbegin();
  }

  auto end(const PropId &prop) const
  {
    return this->cend();
  }
#endif

  static std::shared_ptr<const self_t>
  const_create(const std::shared_ptr<const SupFunc> &func,
               const int s_id,
               const SubGridMap &sub_grid_elem_map,
               const std::shared_ptr<const GridType> &grid)
  {
    return std::make_shared<self_t>(SharedPtrConstnessHandler<SupFunc>(func),
                                    s_id,
                                    sub_grid_elem_map,
                                    SharedPtrConstnessHandler<GridType>(grid));
  }

  static std::shared_ptr<self_t>
  create(const std::shared_ptr<const SupFunc> &func,
         const int s_id,
         const SubGridMap &sub_grid_elem_map,
         const std::shared_ptr<GridType> &grid)
  {
    return std::make_shared<self_t>(SharedPtrConstnessHandler<SupFunc>(func),
                                    s_id,
                                    sub_grid_elem_map,
                                    SharedPtrConstnessHandler<GridType>(grid));
  }

  virtual std::unique_ptr<GridFunctionHandler<sdim,space_dim> >
  create_cache_handler() const override
  {
    return std::make_unique<SubGridFunctionHandler<sdim,dim,space_dim>>(
             std::dynamic_pointer_cast<const self_t>(this->shared_from_this())
           );
  }

  virtual std::unique_ptr<GridFunctionElement<sdim,space_dim> >
  create_element(const ListIt &index, const PropId &prop) const override
  {
    using Elem = SubGridFunctionElement<sdim,dim,space_dim>;
    auto elem = std::make_unique<Elem>(
                  std::dynamic_pointer_cast<const self_t>(this->shared_from_this()),
                  index, prop);
    Assert(elem != nullptr, ExcNullPtr());

    return elem;
  }


  void rebuild_after_insert_knots(
    const SafeSTLArray<SafeSTLVector<double>, sdim> &new_knots,
    const GridType &old_grid) override final
  {
    AssertThrow(false,ExcNotImplemented());
  }

  void print_info(LogStream &out) const override final
  {
    out.begin_item("SubGridFunction<" + std::to_string(sdim) + ">");

    out.begin_item("Sup. function:");
    sup_func_->print_info(out);
    out.end_item();


    out << "Sub-element topology ID: " << s_id_ << std::endl;

    out.begin_item("Sub elements ID:");
    id_elems_sub_grid_.print_info(out);
    out.end_item();

    out.end_item();
  }


  std::shared_ptr<const SupFunc> get_sup_grid_function() const
  {
    return sup_func_.get_ptr_const_data();
  }


  const SafeSTLSet<typename Grid<sdim>::IndexType> &
  get_id_elems_sub_grid() const
  {
    return id_elems_sub_grid_;
  }

  const SafeSTLSet<typename Grid<dim>::IndexType> &
  get_id_elems_sup_grid() const
  {
    return id_elems_sup_grid_;
  }


  const typename Grid<dim>::IndexType &
  get_sup_element_id(const typename Grid<sdim>::IndexType &sub_elem_id) const
  {
    return sub_grid_elem_map_.at(sub_elem_id);
  }

  const SubGridMap &get_sub_grid_elem_map() const
  {
    return sub_grid_elem_map_;
  }

  virtual const SafeSTLSet<typename GridType::IndexType> &
  get_elements_with_property(const PropId &elems_property) const override
  {
    //TODO: (martinelli, Nov 16,2015): the property is not used!
    return id_elems_sub_grid_;
  }

private:
  SharedPtrConstnessHandler<SupFunc> sup_func_;
  const int s_id_;

  const PropId elems_property_;

  const SafeSTLArray<int,dim-sdim> constant_directions_;

  SafeSTLArray<Real,dim-sdim> constant_coordinates_;

  const SubGridMap sub_grid_elem_map_;

  SafeSTLSet<typename Grid<sdim>::IndexType> id_elems_sub_grid_;
  SafeSTLSet<typename Grid< dim>::IndexType> id_elems_sup_grid_;
  //  typename SupFunc::ElementIterator sup_elem_;

};


IGA_NAMESPACE_CLOSE

#endif

