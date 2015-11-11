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

#ifndef __FORMULA_GRID_FUNCTION_H_
#define __FORMULA_GRID_FUNCTION_H_

#include <igatools/base/value_types.h>
#include <igatools/geometry/grid_function.h>
#include <igatools/geometry/grid_function_element.h>
//#include <igatools/geometry/sub_grid_function.h>

IGA_NAMESPACE_OPEN

template <int, int> class FormulaGridFunctionHandler;

/**
 *
 */
template<int dim, int space_dim>
class FormulaGridFunction :
  public GridFunction<dim, space_dim>
{
private:
  using parent_t =  GridFunction<dim, space_dim>;
  using self_t = FormulaGridFunction<dim, space_dim>;
protected:
  using typename parent_t::GridType;
  using ElementHandler = FormulaGridFunctionHandler<dim, space_dim>;
public:
  using typename parent_t::Value;
  using typename parent_t::GridPoint;

  template <int order>
  using Derivative = typename parent_t::template Derivative<order>;

  FormulaGridFunction(const SharedPtrConstnessHandler<GridType> &grid);

  virtual ~FormulaGridFunction() = default;

  std::unique_ptr<typename parent_t::ElementHandler>
  create_cache_handler() const;


  template <int sdim>
  using SubGridElemMap = typename Grid<dim>::template SubGridMap<sdim>;

  template <int sdim>
  std::shared_ptr<const FormulaGridFunction<sdim,space_dim> >
  get_sub_function(const int s_id,
                   const SubGridElemMap<sdim> &sub_grid_elem_map,
                   const std::shared_ptr<const Grid<sdim>> &sub_grid) const
  {
    static_assert(sdim == 0 || (sdim > 0 && sdim < dim),
                  "The dimensionality of the sub_grid is not valid.");

    auto sub_func = SubFormulaGridFunction<sdim>::const_create(
                      this->shared_from_this(),s_id,sub_grid_elem_map,sub_grid);
    AssertThrow(false,ExcNotImplemented());
    /*
        typename RefBasis::template InterSpaceMap<sdim> dof_map;
        auto sub_ref_space = ref_basis_->template get_ref_sub_space<sdim>(s_id,dof_map,sub_grid);

        IgCoefficients sub_coeffs;
        const int n_sub_dofs = dof_map.size();
        for (int sub_dof = 0 ; sub_dof < n_sub_dofs ; ++ sub_dof)
          sub_coeffs[sub_dof] = coeffs_[dof_map[sub_dof]];

        auto sub_func = IgGridFunction<sdim,space_dim>::const_create(sub_ref_space,sub_coeffs);

        return sub_func;
      //*/
    return nullptr;
  }

public:

  virtual void evaluate_0(const ValueVector<GridPoint> &points,
                          ValueVector<Value> &values) const = 0;

  virtual void evaluate_1(const ValueVector<GridPoint> &points,
                          ValueVector<Derivative<1>> &values) const = 0;

  virtual void evaluate_2(const ValueVector<GridPoint> &points,
                          ValueVector<Derivative<2>> &values) const = 0;


private:

  template<int sdim>
  class SubFormulaGridFunction :
    public GridFunction<sdim,space_dim>
  {
  private:
    using self_t = SubFormulaGridFunction<sdim>;


  public:
    using base_t  = GridFunction<sdim,space_dim>;
    using SupFunc = GridFunction<dim,space_dim>;

    using GridType = Grid<sdim>;
    using SuperGrid = Grid<dim>;

    using typename base_t::ElementAccessor;


    //    template <int j>
    using SubGridMap = typename SuperGrid::template SubGridMap<sdim>;


  public:

    SubFormulaGridFunction(const SharedPtrConstnessHandler<SupFunc> &func,
                           const int s_id,
                           const SubGridMap &sub_grid_elem_map,
                           const SharedPtrConstnessHandler<GridType> &grid)
      :
      base_t(grid),
      func_(func),
      s_id_(s_id),
      elems_property_("boundary"),
      constant_directions_(UnitElement<dim>::template get_elem<sdim>(s_id).constant_directions),
      sub_grid_elem_map_(sub_grid_elem_map)
    {
      LogStream out;
      out.begin_item("Grid:");
      this->get_grid()->print_info(out);
      out.end_item();


      out.begin_item("Sup. Grid:");
      auto sup_grid = std::const_pointer_cast<SuperGrid>(func->get_grid());
      sup_grid->print_info(out);
      out.end_item();


      sup_grid->add_property(elems_property_);
      for (const auto &elem_id : sub_grid_elem_map)
        sup_grid->set_property_status_elem(elems_property_,elem_id.second,true);


      sup_grid->print_info(out);


      auto func_elem = this->begin();
      auto func_elem_end = this->end();
      int elem_id = 0;
      for (; func_elem != func_elem_end ; ++func_elem, ++elem_id)
      {
        out.begin_item("Element " + std::to_string(elem_id));

        out << "Element ID: " << func_elem->get_index() << std::endl;

        out.end_item();
      }



      const auto &constant_values = UnitElement<dim>::template get_elem<sdim>(s_id).constant_values;
      int i = 0;
      for (const auto &v : constant_values)
      {
        const auto dir = constant_directions_[i];
        const auto &knots = grid->get_knot_coordinates(dir);

        constant_coordinates_[i] = (v == 0)? knots.front() : knots.back();
        ++i;
      }

      out.begin_item("Constant coordinates:");
      constant_coordinates_.print_info(out);
      out.end_item();

      AssertThrow(false,ExcNotImplemented());
    }

    auto begin() const
    {
      return func_->begin(elems_property_);
    }

    auto end() const
    {
      return func_->end(elems_property_);
    }


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


    void rebuild_after_insert_knots(
      const SafeSTLArray<SafeSTLVector<double>, sdim> &new_knots,
      const GridType &old_grid) override final
    {
      AssertThrow(false,ExcNotImplemented());
    }

    void print_info(LogStream &out) const override final
    {
      out.begin_item("SubFormulaGridFunction<" + std::to_string(sdim) + ">");

      out.begin_item("Sup. function:");
      func_->print_info(out);
      out.end_item();

      out << "Sub-element ID: " << s_id_ << std::endl;

      out.end_item();
    }

  private:
    SharedPtrConstnessHandler<SupFunc> func_;
    const int s_id_;

    const PropId elems_property_;

    const SafeSTLArray<int,dim-sdim> constant_directions_;

    SafeSTLArray<Real,dim-sdim> constant_coordinates_;

    const SubGridMap sub_grid_elem_map_;

    //  typename SupFunc::ElementIterator sup_elem_;

  };

};

IGA_NAMESPACE_CLOSE

#endif
