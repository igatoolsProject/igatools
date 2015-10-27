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

public:

  virtual void evaluate_0(const ValueVector<GridPoint> &points,
                          ValueVector<Value> &values) const = 0;

  virtual void evaluate_1(const ValueVector<GridPoint> &points,
                          ValueVector<Derivative<1>> &values) const = 0;

  virtual void evaluate_2(const ValueVector<GridPoint> &points,
                          ValueVector<Derivative<2>> &values) const = 0;

#ifdef MESH_REFINEMENT
protected:
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
    const Grid<dim> &old_grid)
  {
    Assert(false,ExcMessage("this function must be specialized in the derved class."));
  }


  void create_connection_for_insert_knots(const std::shared_ptr<self_t> &grid_function);
#endif // MESH_REFINEMENT


};

IGA_NAMESPACE_CLOSE

#endif
