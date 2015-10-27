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

#include <igatools/geometry/formula_grid_function.h>
#include <igatools/geometry/formula_grid_function_handler.h>

IGA_NAMESPACE_OPEN

template<int dim, int space_dim>
FormulaGridFunction<dim, space_dim>::
FormulaGridFunction(const SharedPtrConstnessHandler<GridType> &grid)
  :
  parent_t::GridFunction(grid)
{}



template<int dim, int space_dim>
auto
FormulaGridFunction<dim, space_dim>::
create_cache_handler() const
-> std::unique_ptr<typename parent_t::ElementHandler>
{
  return std::make_unique<ElementHandler>(
    std::dynamic_pointer_cast<const self_t>(this->shared_from_this()));
}




#ifdef MESH_REFINEMENT
template<int dim, int space_dim>
void
FormulaGridFunction<dim,space_dim>::
create_connection_for_insert_knots(const std::shared_ptr<self_t> &grid_function)
{
  Assert(grid_function != nullptr, ExcNullPtr());
  Assert(&(*grid_function) == &(*this), ExcMessage("Different objects."));

  auto func_to_connect =
    std::bind(&self_t::rebuild_after_insert_knots,
              grid_function.get(),
              std::placeholders::_1,
              std::placeholders::_2);

  using SlotType = typename Grid<dim>::SignalInsertKnotsSlot;
  this->get_grid()
  ->connect_insert_knots(SlotType(func_to_connect).track_foreign(grid_function));
}


#endif // MESH_REFINEMENT

IGA_NAMESPACE_CLOSE

#include <igatools/geometry/formula_grid_function.inst>

