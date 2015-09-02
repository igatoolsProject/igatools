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

#include <igatools/geometry/physical_domain.h>
#include <igatools/geometry/physical_domain_element.h>
#include <igatools/functions/identity_function_handler.h>
#include <igatools/functions/function_element.h>
#include <igatools/functions/function_cache_handler.h>

using std::shared_ptr;

IGA_NAMESPACE_OPEN

//template<int dim,int space_dim>
//auto
//IdentityFunctionElementHandler<dim,space_dim>::
//create(std::shared_ptr<GridType> grid) -> std::shared_ptr<parent_t>
//{
//    auto identity_function = std::make_shared<self_t>(grid);
//
//#ifdef MESH_REFINEMENT
//    identity_function->create_connection_for_insert_knots(identity_function);
//#endif // MESH_REFINEMENT
//
//    return identity_function;
//}



//template<int dim,int space_dim>
//auto
//IdentityFunctionElementHandler<dim,space_dim>::
//clone() const -> std::shared_ptr<parent_t>
//{
//
//    return std::make_shared<self_t>(*this);
//}

template<int dim,int space_dim>
IdentityFunctionElementHandler<dim,space_dim>::
IdentityFunctionElementHandler(std::shared_ptr<const GridType> grid)
  : grid_handler_(GridHandlerType::create(grid))
{}


template<int dim,int space_dim>
void
IdentityFunctionElementHandler<dim,space_dim>::
set_flags(const topology_variant &sdim,
          const Flags &flag)
{
//  grid_handler_->set_flags(sdim, grid_flag);
}

template<int dim,int space_dim>
void
IdentityFunctionElementHandler<dim,space_dim>::
init_cache(ConstElementAccessor &elem,
           const eval_pts_variant &quad) const
{
  grid_handler_->init_cache(elem.get_domain_element().get_grid_element(), quad);
}



template<int dim,int space_dim>
auto
IdentityFunctionElementHandler<dim,space_dim>::
fill_cache(const topology_variant &sdim,
           ConstElementAccessor &elem,
           const int s_id) const -> void
{
  grid_handler_->
  fill_cache(sdim, elem.get_domain_element().get_grid_element(),
             s_id);
//
//
//        ->template get_points<sdim>(s_id_);
//  auto fill_dispatcher = FillCacheDispatcher(s_id, *this, elem);
//  boost::apply_visitor(fill_dispatcher, sdim);
}



#ifdef MESH_REFINEMENT
template<int dim,int space_dim>
void
IdentityFunctionElementHandler<dim,space_dim>::
rebuild_after_insert_knots(
  const SafeSTLArray<SafeSTLVector<Real>,dim> &knots_to_insert,
  const CartesianGrid<dim> &grid_old)
{
  using std::const_pointer_cast;

  auto previous_grid = this->get_grid()->get_grid_pre_refinement();
  Assert(previous_grid != nullptr,ExcNullPtr());

  this->function_previous_refinement_ =
    IdentityFunctionElementHandler<dim,space_dim>::create(previous_grid);
}

template<int dim,int space_dim>
void
IdentityFunctionElementHandler<dim,space_dim>::
create_connection_for_insert_knots(std::shared_ptr<self_t> &identity_function)
{
  Assert(identity_function != nullptr, ExcNullPtr());
  Assert(&(*identity_function) == &(*this), ExcMessage("Different objects."));

  using SlotType = typename CartesianGrid<dim>::SignalInsertKnotsSlot;

  auto func_to_connect =
    std::bind(&self_t::rebuild_after_insert_knots,
              identity_function.get(),
              std::placeholders::_1,
              std::placeholders::_2);
  std::const_pointer_cast<CartesianGrid<dim>>(this->get_grid())->connect_insert_knots(
                                             SlotType(func_to_connect).track_foreign(identity_function));
}
#endif // MESH_REFINEMENT


IGA_NAMESPACE_CLOSE

#include <igatools/functions/identity_function_handler.inst>

