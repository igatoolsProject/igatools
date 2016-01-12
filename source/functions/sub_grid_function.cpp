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

#include <igatools/functions/sub_grid_function.h>
#include <igatools/functions/sub_grid_function_element.h>
#include <igatools/functions/sub_grid_function_handler.h>

using std::shared_ptr;

IGA_NAMESPACE_OPEN

template<int sdim,int dim,int range>
SubGridFunction<sdim,dim,range>::
SubGridFunction(const SharedPtrConstnessHandler<SupFunc> &sup_func,
                const int s_id,
                const SubGridMap &sub_grid_elem_map,
                const SharedPtrConstnessHandler<GridType> &grid)
  :
  base_t(grid),
  sup_func_(sup_func),
  s_id_(s_id)//,
//  elems_property_("boundary"),
//  sub_grid_elem_map_(sub_grid_elem_map)
{
//    LogStream out;
  for (const auto &elems_id : sub_grid_elem_map)
  {
    id_elems_sub_grid_.emplace_back(elems_id.first);
    id_elems_sup_grid_.emplace_back(elems_id.second);
//      out << "Sub elem ID: " << elems_id.first << "    Sup elem ID: " << elems_id.second << std::endl;
  }
}


template<int sdim,int dim,int range>
auto
SubGridFunction<sdim,dim,range>::
create_element_begin(const PropId &prop) const
-> std::unique_ptr<GridFunctionElement<sdim,range>>
{
  using Elem = SubGridFunctionElement<sdim,dim,range>;
  return std::unique_ptr<Elem>(new Elem(
    std::dynamic_pointer_cast<const self_t>(this->shared_from_this()),
    std::move(this->get_grid()->create_element_begin(prop)),
    std::move(sup_func_->create_element_begin(prop))));
}

template<int sdim,int dim,int range>
auto
SubGridFunction<sdim,dim,range>::
create_element_end(const PropId &prop) const
-> std::unique_ptr<GridFunctionElement<sdim,range>>
{
  using Elem = SubGridFunctionElement<sdim,dim,range>;
  return std::unique_ptr<Elem>(new Elem(
    std::dynamic_pointer_cast<const self_t>(this->shared_from_this()),
    std::move(this->get_grid()->create_element_end(prop)),
    std::move(sup_func_->create_element_end(prop))));
}

template<int sdim,int dim,int range>
GridIterator<GridFunctionElement<sdim,range> >
SubGridFunction<sdim,dim,range>::
cbegin(const PropId &prop) const
{
  using ElemIt = GridIterator<GridFunctionElement<sdim,range> >;
  return ElemIt(this->create_element_begin(prop));
}

template<int sdim,int dim,int range>
GridIterator<GridFunctionElement<sdim,range> >
SubGridFunction<sdim,dim,range>::
cend(const PropId &prop) const
{
  using ElemIt = GridIterator<GridFunctionElement<sdim,range> >;
  return ElemIt(this->create_element_end(prop));
}


template<int sdim,int dim,int range>
auto
SubGridFunction<sdim,dim,range>::
const_create(const std::shared_ptr<const SupFunc> &func,
             const int s_id,
             const SubGridMap &sub_grid_elem_map,
             const std::shared_ptr<const GridType> &grid)
-> std::shared_ptr<const self_t>
{
  return std::make_shared<self_t>(SharedPtrConstnessHandler<SupFunc>(func),
  s_id,
  sub_grid_elem_map,
  SharedPtrConstnessHandler<GridType>(grid));
}

template<int sdim,int dim,int range>
auto
SubGridFunction<sdim,dim,range>::
create(const std::shared_ptr<const SupFunc> &func,
       const int s_id,
       const SubGridMap &sub_grid_elem_map,
       const std::shared_ptr<GridType> &grid)
-> std::shared_ptr<self_t>
{
  return std::make_shared<self_t>(SharedPtrConstnessHandler<SupFunc>(func),
  s_id,
  sub_grid_elem_map,
  SharedPtrConstnessHandler<GridType>(grid));
}

template<int sdim,int dim,int range>
std::unique_ptr<GridFunctionHandler<sdim,range> >
SubGridFunction<sdim,dim,range>::
create_cache_handler() const
{
  using Handler = SubGridFunctionHandler<sdim,dim,range>;
  return std::unique_ptr<Handler>(
           new Handler(std::dynamic_pointer_cast<const self_t>(this->shared_from_this())
                      ));
}

#if 0
template<int sdim,int dim,int range>
std::unique_ptr<GridFunctionElement<sdim,range> >
SubGridFunction<sdim,dim,range>::
create_element(const ListIt &index, const PropId &prop) const
{
#if 0
  using Elem = SubGridFunctionElement<sdim,dim,range>;
  return std::unique_ptr<Elem>(
           new Elem(std::dynamic_pointer_cast<const self_t>(this->shared_from_this()),
                    index, prop));
#endif

  Assert(false,ExcNotImplemented());
  return nullptr;
}
#endif

#ifdef MESH_REFINEMENT
template<int sdim,int dim,int range>
void
SubGridFunction<sdim,dim,range>::
rebuild_after_insert_knots(
  const SafeSTLArray<SafeSTLVector<double>, sdim> &new_knots,
  const GridType &old_grid)
{
  AssertThrow(false,ExcNotImplemented());
}
#endif


template<int sdim,int dim,int range>
void
SubGridFunction<sdim,dim,range>::
print_info(LogStream &out) const
{
  out.begin_item("SubGridFunction<" + std::to_string(sdim) + ">");

  out.begin_item("Sup. function:");
  sup_func_->print_info(out);
  out.end_item();


  out << "Sub-element topology ID: " << s_id_ << std::endl;
  /*
    out.begin_item("Sub-Grid Element Map:");
    sub_grid_elem_map_.print_info(out);
    out.end_item();
  //*/
  out.begin_item("Sub elements ID:");
  id_elems_sub_grid_.print_info(out);
  out.end_item();

  out.begin_item("Sup elements ID:");
  id_elems_sup_grid_.print_info(out);
  out.end_item();

  out.end_item();
}


template<int sdim,int dim,int range>
auto
SubGridFunction<sdim,dim,range>::
get_sup_grid_function() const -> std::shared_ptr<const SupFunc>
{
  return sup_func_.get_ptr_const_data();
}


template<int sdim,int dim,int range>
const SafeSTLVector<typename Grid<sdim>::IndexType> &
SubGridFunction<sdim,dim,range>::
get_id_elems_sub_grid() const
{
  return id_elems_sub_grid_;
}

template<int sdim,int dim,int range>
const SafeSTLVector<typename Grid<dim>::IndexType> &
SubGridFunction<sdim,dim,range>::
get_id_elems_sup_grid() const
{
  return id_elems_sup_grid_;
}

template<int sdim,int dim,int range>
const typename Grid<dim>::IndexType &
SubGridFunction<sdim,dim,range>::
get_sup_element_id(const typename Grid<sdim>::IndexType &sub_elem_id) const
{
  const auto begin = id_elems_sub_grid_.cbegin();
  const auto end   = id_elems_sub_grid_.cend();

  const auto it = std::find(begin,end,sub_elem_id);
  Assert(it != end,ExcMessage("Index not found."));

  return id_elems_sup_grid_[it-begin];

//  return sub_grid_elem_map_.at(sub_elem_id);
}

/*
template<int sdim,int dim,int range>
auto
SubGridFunction<sdim,dim,range>::
get_sub_grid_elem_map() const -> const SubGridMap &
{
  return sub_grid_elem_map_;
}
//*/

IGA_NAMESPACE_CLOSE

#include <igatools/functions/sub_grid_function.inst>

