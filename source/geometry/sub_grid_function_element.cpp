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

#include <igatools/geometry/sub_grid_function_element.h>

IGA_NAMESPACE_OPEN


template <int sdim,int dim,int range>
SubGridFunctionElement<sdim,dim,range>::
SubGridFunctionElement(const std::shared_ptr<ContainerType> &sub_grid_function,
                       std::unique_ptr<GridElement<sdim>> &&sub_grid_element,
                       std::unique_ptr<GridFunctionElement<dim,range>> &&sup_grid_func_element)
  :
  parent_t(sub_grid_function,std::move(sub_grid_element)),
  sup_grid_func_element_(std::move(sup_grid_func_element))
{}

template <int sdim,int dim,int range>
bool
SubGridFunctionElement<sdim,dim,range>::
operator==(const parent_t &elem) const
{
  const self_t &sub_elem = dynamic_cast<const self_t &>(elem);
  Assert(this->same_grid_function_of(elem) &&
         sup_grid_func_element_->same_grid_function_of(*(sub_elem.sup_grid_func_element_)),
         ExcMessage("Cannot compare elements on different GridFunction."));

  return this->parent_t::operator==(elem);
}

template <int sdim,int dim,int range>
bool
SubGridFunctionElement<sdim,dim,range>::
operator!=(const parent_t &elem) const
{
  const self_t &sub_elem = dynamic_cast<const self_t &>(elem);
  Assert(this->same_grid_function_of(elem) &&
         sup_grid_func_element_->same_grid_function_of(*(sub_elem.sup_grid_func_element_)),
         ExcMessage("Cannot compare elements on different GridFunction."));
  return this->parent_t::operator!=(elem);
}



template <int sdim,int dim,int range>
void
SubGridFunctionElement<sdim,dim,range>::
operator++()
{
  parent_t::operator++();

  using SubGridFunc = SubGridFunction<sdim,dim,range>;
  const auto grid_func =
    std::dynamic_pointer_cast<const SubGridFunc>(this->grid_function_);

#if 0
  LogStream out;
  out.begin_item("operator++");

  out.begin_item("Sub-Grid");
  this->get_grid_element().get_grid()->print_info(out);
  out.end_item();

  out.begin_item("Sup-Grid");
  sup_grid_func_element_->get_grid_element().get_grid()->print_info(out);
  out.end_item();

  out.end_item();
#endif
  if (*this->get_grid_element().get_index_iterator() != *grid_func->get_id_elems_sub_grid().end())
  {
    const auto &sub_elem_id = this->get_index();
    const auto &sup_elem_id = grid_func->get_sup_element_id(sub_elem_id);

    sup_grid_func_element_->move_to(sup_elem_id);
//      out << "Sub elem ID: " << sub_elem_id << "    Sup elem ID: " << sup_elem_id << std::endl;
  }
  else
  {
//      const auto & sup_elem_id = *(grid_func->get_id_elems_sup_grid().end());
//      sup_grid_func_element_->move_to(sup_elem_id);
    sup_grid_func_element_->move_to(*(--grid_func->get_id_elems_sup_grid().end()));
    ++(*sup_grid_func_element_);
//      out << "Sub elem ID: " << this->get_index() << "    Sup elem ID: " << sup_grid_func_element_->get_index() << std::endl;
  }

//    Assert(false,ExcNotImplemented());
}


template <int sdim,int dim,int range>
void
SubGridFunctionElement<sdim,dim,range>::
move_to(const IndexType &elem_id)
{
  parent_t::move_to(elem_id);

  using SubGridFunc = SubGridFunction<sdim,dim,range>;
  const auto grid_func =
    std::dynamic_pointer_cast<const SubGridFunc>(this->grid_function_);

  if (*this->get_grid_element().get_index_iterator() != *grid_func->get_id_elems_sub_grid().end())
  {
    const auto &sub_elem_id = this->get_index();
    const auto &sup_elem_id = grid_func->get_sup_element_id(sub_elem_id);

    sup_grid_func_element_->move_to(sup_elem_id);
  }
  else
  {
    sup_grid_func_element_->move_to(*(--grid_func->get_id_elems_sup_grid().end()));
    ++(*sup_grid_func_element_);
  }
}



template <int sdim,int dim,int range>
void
SubGridFunctionElement<sdim,dim,range>::
print_info(LogStream &out) const
{
  using std::to_string;
  out.begin_item("SubGridFunctionElement<" +
                 to_string(sdim) + "," +
                 to_string(dim) + "," +
                 to_string(range) + ">");

  out.begin_item("GridFunctionElement<" + to_string(sdim) + "," + to_string(range) + ">");
  parent_t::print_info(out);
  out.end_item();

  out.begin_item("Sup-GridFunctionElement<" + to_string(dim) + "," + to_string(range) + ">");
  sup_grid_func_element_->print_info(out);
  out.end_item();

  out.end_item();
}



template <int sdim,int dim,int range>
GridFunctionElement<dim,range> &
SubGridFunctionElement<sdim,dim,range>::
get_sup_grid_function_element()
{
  return *sup_grid_func_element_;
}

IGA_NAMESPACE_CLOSE

#include <igatools/geometry/sub_grid_function_element.inst>
