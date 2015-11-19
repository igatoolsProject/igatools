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

#include <igatools/geometry/grid_function_element.h>

IGA_NAMESPACE_OPEN

template<int dim_,int space_dim_>
GridFunctionElement<dim_,space_dim_>::
GridFunctionElement(const std::shared_ptr<ContainerType> &grid_function,
                    const ListIt &index,
                    const PropId &prop)
  :
  grid_function_(grid_function),
  grid_elem_(grid_function_->get_grid()->create_element(index,prop))
{}


template<int dim_,int space_dim_>
bool
GridFunctionElement<dim_,space_dim_>::
same_grid_function_of(const self_t &elem) const
{
  return (grid_function_ == elem.grid_function_);
}


template<int dim_,int space_dim_>
bool
GridFunctionElement<dim_,space_dim_>::
operator ==(const self_t &elem) const
{
  Assert(this->same_grid_function_of(elem),
         ExcMessage("Cannot compare elements on different GridFunction."));
  return (*grid_elem_ == *(elem.grid_elem_));
}



template<int dim_,int space_dim_>
bool
GridFunctionElement<dim_,space_dim_>::
operator !=(const self_t &elem) const
{
  Assert(this->same_grid_function_of(elem),
         ExcMessage("Cannot compare elements on different GridFunction."));
  return (*grid_elem_ != *(elem.grid_elem_));
}



template<int dim_,int space_dim_>
bool
GridFunctionElement<dim_,space_dim_>::
operator <(const self_t &elem) const
{
  Assert(this->same_grid_function_of(elem),
         ExcMessage("Cannot compare elements on different GridFunction."));
  return (*grid_elem_ < *(elem.grid_elem_));
}



template<int dim_,int space_dim_>
bool
GridFunctionElement<dim_,space_dim_>::
operator >(const self_t &elem) const
{
  Assert(this->same_grid_function_of(elem),
         ExcMessage("Cannot compare elements on different GridFunction."));
  return (*grid_elem_ > *(elem.grid_elem_));
}



template<int dim_,int space_dim_>
void
GridFunctionElement<dim_,space_dim_>::
operator++()
{
  ++(*grid_elem_);
}


template<int dim_,int space_dim_>
void
GridFunctionElement<dim_,space_dim_>::
move_to(const IndexType &elem_id)
{
  grid_elem_->move_to(elem_id);
}

template<int dim_,int space_dim_>
auto
GridFunctionElement<dim_,space_dim_>::
get_grid_element() const -> const GridElem &
{
  return *grid_elem_;
}

template<int dim_,int space_dim_>
auto
GridFunctionElement<dim_,space_dim_>::
get_grid_element() -> GridElem &
{
  return *grid_elem_;
}

template<int dim_,int space_dim_>
auto
GridFunctionElement<dim_,space_dim_>::
get_index() const -> const IndexType &
{
  return grid_elem_->get_index();
}


template<int dim_,int space_dim_>
void
GridFunctionElement<dim_,space_dim_>::
print_info(LogStream &out) const
{
  using std::to_string;
  out.begin_item("GridElement<" + to_string(dim_) + "," + to_string(space_dim_) +">");
  grid_elem_->print_info(out);
  out.end_item();
}

template<int dim_,int space_dim_>
void
GridFunctionElement<dim_,space_dim_>::
print_cache_info(LogStream &out) const
{
  out.begin_item("GridElement's cache");
  grid_elem_->print_cache_info(out);
  out.end_item();

  local_cache_.print_info(out);
}


template<int dim_,int space_dim_>
auto
GridFunctionElement<dim_,space_dim_>::
get_element_values_D0() const -> const ValueVector<Value> &
{
  using _D0 = grid_function_element::_D<0>;
  return this->template get_values_from_cache<_D0,dim_>(0);
}

template<int dim_,int space_dim_>
auto
GridFunctionElement<dim_,space_dim_>::
get_element_values_D1() const -> const ValueVector<Derivative<1>> &
{
  using _D1 = grid_function_element::_D<1>;
  return this->template get_values_from_cache<_D1,dim_>(0);
}

template<int dim_,int space_dim_>
auto
GridFunctionElement<dim_,space_dim_>::
get_element_values_D2() const -> const ValueVector<Derivative<2>> &
{
  using _D2 = grid_function_element::_D<2>;
  return this->template get_values_from_cache<_D2,dim_>(0);
}


template<int dim_,int space_dim_>
const ValueVector<Real> &
GridFunctionElement<dim_,space_dim_>::
get_element_weights() const
{
  return grid_elem_->get_element_weights();
}

template<int dim_,int space_dim_>
const ValueVector<Points<dim_>> &
                             GridFunctionElement<dim_,space_dim_>::
                             get_element_points() const
{
  return grid_elem_->get_element_points();
}

#if 0
template<int dim_, int space_dim_>
auto
GridFunctionElement<dim_, space_dim_>::
get_exterior_normals() const -> ValueVector<SafeSTLArray<Value, space_dim_> >
{
  const int sdim = dim_;
  const int s_id = 0;
  Assert(space_dim_ == 1, ExcNotImplemented());
  ValueVector<SafeSTLArray<Value, space_dim_>> res;
  const auto &DF = this->template get_values_from_cache<_Gradient, sdim>(s_id);
  const auto n_points = DF.get_num_points();
  res.resize(n_points);

  for (int pt = 0; pt < n_points; ++pt)
  {
    res[0][pt] = cross_product<dim_, space_dim_>(DF[pt]);
    res[0][pt] /= res[0][pt].norm();
  }

  return res;
}
#endif




IGA_NAMESPACE_CLOSE

#include <igatools/geometry/grid_function_element.inst>

