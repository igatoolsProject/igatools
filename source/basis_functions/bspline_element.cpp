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


#include <igatools/basis_functions/bspline.h>
#include <igatools/basis_functions/bspline_element.h>
#include <igatools/basis_functions/bspline_handler.h>
#include <igatools/basis_functions/bernstein_basis.h>

#include <igatools/utils/multi_array_utils.h>

#include <algorithm>
#include <numeric>
#include <memory>

using std::reverse;
using std::accumulate;
using std::sort;

using std::shared_ptr;
using std::make_shared;



IGA_NAMESPACE_OPEN



template <int dim, int range, int rank>
BSplineElement<dim, range, rank>::
BSplineElement(const std::shared_ptr<ContainerType> &basis,
               std::unique_ptr<GridElement<dim>> &&grid_elem)
  :
  parent_t(basis),
  grid_elem_(std::move(grid_elem))
{}






template <int dim, int range, int rank>
auto
BSplineElement<dim, range, rank>::
get_bspline_basis() const -> std::shared_ptr<const Basis>
{
  auto bsp_basis = std::dynamic_pointer_cast<const Basis>(this->basis_);
  Assert(bsp_basis != nullptr,ExcNullPtr());
  return bsp_basis;
}



template <int dim, int range, int rank>
void
BSplineElement<dim, range, rank>::
print_info(LogStream &out) const
{
  using std::to_string;

  out.begin_item("ReferenceBasisElement<" +
                 to_string(dim) + "," +
                 to_string(range) + "," +
                 to_string(rank) + ">");
  parent_t::print_info(out);
  out.end_item();

  out.begin_item("GridElement<" + to_string(dim) + ">");
  grid_elem_->print_info(out);
  out.end_item();
}


template <int dim, int range, int rank>
void
BSplineElement<dim, range, rank>::
print_cache_info(LogStream &out) const
{
  using std::to_string;
  out.begin_item("BSplineElement<" +
                 to_string(dim) + "," +
                 to_string(range) + "," +
                 to_string(rank) + "> cache:");

  out.begin_item("Splines 1D table:");
  for (int sdim = 0 ; sdim <= dim ; ++sdim)
  {
    out.begin_item("Sub-element dimension: " + to_string(sdim));
    all_splines_1D_table_[sdim].print_info(out);
    out.end_item();
  }
  out.end_item();
//*/
  out.begin_item("BasisElement's cache:");
  BasisElement<dim,0,range,rank>::print_cache_info(out);
  out.end_item();

  out.end_item();
}



template <int dim, int range, int rank>
auto
BSplineElement<dim, range, rank>::
get_grid_element() -> GridElem &
{
  return *grid_elem_;
}

template <int dim, int range, int rank>
auto
BSplineElement<dim, range, rank>::
get_grid_element() const -> const GridElem &
{
  return *grid_elem_;
}

template <int dim, int range, int rank>
void
BSplineElement<dim, range, rank>::
operator++()
{
  ++(*grid_elem_);
}

template <int dim, int range, int rank>
void
BSplineElement<dim, range, rank>::
move_to(const IndexType &elem_id)
{
  grid_elem_->move_to(elem_id);
}


IGA_NAMESPACE_CLOSE

#include <igatools/basis_functions/bspline_element.inst>


