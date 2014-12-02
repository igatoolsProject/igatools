//-+--------------------------------------------------------------------
// Igatools a general purpose Isogeometric analysis library.
// Copyright (C) 2012-2014  by the igatools authors (see authors.txt).
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


#include <igatools/basis_functions/nurbs_element.h>
#include <igatools/basis_functions/nurbs_space.h>

#ifdef NURBS

IGA_NAMESPACE_OPEN


template <int dim, int range, int rank>
NURBSElement<dim, range, rank>::
NURBSElement(const std::shared_ptr<ContainerType> space,
             const Index index)
    :
    parent_t(space,index),
    bspline_elem_(space->get_spline_space(),index),
    weight_elem_(space->weight_func_->get_grid(),index)
{}



template <int dim, int range, int rank>
NURBSElement<dim, range, rank>::
NURBSElement(const std::shared_ptr<ContainerType> space,
             const TensorIndex<dim> &index)
    :
    parent_t(space,index),
    bspline_elem_(space->get_spline_space(),index),
    weight_elem_(space->weight_func_->get_grid(),index)
{}


template <int dim, int range, int rank>
void
NURBSElement<dim, range, rank>::
operator++()
{
    parent_t::operator++();
    ++bspline_elem_;
    ++weight_elem_;
}

template <int dim, int range, int rank>
bool
NURBSElement<dim, range, rank>::
jump(const TensorIndex<dim> &increment)
{
    const bool    grid_elem_active =     parent_t::jump(increment);
    const bool bspline_elem_active = bspline_elem_.jump(increment);
    const bool  weight_elem_active =  weight_elem_.jump(increment);

    return grid_elem_active && bspline_elem_active && weight_elem_active;
}

template <int dim, int range, int rank>
void
NURBSElement<dim, range, rank>::
move_to(const Index flat_index)
{
    parent_t::move_to(flat_index);
    bspline_elem_.move_to(flat_index);
    weight_elem_.move_to(flat_index);
}


template <int dim, int range, int rank>
void
NURBSElement<dim, range, rank>::
move_to(const TensorIndex<dim> &tensor_index)
{
    parent_t::move_to(tensor_index);
    bspline_elem_.move_to(tensor_index);
    weight_elem_.move_to(tensor_index);
}

IGA_NAMESPACE_CLOSE

#endif // #ifdef NURBS

//#include <igatools/basis_functions/bspline_element.inst>


