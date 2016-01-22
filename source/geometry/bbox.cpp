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


#include <igatools/geometry/bbox.h>
//#include <igatools/base/exceptions.h>
#include <igatools/geometry/unit_element.h>

IGA_NAMESPACE_OPEN

template<int dim>
BBox<dim>::
BBox()
{
  for (auto &bounds_dir : (*this))
  {
    bounds_dir[0] = 0.0;
    bounds_dir[1] = 1.0;
  }
}

template<int dim>
void
BBox<dim>::
translate(const Points<dim> &translation_amount)
{
  for (int i = 0 ; i < dim ; ++i)
  {
    (*this)[i][0] += translation_amount[i];
    (*this)[i][1] += translation_amount[i];
  }
}

template<int dim>
void
BBox<dim>::
dilate(const Points<dim> &dilation_factor)
{
  for (int i = 0 ; i < dim ; ++i)
  {
    Assert(dilation_factor[i] > 0., ExcMessage("Dilation factor must be positive."));
    (*this)[i][0] *= dilation_factor[i];
    (*this)[i][1] *= dilation_factor[i];
  }
}

template<int dim>
bool
BBox<dim>::
is_unit() const
{
  return std::all_of(this->begin(),this->end(),
                     [](const auto &interv)
  {
    return (interv[0] == 0.0 && interv[1] == 1.0);
  });
}

template<int dim>
SafeSTLArray<Real,dim>
BBox<dim>::
get_side_lengths() const
{
  SafeSTLArray<Real,dim> side_lengths;
  for (int i = 0 ; i < dim ; ++i)
  {
    side_lengths[i] = (*this)[i][1] - (*this)[i][0];
  }
  return side_lengths;
}


template<int dim>
bool
BBox<dim>::
is_point_inside(const Points<dim> &point,const Real eps) const
{
#ifndef NDEBUG
  const auto lengths = this->get_side_lengths();
  Assert(eps >= 0.0,ExcLowerRange(eps,0.0));
  for (int i = 0 ; i < dim ; ++i)
  {
    const auto half_length = lengths[i]/2.0;
    Assert(eps < half_length,
           ExcMessage("The \"eps\" parameter must be smaller of " + std::to_string(half_length) +
                      " (problem occurred along direction " + std::to_string(i) + ")"));
  }
#endif

  bool is_point_inside = true;
  for (int i = 0 ; i < dim ; ++i)
  {
    if (point[i] < ((*this)[i][0]-eps) || point[i] > ((*this)[i][1])+eps)
    {
      is_point_inside = false;
      break;
    }
  }
  return is_point_inside;

}


IGA_NAMESPACE_CLOSE



#include <igatools/geometry/bbox.inst>

