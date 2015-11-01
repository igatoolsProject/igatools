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

#include <igatools/basis_functions/reference_element_handler.h>

#include <igatools/basis_functions/bspline.h>
#include <igatools/basis_functions/bspline_element_handler.h>

#include <igatools/basis_functions/nurbs.h>
#include <igatools/basis_functions/nurbs_element_handler.h>

using std::shared_ptr;

IGA_NAMESPACE_OPEN


template<int dim, int range , int rank>
ReferenceElementHandler<dim, range, rank>::
ReferenceElementHandler(const shared_ptr<const Space> &space)
  :
  base_t(space),
  grid_handler_(space->get_grid())
{};


#if 0
template<int dim, int range , int rank>
shared_ptr<ReferenceElementHandler<dim,range,rank> >
ReferenceElementHandler<dim, range, rank>::
create(const shared_ptr<const Space> &space)
{
  std::shared_ptr<ReferenceElementHandler<dim,range,rank> > elem_handler = nullptr;
  if (space->is_bspline())
  {
    using BSplineSp = const BSpline<dim,range,rank>;
    auto bsp_space = std::dynamic_pointer_cast< BSplineSp >(space);
    elem_handler = BSplineElementHandler<dim,range,rank>::create(bsp_space);
  }
  else
  {
#ifdef USE_NURBS
    using NURBSSp = const NURBS<dim,range,rank>;
    auto nrb_space = std::dynamic_pointer_cast< NURBSSp >(space);
    elem_handler = NURBSElementHandler<dim,range,rank>::create(nrb_space);
#else
    Assert(false,ExcMessage("NURBS support disabled from configuration cmake parameters."));
    AssertThrow(false,ExcMessage("NURBS support disabled from configuration cmake parameters."));
#endif
  }

  return elem_handler;
}
#endif

#if 0
template<int dim, int range , int rank>
const GridHandler<dim> &
ReferenceElementHandler<dim, range, rank>::
get_grid_handler() const
{
  return this->grid_handler_;
}
#endif

#if 0
template<int dim, int range , int rank>
void
ReferenceElementHandler<dim, range, rank>::
init_cache(SpaceElement<dim,0,range,rank,Transformation::h_grad> &space_elem,
           const topology_variant &topology)
{
  auto &ref_elem = dynamic_cast<ElementAccessor &>(space_elem);
  this->init_ref_elem_cache(ref_elem,topology);
}

template<int dim, int range , int rank>
void
ReferenceElementHandler<dim, range, rank>::
fill_cache(SpaceElement<dim,0,range,rank,Transformation::h_grad> &space_elem,
           const topology_variant &topology,
           const int sub_elem_id)
{
  auto &ref_elem = dynamic_cast<ElementAccessor &>(space_elem);
  this->fill_ref_elem_cache(ref_elem,topology,sub_elem_id);
}
#endif

#if 0
template<int dim, int range , int rank>
void
ReferenceElementHandler<dim, range, rank>::
reset_one_element(
  const ValueFlags &flag,
  const eval_pts_variant &eval_points,
  const int elem_flat_id)
{
  this->reset_selected_elements(flag,eval_points,SafeSTLVector<int>(1,elem_flat_id));
}
#endif


IGA_NAMESPACE_CLOSE

#include <igatools/basis_functions/reference_element_handler.inst>
