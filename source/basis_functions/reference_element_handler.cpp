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

#include <igatools/basis_functions/reference_element_handler.h>

#include <igatools/basis_functions/bspline_element.h>
#include <igatools/basis_functions/bspline_element_handler.h>

#include <igatools/basis_functions/nurbs_space.h>
#include <igatools/basis_functions/nurbs_element_handler.h>

using std::shared_ptr;

IGA_NAMESPACE_OPEN


template<int dim, int range , int rank>
ReferenceElementHandler<dim, range, rank>::
ReferenceElementHandler(shared_ptr<const Space> space)
    :
    grid_handler_(space->get_grid()),
    space_(space)
{
    Assert(space != nullptr, ExcNullPtr());
};


template<int dim, int range , int rank>
shared_ptr<ReferenceElementHandler<dim,range,rank> >
ReferenceElementHandler<dim, range, rank>::
create(shared_ptr<const Space> space)
{
    std::shared_ptr<ReferenceElementHandler<dim,range,rank> > elem_handler = nullptr;
    if (space->is_bspline())
    {
        using BSplineSp = const BSplineSpace<dim,range,rank>;
        auto bsp_space = std::dynamic_pointer_cast< BSplineSp >(space);
        elem_handler = BSplineElementHandler<dim,range,rank>::create(bsp_space);
    }
    else
    {
#ifdef NURBS
        using NURBSSp = const NURBSSpace<dim,range,rank>;
        auto nrb_space = std::dynamic_pointer_cast< NURBSSp >(space);
        elem_handler = NURBSElementHandler<dim,range,rank>::create(nrb_space);
#else
        Assert(false,ExcMessage("NURBS support disabled from configuration cmake parameters."));
        AssertThrow(false,ExcMessage("NURBS support disabled from configuration cmake parameters."));
#endif
    }

    return elem_handler;
}


template<int dim, int range , int rank>
auto
ReferenceElementHandler<dim, range, rank>::
get_space() const -> shared_ptr<const Space>
{
    Assert(space_ != nullptr,ExcNullPtr());
    return space_;
}

template<int dim, int range , int rank>
const GridElementHandler<dim> &
ReferenceElementHandler<dim, range, rank>::
get_grid_handler() const
{
    return this->grid_handler_;
}


template<int dim, int range , int rank>
void
ReferenceElementHandler<dim, range, rank>::
reset(const ValueFlags &flag, const eval_pts_variant &eval_pts)
{
    using ElemProperty = typename CartesianGrid<dim>::ElementProperty;
    const std::set<int> active_elems_id =
        this->get_space()->get_grid()->get_elements_id_same_property(ElemProperty::active);

    this->reset_selected_elements(
        flag,
        eval_pts,
        vector<int>(active_elems_id.begin(),active_elems_id.end()));
}



template<int dim, int range , int rank>
void
ReferenceElementHandler<dim, range, rank>::
reset_one_element(
    const ValueFlags &flag,
    const eval_pts_variant &eval_points,
    const int elem_flat_id)
{
    this->reset_selected_elements(flag,eval_points,vector<int>(1,elem_flat_id));
}



IGA_NAMESPACE_CLOSE

#include <igatools/basis_functions/reference_element_handler.inst>
