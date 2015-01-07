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

#include <igatools/basis_functions/bspline_element.h>
#include <igatools/basis_functions/bspline_element_handler.h>
#include <igatools/basis_functions/bernstein_basis.h>
#include <igatools/utils/multi_array_utils.h>

#include <igatools/basis_functions/nurbs_space.h>
#include <igatools/basis_functions/nurbs_element_handler.h>

#include <algorithm>
using std::shared_ptr;

IGA_NAMESPACE_OPEN



template<int dim_, int range_ , int rank_>
std::shared_ptr<ReferenceElementHandler<dim_,range_,rank_> >
ReferenceElementHandler<dim_, range_, rank_>::
create(std::shared_ptr<const Space> space)
{
    using BSplineSp = const BSplineSpace<dim_,range_,rank_>;
    auto bsp_space = std::dynamic_pointer_cast< BSplineSp >(space);

    using NURBSSp = const NURBSSpace<dim_,range_,rank_>;
    auto nrb_space = std::dynamic_pointer_cast< NURBSSp >(space);

    std::shared_ptr<ReferenceElementHandler<dim_,range_,rank_> > elem_handler = nullptr;
    if (bsp_space)
    {
        elem_handler = BSplineElementHandler<dim_,range_,rank_>::create(bsp_space);
    }
    else if (nrb_space)
    {
        elem_handler = NURBSElementHandler<dim_,range_,rank_>::create(nrb_space);
    }
    else
    {
        Assert(false,ExcInvalidState());
        AssertThrow(false,ExcInvalidState());
    }
    return elem_handler;
}



IGA_NAMESPACE_CLOSE

#include <igatools/basis_functions/reference_element_handler.inst>
