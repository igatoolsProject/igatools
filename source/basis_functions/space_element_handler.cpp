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

#include <igatools/basis_functions/space_element_handler.h>


using std::shared_ptr;


IGA_NAMESPACE_OPEN

template<int dim,int codim,int range,int rank>
SpaceElementHandler<dim,codim,range,rank>::
SpaceElementHandler(std::shared_ptr<const Space<dim,codim,range,rank>> space)
    :
    space_(space)
{
    Assert(space != nullptr,ExcNullPtr());
}


template<int dim,int codim,int range,int rank>
void
SpaceElementHandler<dim,codim,range,rank>::
reset(const ValueFlags &flag, const eval_pts_variant &eval_pts)
{
    const std::set<int> elems_id =
        space_->get_grid()->get_elements_id();

    this->reset_selected_elements(
        flag,
        eval_pts,
        SafeSTLVector<int>(elems_id.begin(),elems_id.end()));
}

template<int dim,int codim,int range,int rank>
std::shared_ptr<const Space<dim,codim,range,rank> >
SpaceElementHandler<dim,codim,range,rank>::
get_space() const
{
    return space_;
}


/*
#ifdef SERIALIZATION
template<int dim,int codim,int range,int rank>
template<class Archive>
void
SpaceElementHandler<dim,codim,range,rank>::
serialize(Archive &ar, const unsigned int version)
{
    auto non_const_space = std::const_pointer_cast<Space<dim,codim,range,rank>>(space_);
    ar &boost::serialization::make_nvp("space_", non_const_space);
    space_ = non_const_space;
    Assert(space_ != nullptr,ExcNullPtr());
}
#endif // SERIALIZATION
*/

IGA_NAMESPACE_CLOSE

#include <igatools/basis_functions/space_element_handler.inst>
