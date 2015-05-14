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


#include <igatools/basis_functions/space_element_base.h>



IGA_NAMESPACE_OPEN

template <int dim>
SpaceElementBase<dim>::
SpaceElementBase(const std::shared_ptr<const SpaceBase<dim>> space,
                 const Index elem_index)
    :
    base_t(space->get_grid(), elem_index),
    space_(space)
{
    Assert(space_ != nullptr, ExcNullPtr());
}

template <int dim>
SpaceElementBase<dim>::
SpaceElementBase(const self_t &elem,
                 const CopyPolicy &copy_policy)
    :
    base_t(elem,copy_policy),
    space_(elem.space_)
{}

template <int dim>
CartesianGridElement<dim> &
SpaceElementBase<dim>::
as_cartesian_grid_element_accessor()
{
    return static_cast<CartesianGridElement<dim> &>(*this);
}

template <int dim>
const CartesianGridElement<dim> &
SpaceElementBase<dim>::
as_cartesian_grid_element_accessor() const
{
    return static_cast<const CartesianGridElement<dim> &>(*this);
}

template <int dim>
void
SpaceElementBase<dim>::
print_info(LogStream &out) const
{
    base_t::print_info(out);

    out.begin_item("Element global connectivity (property=\"" + DofProperties::active + "\"):");
    const auto glob_dofs = this->get_local_to_global(DofProperties::active);
    glob_dofs.print_info(out);
    out.end_item();
}


template <int dim>
void
SpaceElementBase<dim>::
print_cache_info(LogStream &out) const
{
    base_t::print_cache_info(out);
}

template <int dim>
void
SpaceElementBase<dim>::
copy_from(const SpaceElementBase<dim> &elem,
          const CopyPolicy &copy_policy)
{
    if (this != &elem)
    {
        CartesianGridElement<dim>::copy_from(elem,copy_policy);
        space_ = elem.space_;
    }
}

template <int dim>
SafeSTLVector<Index>
SpaceElementBase<dim>::
get_local_to_global(const std::string &dofs_property) const
{
    SafeSTLVector<Index> dofs_global;
    SafeSTLVector<Index> dofs_loc_to_patch;
    SafeSTLVector<Index> dofs_loc_to_elem;
    this->space_->get_element_dofs(
        this->get_flat_index(),
        dofs_global,dofs_loc_to_patch,dofs_loc_to_elem,dofs_property);

    return dofs_global;
}

template <int dim>
SafeSTLVector<Index>
SpaceElementBase<dim>::
get_local_to_patch(const std::string &dofs_property) const
{
    SafeSTLVector<Index> dofs_global;
    SafeSTLVector<Index> dofs_loc_to_patch;
    SafeSTLVector<Index> dofs_loc_to_elem;
    this->space_->get_element_dofs(
        this->get_flat_index(),
        dofs_global,dofs_loc_to_patch,dofs_loc_to_elem,dofs_property);

    return dofs_loc_to_patch;
}

template <int dim>
SafeSTLVector<Index>
SpaceElementBase<dim>::
get_local_dofs(const std::string &dofs_property) const
{
    SafeSTLVector<Index> dofs_global;
    SafeSTLVector<Index> dofs_loc_to_patch;
    SafeSTLVector<Index> dofs_loc_to_elem;
    this->space_->get_element_dofs(
        this->get_flat_index(),
        dofs_global,dofs_loc_to_patch,dofs_loc_to_elem,dofs_property);

    return dofs_loc_to_elem;
}

template <int dim>
Size
SpaceElementBase<dim>::
get_num_basis(const std::string &dofs_property) const
{
    const auto dofs_global = this->get_local_to_global(dofs_property);
    return dofs_global.size();
}

template <int dim>
bool
SpaceElementBase<dim>::
operator==(const self_t &a) const
{
    Assert(space_ == a.space_,
           ExcMessage("Comparison between elements defined on different spaces"));
    return this->as_cartesian_grid_element_accessor() == a.as_cartesian_grid_element_accessor();
}

template <int dim>
bool
SpaceElementBase<dim>::
operator!=(const self_t &a) const
{
    Assert(space_ == a.space_,
           ExcMessage("Comparison between elements defined on different spaces"));
    return this->as_cartesian_grid_element_accessor() != a.as_cartesian_grid_element_accessor();
}

template <int dim>
bool
SpaceElementBase<dim>::
operator<(const self_t &a) const
{
    Assert(space_ == a.space_,
           ExcMessage("Comparison between elements defined on different spaces"));
    return this->as_cartesian_grid_element_accessor() < a.as_cartesian_grid_element_accessor();
}

template <int dim>
bool
SpaceElementBase<dim>::
operator>(const self_t &a) const
{
    Assert(space_ == a.space_,
           ExcMessage("Comparison between elements defined on different spaces"));
    return this->as_cartesian_grid_element_accessor() > a.as_cartesian_grid_element_accessor();
}

#ifdef SERIALIZATION
template <int dim>
template<class Archive>
void
SpaceElementBase<dim>::
serialize(Archive &ar, const unsigned int version)
{
    ar &boost::serialization::make_nvp("SpaceElementBase_base_t",
                                       boost::serialization::base_object<CartesianGridElement<dim>>(*this));

    auto non_const_space = std::const_pointer_cast<SpaceBase<dim>>(space_);
    ar &boost::serialization::make_nvp("space_",non_const_space);
    space_ = non_const_space;
    Assert(space_ != nullptr,ExcNullPtr());
}
///@}
#endif // SERIALIZATION

IGA_NAMESPACE_CLOSE


#include <igatools/basis_functions/space_element_base.inst>

