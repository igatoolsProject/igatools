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


#ifndef SPACE_ELEMENT_INLINE_H_
#define SPACE_ELEMENT_INLINE_H_

// TODO (pauletti, Oct 20, 2014): the inline name is missleading,
// there is no need for inlining, it is only ignorance on how to
// provide proper instantiation, needs to be improved

#include <igatools/basis_functions/space_element.h>

IGA_NAMESPACE_OPEN

template<class Space>
inline
SpaceElement<Space>::
SpaceElement(const std::shared_ptr<const Space> space,
             const Index elem_index)
    :
    base_t(space, elem_index),
    space_(space)
{}




template<class Space>
inline
SpaceElement<Space>::
SpaceElement(const SpaceElement<Space> &elem,
             const CopyPolicy &copy_policy)
    :
    base_t(elem,copy_policy),
    space_(elem.space_)
{
    Assert(this->space_ != nullptr,ExcNullPtr());

    if (elem.local_cache_ != nullptr)
    {
        if (copy_policy == CopyPolicy::shallow)
        {
            local_cache_ = elem.local_cache_;
        }
        else
        {
            local_cache_ = std::shared_ptr<LocalCache>(new LocalCache(*elem.local_cache_));
        }
    }
}



template<class Space>
void
SpaceElement<Space>::
copy_from(const SpaceElement<Space> &elem,
          const CopyPolicy &copy_policy)
{
    if (this != &elem)
    {
    	SpaceElementBase<Space::dim>::copy_from(elem,copy_policy);

        this->space_ = elem.space_;

        if (copy_policy == CopyPolicy::deep)
        {
            Assert(elem.local_cache_ != nullptr, ExcNullPtr());
            local_cache_ = std::shared_ptr<LocalCache>(new LocalCache(*elem.local_cache_));
        }
        else if (copy_policy == CopyPolicy::shallow)
        {
            local_cache_ = elem.local_cache_;
        }
        else
        {
            Assert(false,ExcNotImplemented());
            AssertThrow(false,ExcNotImplemented());
        }
    }
}



template<class Space>
void
SpaceElement<Space>::
deep_copy_from(const SpaceElement<Space> &elem)
{
    this->copy_from(elem,CopyPolicy::deep);
}



template<class Space>
void
SpaceElement<Space>::
shallow_copy_from(const SpaceElement<Space> &elem)
{
    this->copy_from(elem,CopyPolicy::shallow);
}



template<class Space>
SpaceElement<Space> &
SpaceElement<Space>::
operator=(const SpaceElement<Space> &element)
{
    this->shallow_copy_from(element);
    return (*this);
}



#if 0
template<class Space>
inline
Size
SpaceElement<Space>::
get_num_basis(const std::string &dofs_property) const
{
    const auto dofs_global = this->get_local_to_global(dofs_property);
    return dofs_global.size();
}
#endif

#if 0
template<class Space>
inline
auto
SpaceElement<Space>::
get_local_to_global(const std::string &dofs_property) const -> vector<Index>
{
//    return space_->get_loc_to_global(*this,dofs_property);

    vector<Index> dofs_global;
    vector<Index> dofs_loc_to_patch;
    vector<Index> dofs_loc_to_elem;
    this->space_->get_element_dofs(
    		this->as_cartesian_grid_element_accessor(),
			dofs_global,dofs_loc_to_patch,dofs_loc_to_elem,dofs_property);

    return dofs_global;
}
#endif

#if 0
template<class Space>
inline
auto
SpaceElement<Space>::
get_local_to_patch(const std::string &dofs_property) const -> vector<Index>
{
//    return space_->get_loc_to_patch(*this,dofs_property);

    vector<Index> dofs_global;
    vector<Index> dofs_loc_to_patch;
    vector<Index> dofs_loc_to_elem;
    this->space_->get_element_dofs(
    		this->as_cartesian_grid_element_accessor(),
			dofs_global,dofs_loc_to_patch,dofs_loc_to_elem,dofs_property);

    return dofs_loc_to_patch;
}
#endif

template<class Space>
inline
auto
SpaceElement<Space>::
get_space() const -> std::shared_ptr<const Space>
{
    Assert(this->space_ != nullptr,ExcNullPtr());
    return this->space_;
}



template<class Space>
inline
void
SpaceElement<Space>::
ValuesCache::
resize(const FunctionFlags &flags_handler,
       const Size total_n_points,
       const Size total_n_basis)
{
    flags_handler_ = flags_handler;

    Assert(total_n_points >= 0, ExcLowerRange(total_n_points,1));
    Assert(total_n_basis > 0, ExcLowerRange(total_n_basis,1));

    if (flags_handler_.fill_values())
        resize_der<0>(total_n_basis,total_n_points);
    if (flags_handler_.fill_gradients())
        resize_der<1>(total_n_basis,total_n_points);
    if (flags_handler_.fill_hessians())
        resize_der<2>(total_n_basis,total_n_points);

#if 0
    if (flags_handler_.fill_divergences())
    {
        Assert(flags_handler_.fill_gradients(),
               ExcMessage("Divergence requires gradient to be filled."));

        if (div_phi_.get_num_points() != total_n_points ||
            div_phi_.get_num_functions() != total_n_basis)
        {
            div_phi_.resize(total_n_basis,total_n_points);
            div_phi_.zero();
        }
    }
    else
    {
        div_phi_.clear();
    }
#endif

    this->set_initialized(true);
}



template<class Space>
inline
auto
SpaceElement<Space>::
ValuesCache::print_info(LogStream &out) const -> void
{
    out.begin_item("Fill flags:");
    flags_handler_.print_info(out);
    out.end_item();

    if (flags_handler_.values_filled())
    {
        out.begin_item("Values:");
        get_der<0>().print_info(out);
        out.end_item();
    }

    if (flags_handler_.gradients_filled())
    {
        out.begin_item("Gradients:");
        get_der<1>().print_info(out);
        out.end_item();
    }

    if (flags_handler_.hessians_filled())
    {
        out.begin_item("Hessians:");
        get_der<2>().print_info(out);
        out.end_item();
    }

//    out.begin_item("Divergeces:");
//    div_phi_.print_info(out);
//    out.end_item();
}

template<class Space>
void
SpaceElement<Space>::print_info(LogStream &out) const
{
    base_t::print_info(out);
}

template<class Space>
void
SpaceElement<Space>::
print_cache_info(LogStream &out) const
{
    base_t::print_cache_info(out);

    Assert(local_cache_ != nullptr,ExcNullPtr());
    local_cache_->print_info(out);
}


template<class Space>
void
SpaceElement<Space>::LocalCache::
print_info(LogStream &out) const
{
    cacheutils::print_caches(values_, out);

//    out.begin_item("Element Cache:");
//    get_value_cache<0>(0).print_info(out);
//    out.end_item();
//
//    for (int i = 0 ; i < n_faces ; ++i)
//    {
//        out.begin_item("Face "+ std::to_string(i) + " Cache:");
//        get_value_cache<1>(i).print_info(out);
//        out.end_item();
//    }

}


IGA_NAMESPACE_CLOSE

#endif
