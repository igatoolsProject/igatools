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

#include <igatools/functions/function.h>
#include <igatools/functions/function_element.h>
#include <igatools/utils/unique_id_generator.h>
#include <igatools/geometry/mapping.h>

IGA_NAMESPACE_OPEN

template<int dim_, int codim_, int range_, int rank_ >
Function<dim_, codim_, range_, rank_ >::
Function(std::shared_ptr<GridType> grid)
    :
    GridElementHandler<dim_>(grid),
    object_id_(UniqueIdGenerator::get_unique_id()),
    grid_(std::const_pointer_cast<CartesianGrid<dim_>>(grid))
{
    Assert(grid != nullptr,ExcNullPtr());
}

template<int dim_, int codim_, int range_, int rank_>
Function<dim_, codim_, range_, rank_ >::
Function(const self_t &func)
    :
    GridElementHandler<dim_>(func),
    object_id_(UniqueIdGenerator::get_unique_id()),
    grid_(func.grid_)
{}


template<int dim_, int codim_, int range_, int rank_>
void
Function<dim_, codim_, range_, rank_ >::
set_flags(const topology_variant &sdim,
          const typename ElementAccessor::Flags &flag)
{
    auto set_flags_dispatcher = SetFlagsDispatcher(flag,  phys_domain_, flags_);
    boost::apply_visitor(set_flags_dispatcher, sdim);
}



template<int dim_, int codim_, int range_, int rank_>
void
Function<dim_, codim_, range_, rank_ >::
init_cache(ElementAccessor &elem,
           const eval_pts_variant &quad) const
{
    auto init_dispatcher = InitCacheDispatcher(*this, elem);
    boost::apply_visitor(init_dispatcher, quad);
}



template<int dim_, int codim_, int range_, int rank_>
void
Function<dim_, codim_, range_, rank_ >::
fill_cache(const topology_variant &sdim,
           ElementAccessor &elem,
           const int s_id) const
{
    auto fill_dispatcher = FillCacheDispatcher(s_id, *this, func_elem);
    boost::apply_visitor(fill_dispatcher, sdim);
}



#if 0
template<int dim_, int codim_, int range_, int rank_>
void
Function<dim_, codim_, range_, rank_ >::
init_cache(ElementIterator &elem, const topology_variant &sdim) const
{
    init_cache(*elem, k);
}


template<int dim_, int codim_, int range_, int rank_>
void
Function<dim_, codim_, range_, rank_ >::
init_element_cache(ElementAccessor &elem) const
{
    this->init_cache(elem, Topology<dim_>());
}


template<int dim_, int codim_, int range_, int rank_>
void
Function<dim_, codim_, range_, rank_ >::
init_element_cache(ElementIterator &elem) const
{
    this->init_cache(*elem, Topology<dim_>());
}
#endif



#if 0
template<int dim_, int codim_, int range_, int rank_>
void
Function<dim_, codim_, range_, rank_ >::
fill_cache(ElementIterator &elem, const topology_variant &sdim, const int j) const
{
    this->fill_cache(*elem, k, j);
}


template<int dim_, int codim_, int range_, int rank_>
void
Function<dim_, codim_, range_, rank_ >::
fill_element_cache(ElementAccessor &elem) const
{
    this->fill_cache(elem, Topology<dim_>(),0);
}


template<int dim_, int codim_, int range_, int rank_>
void
Function<dim_, codim_, range_, rank_ >::
fill_element_cache(ElementIterator &elem) const
{
    this->fill_cache(*elem, Topology<dim_>(),0);
}
#endif

/*
template<int dim_, int codim_, int range_, int rank_>
auto
Function<dim_, codim_, range_, rank_ >::
get_cache(ElementAccessor &elem)
-> std::shared_ptr<typename ElementAccessor::CacheType> &
{
    Assert(elem.all_sub_elems_cache_ != nullptr,ExcNullPtr());
    return elem.all_sub_elems_cache_;
}
//*/

template<int dim_, int codim_, int range_, int rank_>
auto
Function<dim_, codim_, range_, rank_ >::
begin(const PropId &prop) -> ElementIterator
{
    return ElementIterator(this->shared_from_this(),
    grid_->get_element_property(prop).begin(),
    prop);
}



template<int dim_, int codim_, int range_, int rank_>
auto
Function<dim_, codim_, range_, rank_ >::
end(const PropId &prop) -> ElementIterator
{
    return ElementIterator(this->shared_from_this(),
    grid_->get_element_property(prop).end(),
    prop);
}



template<int dim_, int codim_, int range_, int rank_>
template<int dim_, int codim_, int range_, int rank_>
void
Function<dim_, codim_, range_, rank_ >::
print_info(LogStream &out) const
{
    parent_t::print_info(out);
}



#ifdef SERIALIZATION
template<int dim_, int codim_, int range_, int rank_>
auto
Function<dim_, codim_, range_, rank_ >::
create_element(const ListIt &index, const PropId &property) const -> std::shared_ptr<ElementAccessor>
{
    auto elem = std::make_shared<ElementAccessor>(this->shared_from_this(),index,property);
    Assert(elem != nullptr,ExcNullPtr());

    return elem;
}


template<int dim_, int codim_, int range_, int rank_>
Index
Function<dim_, codim_, range_, rank_ >::
get_object_id() const
{
    return object_id_;
}


template<int dim_, int codim_, int range_, int rank_>
const std::string &
Function<dim_, codim_, range_, rank_ >::
get_name() const
{
    return name_;
}

template<int dim_, int codim_, int range_, int rank_>
void
Function<dim_, codim_, range_, rank_ >::
set_name(const std::string &name)
{
    name_ = name;
}


template<class Archive>
void
Function<dim_, codim_, range_, rank_ >::
serialize(Archive &ar, const unsigned int version)
{
    ar &boost::serialization::make_nvp("grid_elem_handler_",
                                       boost::serialization::base_object<GridElementHandler<dim_>>(*this));

    ar &boost::serialization::make_nvp("object_id_",object_id_);
    ar &boost::serialization::make_nvp("name_",name_);

    ar &boost::serialization::make_nvp("flags_",flags_);

    ar &boost::serialization::make_nvp("grid_",grid_);

#ifdef MESH_REFINEMENT
    ar &boost::serialization::make_nvp("function_previous_refinement_",function_previous_refinement_);
#endif // MESH_REFINEMENT
}
#endif // SERIALIZATION


IGA_NAMESPACE_CLOSE

#include <igatools/functions/function.inst>

