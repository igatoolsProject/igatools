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
Index
Function<dim_, codim_, range_, rank_ >::
get_object_id() const
{
    return object_id_;
}

template<int dim_, int codim_, int range_, int rank_>
void
Function<dim_, codim_, range_, rank_ >::
reset(const ValueFlags &flag, const eval_pts_variant &quad)
{
    auto reset_dispatcher = ResetDispatcher(flag,*this,flags_);
    boost::apply_visitor(reset_dispatcher, quad);
}

template<int dim_, int codim_, int range_, int rank_>
void
Function<dim_, codim_, range_, rank_ >::
init_cache(ElementAccessor &func_elem, const topology_variant &k) const
{
    auto init_cache_dispatcher = InitCacheDispatcher(*this,flags_,func_elem);

    boost::apply_visitor(init_cache_dispatcher, k);
}

template<int dim_, int codim_, int range_, int rank_>
void
Function<dim_, codim_, range_, rank_ >::
init_cache(ElementIterator &elem, const topology_variant &k) const
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

template<int dim_, int codim_, int range_, int rank_>
void
Function<dim_, codim_, range_, rank_ >::
fill_cache(ElementAccessor &func_elem, const topology_variant &k,const int sub_elem_id) const
{
    auto fill_cache_dispatcher = FillCacheDispatcher(sub_elem_id,*this,func_elem);
    boost::apply_visitor(fill_cache_dispatcher, k);
}

template<int dim_, int codim_, int range_, int rank_>
void
Function<dim_, codim_, range_, rank_ >::
fill_cache(ElementIterator &elem, const topology_variant &k, const int j) const
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
create_element(const Index flat_index) const -> std::shared_ptr<ElementAccessor>
{
    auto elem = std::make_shared<ElementAccessor>(this->shared_from_this(),flat_index);
    Assert(elem != nullptr,ExcNullPtr());

    return elem;
}


template<int dim_, int codim_, int range_, int rank_>
void
Function<dim_, codim_, range_, rank_ >::
print_info(LogStream &out) const
{
    parent_t::print_info(out);
}


template<int dim_, int codim_, int range_, int rank_>
auto
Function<dim_, codim_, range_, rank_ >::
begin() const -> ElementIterator
{
    return ElementIterator(this->create_element(0),ElementProperties::active);
}

template<int dim_, int codim_, int range_, int rank_>
auto
Function<dim_, codim_, range_, rank_ >::
end() const -> ElementIterator
{
    return ElementIterator(this->create_element(IteratorState::pass_the_end),ElementProperties::active);
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


#ifdef SERIALIZATION
template<int dim_, int codim_, int range_, int rank_>
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
    ar &boost::serialization::make_nvp("function_previous_refinement_",function_previous_refinement_);
}
#endif // SERIALIZATION


IGA_NAMESPACE_CLOSE

#include <igatools/functions/function.inst>

