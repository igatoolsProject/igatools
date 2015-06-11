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

IGA_NAMESPACE_OPEN

template<int dim_, int codim_, int range_, int rank_ >
Function<dim_, codim_, range_, rank_ >::
Function(std::shared_ptr<GridType> grid)
    :
    GridElementHandler<dim_>(grid)
#ifdef MESH_REFINEMENT
    ,
    functions_knots_refinement_(std::const_pointer_cast<CartesianGrid<dim_>>(grid))
#endif
{}


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
init_cache(ElementAccessor &func_elem, const topology_variant &k)
{
    auto init_cache_dispatcher = InitCacheDispatcher(*this,flags_,func_elem);

    boost::apply_visitor(init_cache_dispatcher, k);
}

template<int dim_, int codim_, int range_, int rank_>
void
Function<dim_, codim_, range_, rank_ >::
init_cache(ElementIterator &elem, const topology_variant &k)
{
    init_cache(*elem, k);
}


template<int dim_, int codim_, int range_, int rank_>
void
Function<dim_, codim_, range_, rank_ >::
init_element_cache(ElementAccessor &elem)
{
    this->init_cache(elem, Topology<dim_>());
}


template<int dim_, int codim_, int range_, int rank_>
void
Function<dim_, codim_, range_, rank_ >::
init_element_cache(ElementIterator &elem)
{
    this->init_cache(*elem, Topology<dim_>());
}

template<int dim_, int codim_, int range_, int rank_>
void
Function<dim_, codim_, range_, rank_ >::
fill_cache(ElementAccessor &func_elem, const topology_variant &k,const int sub_elem_id)
{
    auto fill_cache_dispatcher = FillCacheDispatcher(sub_elem_id,*this,func_elem);
    boost::apply_visitor(fill_cache_dispatcher, k);
}

template<int dim_, int codim_, int range_, int rank_>
void
Function<dim_, codim_, range_, rank_ >::
fill_cache(ElementIterator &elem, const topology_variant &k, const int j)
{
    fill_cache(*elem, k, j);
}


template<int dim_, int codim_, int range_, int rank_>
void
Function<dim_, codim_, range_, rank_ >::
fill_element_cache(ElementAccessor &elem)
{
    this->fill_cache(elem, Topology<dim_>(),0);
}


template<int dim_, int codim_, int range_, int rank_>
void
Function<dim_, codim_, range_, rank_ >::
fill_element_cache(ElementIterator &elem)
{
    this->fill_cache(*elem, Topology<dim_>(),0);
}

template<int dim_, int codim_, int range_, int rank_>
auto
Function<dim_, codim_, range_, rank_ >::
get_cache(ElementAccessor &elem)
-> std::shared_ptr<typename ElementAccessor::CacheType> &
{
    Assert(elem.all_sub_elems_cache_ != nullptr,ExcNullPtr());
    return elem.all_sub_elems_cache_;
}

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
    return ElementIterator(this->create_element(0),ElementProperties::none);
}

template<int dim_, int codim_, int range_, int rank_>
auto
Function<dim_, codim_, range_, rank_ >::
end() const -> ElementIterator
{
    return ElementIterator(this->create_element(IteratorState::pass_the_end),ElementProperties::none);
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

    ar &boost::serialization::make_nvp("flags_",flags_);

#ifdef MESH_REFINEMENT
    ar &boost::serialization::make_nvp("functions_knots_refinement_",functions_knots_refinement_);
#endif // MESH_REFINEMENT
}
#endif // SERIALIZATION


IGA_NAMESPACE_CLOSE

#include <igatools/functions/function.inst>

