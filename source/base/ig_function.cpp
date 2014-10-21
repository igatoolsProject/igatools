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

#include <igatools/base/ig_function.h>
#include <igatools/base/function_element.h>

using std::shared_ptr;

IGA_NAMESPACE_OPEN

template<class Space>
IgFunction<Space>::
IgFunction(const NewValueFlags &flag, const Quadrature<dim> &quad,
           std::shared_ptr<const Space> space,
           const CoeffType &coeff)
    :
    parent_t::NewFunction(space->get_grid(), flag, quad),
    flag_(flag),
    quad_(quad),
    space_(space),
    coeff_(coeff),
    elem_(space_->begin()),
    space_filler_(space_, flag, quad)
{}



template<class Space>
auto
IgFunction<Space>::
init_elem(ElementAccessor &elem) -> void
{
    GridElementHandler<dim>::init_element_cache(elem);
    auto &cache = this->get_cache(elem);
    if (cache == nullptr)
    {
        using Cache = typename ElementAccessor::CacheType;
        cache = shared_ptr<Cache>(new Cache);
    }
    cache->resize(flag_, quad_.get_num_points());

    space_filler_.init_element_cache(elem_);
}



template<class Space>
auto
IgFunction<Space>::
fill_elem(ElementAccessor &elem) -> void
{
    GridElementHandler<dim>::fill_element_cache(elem);
    space_filler_.fill_element_cache(elem_);
    auto &cache = this->get_cache(elem);

    elem_.move_to(elem.get_flat_index());
    if (flag_.fill_points())
        Assert(false, ExcNotImplemented());//cache->points_ = elem_->get_points();
    const auto loc_coeff = coeff_.get_local_coefs(elem_->get_local_to_global());

    if (flag_.fill_values())
        cache->values_ = elem_->template eval_field_ders<0, 0>(0, loc_coeff);
   if (flag_.fill_gradients())
        std::get<1>(cache->derivatives_) =
                elem_->template eval_field_ders<0, 1>(0, loc_coeff);
    if (flag_.fill_hessians())
        std::get<2>(cache->derivatives_) =
                elem_->template eval_field_ders<0, 2>(0, loc_coeff);
}

IGA_NAMESPACE_CLOSE

#include <igatools/base/ig_function.inst>
