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

#include <igatools/base/tensor.h>
#include <igatools/geometry/new_mapping.h>
#include <igatools/geometry/new_mapping_element_accessor.h>

using std::shared_ptr;

IGA_NAMESPACE_OPEN

template<int dim, int codim>
auto
NewMapping<dim, codim>::
get_cache(ElementAccessor &elem)
-> std::shared_ptr<typename ElementAccessor::CacheType> &
{
    return elem.elem_cache_;
}



template<int dim, int codim>
NewMapping<dim, codim>::
NewMapping(std::shared_ptr<FuncType> F,
           const ValueFlags flag,
           const Quadrature<dim> &quad)
    :
    F_(F),
    flag_(flag),
    quad_(quad)
{}



template<int dim, int codim>
NewMapping<dim, codim>::
~NewMapping()
{}



template<int dim, int codim>
auto
NewMapping<dim, codim>::
init_element(ElementIterator &elem) -> void
{
    auto &el = elem.get_accessor();
    F_->init_elem(el);
    auto &cache = this->get_cache(el);
    if (cache == nullptr)
    {
        using Cache = typename ElementAccessor::CacheType;
        cache = shared_ptr<Cache>(new Cache);
    }
    cache->resize(flag_, quad_.get_num_points());
}



template<int dim, int codim>
auto
NewMapping<dim, codim>::
fill_element(ElementIterator &elem) -> void
{
    auto &el    = elem.get_accessor();
    F_->fill_elem(el);
    //const auto points = el.CartesianGridElement<dim>::get_points();
    const auto n_points = quad_.get_num_points();

    auto &cache = this->get_cache(el);
    if (flag_.fill_measures())
    {
        const auto &DF = el.get_gradients();
        for (int i=0; i<n_points; ++i)
            cache->measures_[i] = determinant<dim,space_dim>(DF[i]);

    }
    if (flag_.fill_w_measures())

    {
        const auto &meas = cache->measures_;
        const auto &w = el.CartesianGridElement<dim>::get_w_measures();
        for (int i=0; i<n_points; ++i)
            cache->w_measures_[i] = w[i] * meas[i];
    }
    //    if (flag_.fill_values())
//        this->evaluate_0(cache->points_, cache->values_);
//    if (flag_.fill_gradients())
//        this->evaluate_1(cache->points_, std::get<1>(cache->derivatives_));
//    if (flag_.fill_hessians())
//        this->evaluate_2(cache->points_, std::get<2>(cache->derivatives_));
}

IGA_NAMESPACE_CLOSE

#include <igatools/geometry/new_mapping.inst>

