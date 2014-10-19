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

#include <igatools/base/formula_function.h>
#include <igatools/base/function_element.h>

using std::shared_ptr;

//namespace
//{
//auto function_to_grid_flag(const ValueFlags &fun_flag)
//{
//    ValueFlags grid_flag;
//    if (contain)
//}
//
//}
IGA_NAMESPACE_OPEN

template<int dim, int codim, int range, int rank>
FormulaFunction<dim, codim, range, rank>::
FormulaFunction(shared_ptr<const CartesianGrid<dim>> grid,
                const ValueFlags &flag,
                const Quadrature<dim> &quad)
    :
    parent_t::NewFunction(grid, flag, quad),
    flag_(flag),
    quad_(quad)
{}



template<int dim, int codim, int range, int rank>
void
FormulaFunction<dim, codim, range, rank>::
reset(const ValueFlags &flag, const Quadrature<dim> &quad)
{
	parent_t::reset(flag, quad);
    flag_ = flag;
    quad_ = quad;
}


template<int dim, int codim, int range, int rank>
auto
FormulaFunction<dim, codim, range, rank>::
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
}



template<int dim, int codim, int range, int rank>
auto
FormulaFunction<dim, codim, range, rank>::
fill_elem(ElementAccessor &elem) -> void
{
    GridElementHandler<dim>::fill_element_cache(elem);
    if (!flag_.fill_none())
    {
        const auto points = elem.CartesianGridElement<dim>::get_points();
        auto &cache = this->get_cache(elem);
        if (flag_.fill_points())
            this->parametrization(points, cache->points_);
        if (flag_.fill_values())
            this->evaluate_0(cache->points_, cache->values_);
        if (flag_.fill_gradients())
            this->evaluate_1(cache->points_, std::get<1>(cache->derivatives_));
        if (flag_.fill_hessians())
            this->evaluate_2(cache->points_, std::get<2>(cache->derivatives_));
    }
}

IGA_NAMESPACE_CLOSE

#include <igatools/base/formula_function.inst>
