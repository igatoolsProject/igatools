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

#include <igatools/base/new_function.h>
#include <igatools/base/function_element.h>
IGA_NAMESPACE_OPEN


template<int dim, int range, int rank >
FormulaFunction<dim, range, rank >::
FormulaFunction(std::shared_ptr<const CartesianGrid<dim>> grid,
		const ValueFlags flag,
		const Quadrature<dim> &quad)
		 :
		parent_t::NewFunction(grid, flag, quad),
		flag_(flag),
		quad_(quad)
{}



template<int dim, int range, int rank >
auto
FormulaFunction<dim, range, rank >::
init_element(ElementIterator &elem) -> void
{
	auto &el = elem.get_accessor();
	GridUniformQuadCache<dim>::init_element_cache(el);
	auto &cache = this->get_cache(elem);
	if (cache == nullptr)
	{
		using Cache = typename ElementAccessor::CacheType;
		cache = shared_ptr<Cache>(new Cache);
	}
	cache->resize(flag_, quad_.get_num_points());
}



template<int dim, int range, int rank >
auto
FormulaFunction<dim, range, rank >::
fill_element(ElementIterator &elem) -> void
{
	auto &el    = elem.get_accessor();
	GridUniformQuadCache<dim>::fill_element_cache(el);
	const auto points = el.get_points();
	auto &cache = this->get_cache(elem);
	if (flag_.fill_values())
		this->evaluate_0(points, cache->values_);
	if (flag_.fill_gradients())
		this->evaluate_1(points, std::get<1>(cache->derivatives_));
	if (flag_.fill_hessians())
		this->evaluate_2(points, std::get<2>(cache->derivatives_));
}


template<int dim, int range, int rank >
auto
FormulaFunction<dim, range, rank >::
->

template<int dim, int range, int rank>
auto
NewFunction<dim,range,rank>::
get_cache(NewFunction<dim,range,rank>::ElementIterator &elem)
-> std::shared_ptr<typename ElementAccessor::CacheType> &
{
    return elem.get_accessor().elem_cache_;
}



IGA_NAMESPACE_CLOSE

#include <igatools/base/new_function.inst>

