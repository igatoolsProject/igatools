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

#include <igatools/functions/formula_function.h>
#include <igatools/functions/function_element.h>
#include <igatools/geometry/physical_domain_element.h>

using std::shared_ptr;

IGA_NAMESPACE_OPEN

template<int dim, int codim, int range, int rank>
FormulaFunction<dim, codim, range, rank>::
FormulaFunction(std::shared_ptr<GridType> grid, std::shared_ptr<PhysDomain> map)
    :
    parent_t::Function(grid),
    mapping_(PhysDom::create(map)),
    map_elem_(mapping_->begin())
{}



template<int dim, int codim, int range, int rank>
FormulaFunction<dim, codim, range, rank>::
FormulaFunction(const self_t &func)
    :
    parent_t::Function(func),
    mapping_(func.mapping_),
    map_elem_(func.mapping_->begin())
{}



template<int dim, int codim, int range, int rank>
void
FormulaFunction<dim, codim, range, rank>::
reset(const ValueFlags &flag, const eval_pts_variant &quad)
{
    parent_t::reset(flag, quad);
    mapping_->reset(ValueFlags::value|ValueFlags::point, quad);
}



template<int dim, int codim, int range, int rank>
void
FormulaFunction<dim, codim, range, rank>::
init_cache(ElementAccessor &elem, const topology_variant &k) const
{
    parent_t::init_cache(elem, k);
    using MapElem = typename PhysDomain::ElementAccessor;
    mapping_->init_cache(const_cast<MapElem &>(*map_elem_), k);
}



template<int dim, int codim, int range, int rank>
auto
FormulaFunction<dim, codim, range, rank>::
fill_cache(ElementAccessor &elem, const topology_variant &k, const int sub_elem_id) const  -> void
{
    parent_t::fill_cache(elem,k,sub_elem_id);
    using MapElem = typename PhysDomain::ElementAccessor;

    auto & map_elem_non_const = const_cast<MapElem &>(*map_elem_);
    map_elem_non_const.move_to(elem.get_flat_index());
    mapping_->fill_cache(map_elem_non_const,k,sub_elem_id);

    auto fill_cache_dispatcher = FillCacheDispatcher(sub_elem_id,*this,elem);
    boost::apply_visitor(fill_cache_dispatcher, k);
}

IGA_NAMESPACE_CLOSE

#include <igatools/functions/formula_function.inst>
