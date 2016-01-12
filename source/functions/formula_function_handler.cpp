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

#include <igatools/functions/formula_function_handler.h>

IGA_NAMESPACE_OPEN

template<int dim_, int codim_, int range_, int rank_ >
FormulaFunctionHandler<dim_, codim_, range_, rank_ >::
FormulaFunctionHandler(std::shared_ptr<FuncType> func)
  :
  parent_t::FunctionHandler(func)
{}



template<int dim_, int codim_, int range_, int rank_ >
void
FormulaFunctionHandler<dim_, codim_, range_, rank_ >::
set_flags(const topology_variant &sdim, const Flags &flag)

{
  parent_t::set_flags(sdim, flag);
  this->get_domain_handler()->set_flags(sdim, DomainHandlerType::Flags::point);
}

template<int dim_, int codim_, int range_, int rank_>
void
FormulaFunctionHandler<dim_, codim_, range_, rank_ >::
fill_cache(const topology_variant &sdim,
           ElementAccessor &elem,
           const int s_id) const
{
  this->get_domain_handler()->fill_cache(sdim, elem.get_domain_element(), s_id);

  parent_t::fill_cache(sdim, elem, s_id);

  FuncType &func = *std::dynamic_pointer_cast<FuncType>(this->get_function());
  auto disp = FillCacheDispatcher(func, *this, elem, s_id);
  boost::apply_visitor(disp, sdim);

//  parent_t::fill_cache(sdim, elem, s_id);

}

IGA_NAMESPACE_CLOSE

#include <igatools/functions/formula_function_handler.inst>

