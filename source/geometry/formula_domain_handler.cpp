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

#include <igatools/geometry/formula_domain_handler.h>

IGA_NAMESPACE_OPEN

template<int dim, int codim>
FormulaDomainHandler<dim, codim>::
FormulaDomainHandler(std::shared_ptr<DomainType> domain)
  :
  parent_t::PhysicalDomainElementHandler(domain)
{}



template<int dim, int codim>
auto
FormulaDomainHandler<dim, codim>::
fill_cache(const topology_variant &sdim,
           ConstElementAccessor &elem,
           const int s_id) const  -> void
{
  parent_t::fill_cache(sdim, elem, s_id);

// std::dynamic_pointer_cast<DomainType>(this->get_domain());
  auto disp = FillCacheDispatcher(*std::dynamic_pointer_cast<DomainType>(this->get_domain()), *this, elem, s_id);
  boost::apply_visitor(disp, sdim);
}

IGA_NAMESPACE_CLOSE

#include <igatools/geometry/formula_domain_handler.inst>

