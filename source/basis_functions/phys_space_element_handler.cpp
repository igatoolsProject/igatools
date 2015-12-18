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

#include <igatools/basis_functions/phys_space_element_handler.h>
#include <igatools/basis_functions/physical_space_element.h>

#include <functional>

using std::shared_ptr;



IGA_NAMESPACE_OPEN

template<int dim, int codim>
using MapFunc= Function<dim, 0, dim + codim, 1>;


template<int dim_,int range_,int rank_,int codim_>
PhysSpaceElementHandler<dim_,range_,rank_,codim_>::
PhysSpaceElementHandler(std::shared_ptr<const PhysSpace> space)
  :
  base_t(space),
  ref_space_handler_(space->get_reference_basis()->create_cache_handler()),
  phys_domain_handler_(space->get_physical_domain()->create_cache_handler()),
  phys_space_(space)
{}





template<int dim_,int range_,int rank_,int codim_>
auto
PhysSpaceElementHandler<dim_,range_,rank_,codim_>::
print_info(LogStream &out) const -> void
{
  ref_space_handler_->print_info(out);
  //  PFCache::print_info(out);
}


IGA_NAMESPACE_CLOSE

#include <igatools/basis_functions/phys_space_element_handler.inst>
