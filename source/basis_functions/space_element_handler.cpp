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

#include <igatools/basis_functions/space_element_handler.h>
#include <igatools/basis_functions/space_element.h>


using std::shared_ptr;


IGA_NAMESPACE_OPEN

template<int dim,int codim,int range,int rank>
SpaceElementHandler<dim,codim,range,rank>::
SpaceElementHandler(const std::shared_ptr<const Sp> &space)
  :
  space_(space)
{
  Assert(space != nullptr,ExcNullPtr());
}




template<int dim,int codim,int range,int rank>
auto
SpaceElementHandler<dim,codim,range,rank>::
get_space() const -> std::shared_ptr<const Sp>
{
  return space_;
}


template<int dim,int codim,int range,int rank>
void
SpaceElementHandler<dim,codim,range,rank>::
set_element_flags(const typename space_element::Flags &flag)
{
  this->set_flags_impl(Topology<dim>(),flag);
}


template<int dim,int codim,int range,int rank>
void
SpaceElementHandler<dim,codim,range,rank>::
init_element_cache(ElementAccessor &elem,
                   std::shared_ptr<const Quadrature<dim>> quad) const
{
  init_cache<dim>(elem, quad);
}

template<int dim,int codim,int range,int rank>
void
SpaceElementHandler<dim,codim,range,rank>::
init_element_cache(ElementIterator &elem,
                   std::shared_ptr<const Quadrature<dim>> quad) const
{
  init_element_cache(*elem, quad);
}

template<int dim,int codim,int range,int rank>
void
SpaceElementHandler<dim,codim,range,rank>::
init_face_cache(ElementAccessor &elem,
                std::shared_ptr<const Quadrature<(dim > 0) ? dim-1 : 0>> quad) const
{
  Assert(dim > 0,ExcMessage("No face defined for element with topological dimension 0."));
  init_cache<(dim > 0) ? dim-1 : 0>(elem, quad);
}

template<int dim,int codim,int range,int rank>
void
SpaceElementHandler<dim,codim,range,rank>::
init_face_cache(ElementIterator &elem,
                std::shared_ptr<const Quadrature<(dim > 0) ? dim-1 : 0>> quad) const
{
  init_face_cache(*elem, quad);
}


IGA_NAMESPACE_CLOSE

#include <igatools/basis_functions/space_element_handler.inst>
