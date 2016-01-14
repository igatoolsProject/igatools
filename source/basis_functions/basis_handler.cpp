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

#include <igatools/basis_functions/basis_handler.h>
#include <igatools/basis_functions/basis_element.h>


using std::shared_ptr;


IGA_NAMESPACE_OPEN

template<int dim,int codim,int range,int rank>
BasisHandler<dim,codim,range,rank>::
BasisHandler(const std::shared_ptr<const Bs> &basis)
  :
  basis_(basis)
{
  Assert(basis != nullptr,ExcNullPtr());
}




template<int dim,int codim,int range,int rank>
auto
BasisHandler<dim,codim,range,rank>::
get_basis() const -> std::shared_ptr<const Bs>
{
  return basis_;
}


template<int dim,int codim,int range,int rank>
template<int sdim>
void
BasisHandler<dim,codim,range,rank>::
set_flags(const typename basis_element::Flags &flag)
{
  this->set_flags_impl(Topology<sdim>(),flag);
}

template<int dim,int codim,int range,int rank>
void
BasisHandler<dim,codim,range,rank>::
set_element_flags(const typename basis_element::Flags &flag)
{
  this->set_flags_impl(Topology<dim>(),flag);
}


template<int dim,int codim,int range,int rank>
template<int sdim>
void
BasisHandler<dim,codim,range,rank>::
init_cache(ElementAccessor &elem,
           const std::shared_ptr<const Quadrature<sdim>> &quad) const
{
  this->init_cache_impl(elem,quad);
}

template<int dim,int codim,int range,int rank>
template<int sdim>
void
BasisHandler<dim,codim,range,rank>::
init_cache(ElementIterator &elem,
           const std::shared_ptr<const Quadrature<sdim>> &quad) const
{
  init_cache<sdim>(*elem, quad);
}


template<int dim,int codim,int range,int rank>
void
BasisHandler<dim,codim,range,rank>::
init_element_cache(ElementAccessor &elem,
                   const std::shared_ptr<const Quadrature<dim>> &quad) const
{
  init_cache<dim>(elem, quad);
}

template<int dim,int codim,int range,int rank>
void
BasisHandler<dim,codim,range,rank>::
init_element_cache(ElementIterator &elem,
                   const std::shared_ptr<const Quadrature<dim>> &quad) const
{
  init_element_cache(*elem, quad);
}

template<int dim,int codim,int range,int rank>
void
BasisHandler<dim,codim,range,rank>::
init_face_cache(ElementAccessor &elem,
                const std::shared_ptr<const Quadrature<(dim > 0) ? dim-1 : 0>> &quad) const
{
  Assert(dim > 0,ExcMessage("No face defined for element with topological dimension 0."));
  init_cache<(dim > 0) ? dim-1 : 0>(elem, quad);
}

template<int dim,int codim,int range,int rank>
void
BasisHandler<dim,codim,range,rank>::
init_face_cache(ElementIterator &elem,
                const std::shared_ptr<const Quadrature<(dim > 0) ? dim-1 : 0>> &quad) const
{
  init_face_cache(*elem, quad);
}


template<int dim,int codim,int range,int rank>
template<int sdim>
void
BasisHandler<dim,codim,range,rank>::
fill_cache(ElementAccessor &elem, const int s_id) const
{
  Assert(s_id >= 0 && s_id < UnitElement<dim>::template num_elem<sdim>(),
         ExcIndexRange(s_id,0,UnitElement<dim>::template num_elem<sdim>()));
  this->fill_cache_impl(Topology<sdim>(),elem,s_id);
}


template<int dim,int codim,int range,int rank>
template<int sdim>
void
BasisHandler<dim,codim,range,rank>::
fill_cache(ElementIterator &elem, const int s_id) const
{
  fill_cache<sdim>(*elem, s_id);
}


template<int dim,int codim,int range,int rank>
void
BasisHandler<dim,codim,range,rank>::
fill_element_cache(ElementAccessor &elem) const
{
  fill_cache<dim>(elem,0);
}

template<int dim,int codim,int range,int rank>
void
BasisHandler<dim,codim,range,rank>::
fill_element_cache(ElementIterator &elem) const
{
  fill_element_cache(*elem);
}


template<int dim,int codim,int range,int rank>
void
BasisHandler<dim,codim,range,rank>::
fill_face_cache(ElementAccessor &elem, const int s_id) const
{
  Assert(dim > 0,ExcMessage("No face defined for element with topological dimension 0."));
  fill_cache<(dim > 0) ? dim-1 : 0>(elem,s_id);
}

template<int dim,int codim,int range,int rank>
void
BasisHandler<dim,codim,range,rank>::
fill_face_cache(ElementIterator &elem, const int s_id) const
{
  fill_face_cache(*elem,s_id);
}


template<int dim,int codim,int range,int rank>
auto
BasisHandler<dim,codim,range,rank>::
get_element_cache(ElementAccessor &elem) const -> typename ElementAccessor::CacheType &
{
  return elem.all_sub_elems_cache_;
}




IGA_NAMESPACE_CLOSE

#include <igatools/basis_functions/basis_handler.inst>
