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

#include <igatools/basis_functions/physical_space_element.h>
#include <igatools/base/exceptions.h>


using std::shared_ptr;
using std::make_shared;
using std::make_unique;

IGA_NAMESPACE_OPEN

template<int dim_,int range_,int rank_,int codim_,Transformation type_>
PhysicalSpaceElement<dim_,range_,rank_,codim_,type_>::
PhysicalSpaceElement(const std::shared_ptr<ContainerType> phys_space,
                     const ListIt &index,
                     const PropId &prop)
  :
  parent_t(phys_space,index,prop),
  ref_space_element_(phys_space->get_reference_space()->create_ref_element(index,prop)),
  phys_domain_element_(make_unique<PhysDomainElem>(
                        std::const_pointer_cast<PhysDomain>(phys_space->get_physical_domain()),
                        index, prop))
{
//    push_fwd_element_ = std::make_shared<PfElemAccessor>(phys_space->get_map_func(), index);
  Assert(phys_domain_element_ != nullptr, ExcNullPtr());
}


#if 0
template<int dim_,int range_,int rank_,int codim_,Transformation type_>
PhysicalSpaceElement<dim_,range_,rank_,codim_,type_>::
PhysicalSpaceElement(const PhysicalSpaceElement<dim_,range_,rank_,codim_,type_> &in,
                     const CopyPolicy &copy_policy)
  :
  parent_t(in,copy_policy)
{
  if (copy_policy == CopyPolicy::shallow)
  {
    ref_space_element_ = in.ref_space_element_;
    phys_domain_element_ = in.phys_domain_element_;
  }
  else
  {
    ref_space_element_ = std::dynamic_pointer_cast<RefElemAccessor>(in.ref_space_element_->clone());
    phys_domain_element_ = make_shared<PhysDomainElem>(*in.phys_domain_element_);
  }

  Assert(false,ExcNotTested());
}
#endif



template<int dim_,int range_,int rank_,int codim_,Transformation type_>
void
PhysicalSpaceElement<dim_,range_,rank_,codim_,type_>::
operator++()
{
//  parent_t::operator++();

  ++(*phys_domain_element_);

  ++(*ref_space_element_);
}




template<int dim_,int range_,int rank_,int codim_,Transformation type_>
template <int sdim>
auto
PhysicalSpaceElement<dim_,range_,rank_,codim_,type_>::
get_points(const int s_id) const -> const ValueVector<PhysPoint>
{
//  using _Point = typename PhysDomainElem::_Point;
  return phys_domain_element_->template get_points<sdim>(s_id);
}

template<int dim_,int range_,int rank_,int codim_,Transformation type_>
auto
PhysicalSpaceElement<dim_,range_,rank_,codim_,type_>::
get_element_points() const -> const ValueVector<PhysPoint>
{
  return this->template get_points<dim>(0);
}


template<int dim_,int range_,int rank_,int codim_,Transformation type_>
auto
PhysicalSpaceElement<dim_,range_,rank_,codim_,type_>::
get_element_w_measures() const -> const ValueVector<Real>
{
  return this->template get_w_measures<dim>(0);
}


template<int dim_,int range_,int rank_,int codim_,Transformation type_>
auto
PhysicalSpaceElement<dim_,range_,rank_,codim_,type_>::
get_index() const -> IndexType
{
  return parent_t::get_index();
}

#if 0
template<int dim_,int range_,int rank_,int codim_,Transformation type_>
void
PhysicalSpaceElement<dim_,range_,rank_,codim_,type_>::
move_to(const Index flat_index)
{
  this->get_grid_element().move_to(flat_index);
  ref_space_element_->move_to(flat_index);
  phys_domain_element_->move_to(flat_index);
}
#endif



template<int dim_,int range_,int rank_,int codim_,Transformation type_>
auto
PhysicalSpaceElement<dim_,range_,rank_,codim_,type_>::
get_ref_space_element() const -> const RefElemAccessor &
{
  return dynamic_cast<const RefElemAccessor &>(*ref_space_element_);
}

template<int dim_,int range_,int rank_,int codim_,Transformation type_>
auto
PhysicalSpaceElement<dim_,range_,rank_,codim_,type_>::
get_ref_space_element() -> RefElemAccessor &
{
  return dynamic_cast<RefElemAccessor &>(*ref_space_element_);
}

template<int dim_,int range_,int rank_,int codim_,Transformation type_>
auto
PhysicalSpaceElement<dim_,range_,rank_,codim_,type_>::
get_grid() const -> const std::shared_ptr<const Grid<dim> >
{
  return this->get_ref_space_element().get_grid();
}


template<int dim_,int range_,int rank_,int codim_,Transformation type_>
auto
PhysicalSpaceElement<dim_,range_,rank_,codim_,type_>::
get_physical_domain_element() const -> const PhysDomainElem &
{
  return *phys_domain_element_;
}

template<int dim_,int range_,int rank_,int codim_,Transformation type_>
auto
PhysicalSpaceElement<dim_,range_,rank_,codim_,type_>::
get_physical_domain_element() -> PhysDomainElem &
{
  return *phys_domain_element_;
}
//*/


template<int dim_,int range_,int rank_,int codim_,Transformation type_>
void
PhysicalSpaceElement<dim_,range_,rank_,codim_,type_>::
print_info(LogStream &out) const
{
  using std::to_string;
  out.begin_item("ReferenceElement<" +
                 to_string(RefSpace::dim) + "," +
                 to_string(RefSpace::range) + "," +
                 to_string(RefSpace::rank) + ">");
  ref_space_element_->print_info(out);
  out.end_item();

  out.begin_item("DomainElement<" +
                 to_string(dim_) + "," +
                 to_string(codim_) + ">");
  phys_domain_element_->print_info(out);
  out.end_item();
}

template<int dim_,int range_,int rank_,int codim_,Transformation type_>
void
PhysicalSpaceElement<dim_,range_,rank_,codim_,type_>::
print_cache_info(LogStream &out) const
{
  using std::to_string;
  out.begin_item("ReferenceElement<" +
                 to_string(RefSpace::dim) + "," +
                 to_string(RefSpace::range) + "," +
                 to_string(RefSpace::rank) + "> cache:");
  ref_space_element_->print_cache_info(out);
  out.end_item();

  out.begin_item("DomainElement<" +
                 to_string(dim_) + "," +
                 to_string(codim_) + "> cache:");
  phys_domain_element_->print_cache_info(out);
  out.end_item();
}



template<int dim_,int range_,int rank_,int codim_,Transformation type_>
bool
PhysicalSpaceElement<dim_,range_,rank_,codim_,type_>::
operator==(const parent_t &a) const
{
  return !((*this) != a);
}

template<int dim_,int range_,int rank_,int codim_,Transformation type_>
bool
PhysicalSpaceElement<dim_,range_,rank_,codim_,type_>::
operator!=(const parent_t &a) const
{
  const self_t &elem = dynamic_cast<const self_t &>(a);
  Assert(this->is_comparable_with(elem),
         ExcMessage("The elements are not comparable"));
  return (*ref_space_element_ != *elem.ref_space_element_) ||
         (*phys_domain_element_ != *elem.phys_domain_element_);
}

template<int dim_,int range_,int rank_,int codim_,Transformation type_>
bool
PhysicalSpaceElement<dim_,range_,rank_,codim_,type_>::
operator<(const parent_t &a) const
{
  const self_t &elem = dynamic_cast<const self_t &>(a);
  Assert(this->is_comparable_with(elem),
         ExcMessage("The elements are not comparable"));
  return (*ref_space_element_ < *elem.ref_space_element_);
}

template<int dim_,int range_,int rank_,int codim_,Transformation type_>
bool
PhysicalSpaceElement<dim_,range_,rank_,codim_,type_>::
operator>(const parent_t &a) const
{
  const self_t &elem = dynamic_cast<const self_t &>(a);
  Assert(this->is_comparable_with(elem),
         ExcMessage("The elements are not comparable"));
  return (*ref_space_element_ > *elem.ref_space_element_);
}


template<int dim_,int range_,int rank_,int codim_,Transformation type_>
bool
PhysicalSpaceElement<dim_,range_,rank_,codim_,type_>::
is_comparable_with(const self_t &elem) const
{
//   bool same_ref_space =
//       (ref_space_element_->get_space() == elem.ref_space_element_->get_space());
  const bool same_phys_space =
    (this->get_space() == elem.get_space());
  return (same_phys_space);
}

IGA_NAMESPACE_CLOSE

#include <igatools/basis_functions/physical_space_element.inst>
