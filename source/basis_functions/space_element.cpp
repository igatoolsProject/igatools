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



#include <igatools/basis_functions/space_element.h>
#include <igatools/basis_functions/reference_element.h>
#include <igatools/basis_functions/physical_space_element.h>


IGA_NAMESPACE_OPEN


template<int dim_,int codim_,int range_,int rank_,Transformation type_>
SpaceElement<dim_,codim_,range_,rank_,type_>::
SpaceElement(const std::shared_ptr<const Space<dim_,codim_,range_,rank_,type_>> space,
             const ListIt &index,
             const PropId &prop)
  :
  base_t(space,index,prop),
  space_(space)
{}






template<int dim_,int codim_,int range_,int rank_,Transformation type_>
void
SpaceElement<dim_,codim_,range_,rank_,type_>::
print_info(LogStream &out) const
{
  base_t::print_info(out);
}

template<int dim_,int codim_,int range_,int rank_,Transformation type_>
void
SpaceElement<dim_,codim_,range_,rank_,type_>::
print_cache_info(LogStream &out) const
{
  out.begin_item("SpaceElementBase<" + std::to_string(dim_) + "> cache:");
  base_t::print_cache_info(out);
  out.end_item();

//    Assert(all_sub_elems_cache_ != nullptr,ExcNullPtr());
  if (all_sub_elems_cache_)
    all_sub_elems_cache_->print_info(out);
  else
    out << "Cache not allocated." << std::endl;
}


template<int dim_,int codim_,int range_,int rank_,Transformation type_>
template <int k>
ValueVector<Real>
SpaceElement<dim_,codim_,range_,rank_,type_>::
get_w_measures(const int j) const
{
  ValueVector<Real> w_measures;

  using RefElem = const ReferenceElement<dim_,range_,rank_>;
  RefElem *as_ref_elem = dynamic_cast<RefElem *>(this);
  if (as_ref_elem)
    w_measures = as_ref_elem->template get_w_measures<k>(j);

  using PhysElem = const PhysicalSpaceElement<dim_,range_,rank_,codim_>;
  PhysElem *as_phys_elem = dynamic_cast<PhysElem *>(this);
  if (as_phys_elem)
    w_measures = as_phys_elem->template get_w_measures<k>(j);

  return w_measures;
}

template<int dim_,int codim_,int range_,int rank_,Transformation type_>
auto
SpaceElement<dim_,codim_,range_,rank_,type_>::
get_space() const -> std::shared_ptr<const Sp>
{
  return space_;
}





#ifdef SERIALIZATION
template<int dim_,int codim_,int range_,int rank_,Transformation type_>
template<class Archive>
void
SpaceElement<dim_,codim_,range_,rank_,type_>::
serialize(Archive &ar, const unsigned int version)
{
  ar &boost::serialization::make_nvp("SpaceElement_base_t",
                                     boost::serialization::base_object<SpaceElementBase<dim_>>(*this));

  ar &boost::serialization::make_nvp("all_sub_elems_cache_",all_sub_elems_cache_);

  auto non_const_space = std::const_pointer_cast<Space<dim_,codim_,range_,rank_,type_>>(space_);
  ar &boost::serialization::make_nvp("space_",non_const_space);
  space_ = non_const_space;
  Assert(space_ != nullptr,ExcNullPtr());
}
#endif // SERIALIZATION



IGA_NAMESPACE_CLOSE


#include <igatools/basis_functions/space_element.inst>


