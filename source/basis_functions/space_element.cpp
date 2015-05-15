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


template<int dim_,int codim_,int range_,int rank_>
SpaceElement<dim_,codim_,range_,rank_>::
SpaceElement(const std::shared_ptr<const Space<dim_,codim_,range_,rank_>> space,
             const Index elem_index)
    :
    base_t(space,elem_index),
    space_(space)
{}



template<int dim_,int codim_,int range_,int rank_>
SpaceElement<dim_,codim_,range_,rank_>::
SpaceElement(const self_t &elem,
             const CopyPolicy &copy_policy)
    :
    base_t(elem,copy_policy)
{
    if (elem.all_sub_elems_cache_ != nullptr)
    {
        if (copy_policy == CopyPolicy::shallow)
        {
            all_sub_elems_cache_ = elem.all_sub_elems_cache_;
        }
        else
        {
            all_sub_elems_cache_ = std::make_shared<LocalCache<Cache>>(*elem.all_sub_elems_cache_);
        }
    }
}



template<int dim_,int codim_,int range_,int rank_>
void
SpaceElement<dim_,codim_,range_,rank_>::
copy_from(const self_t &elem,
          const CopyPolicy &copy_policy)
{
    if (this != &elem)
    {
        base_t::copy_from(elem,copy_policy);

        space_ = elem.space_;
        if (copy_policy == CopyPolicy::deep)
        {
            Assert(elem.all_sub_elems_cache_ != nullptr, ExcNullPtr());
            all_sub_elems_cache_ = std::make_shared<LocalCache<Cache>>(*elem.all_sub_elems_cache_);
        }
        else if (copy_policy == CopyPolicy::shallow)
        {
            all_sub_elems_cache_ = elem.all_sub_elems_cache_;
        }
        else
        {
            Assert(false,ExcNotImplemented());
            AssertThrow(false,ExcNotImplemented());
        }
    }
}



template<int dim_,int codim_,int range_,int rank_>
void
SpaceElement<dim_,codim_,range_,rank_>::
deep_copy_from(const self_t &elem)
{
    this->copy_from(elem,CopyPolicy::deep);
}



template<int dim_,int codim_,int range_,int rank_>
void
SpaceElement<dim_,codim_,range_,rank_>::
shallow_copy_from(const self_t &elem)
{
    this->copy_from(elem,CopyPolicy::shallow);
}

template<int dim_,int codim_,int range_,int rank_>
auto
SpaceElement<dim_,codim_,range_,rank_>::
clone() const -> std::shared_ptr<self_t>
{
    Assert(false,ExcMessage("This function must not be called. "
    "You should call the clone() function of a derived base class."));
    return nullptr;
}


template<int dim_,int codim_,int range_,int rank_>
auto
SpaceElement<dim_,codim_,range_,rank_>::
operator=(const self_t &element) -> self_t &
{
    this->shallow_copy_from(element);
    return (*this);
}



template<int dim_,int codim_,int range_,int rank_>
void
SpaceElement<dim_,codim_,range_,rank_>::
print_info(LogStream &out) const
{
    base_t::print_info(out);
}

template<int dim_,int codim_,int range_,int rank_>
void
SpaceElement<dim_,codim_,range_,rank_>::
print_cache_info(LogStream &out) const
{
    base_t::print_cache_info(out);

    Assert(all_sub_elems_cache_ != nullptr,ExcNullPtr());
    all_sub_elems_cache_->print_info(out);
}


template<int dim_,int codim_,int range_,int rank_>
template <int k>
ValueVector<Real>
SpaceElement<dim_,codim_,range_,rank_>::
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



#ifdef SERIALIZATION
template<int dim_,int codim_,int range_,int rank_>
template<class Archive>
void
SpaceElement<dim_,codim_,range_,rank_>::
serialize(Archive &ar, const unsigned int version)
{
    ar &boost::serialization::make_nvp("SpaceElement_base_t",
                                       boost::serialization::base_object<SpaceElementBase<dim_>>(*this));

    ar &boost::serialization::make_nvp("all_sub_elems_cache_",all_sub_elems_cache_);

    auto non_const_space = std::const_pointer_cast<Space<dim_,codim_,range_,rank_>>(space_);
    ar &boost::serialization::make_nvp("space_",non_const_space);
    space_ = non_const_space;
    Assert(space_ != nullptr,ExcNullPtr());
}
#endif // SERIALIZATION



IGA_NAMESPACE_CLOSE


#include <igatools/basis_functions/space_element.inst>


