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

IGA_NAMESPACE_OPEN


template<int dim,int codim,int range,int rank>
SpaceElement<dim,codim,range,rank>::
SpaceElement(const std::shared_ptr<const Space<dim>> space,
             const Index elem_index)
    :
    base_t(space,elem_index)
{}



template<int dim,int codim,int range,int rank>
SpaceElement<dim,codim,range,rank>::
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



template<int dim,int codim,int range,int rank>
void
SpaceElement<dim,codim,range,rank>::
copy_from(const self_t &elem,
          const CopyPolicy &copy_policy)
{
    if (this != &elem)
    {
        base_t::copy_from(elem,copy_policy);


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



template<int dim,int codim,int range,int rank>
void
SpaceElement<dim,codim,range,rank>::
deep_copy_from(const self_t &elem)
{
    this->copy_from(elem,CopyPolicy::deep);
}



template<int dim,int codim,int range,int rank>
void
SpaceElement<dim,codim,range,rank>::
shallow_copy_from(const self_t &elem)
{
    this->copy_from(elem,CopyPolicy::shallow);
}



template<int dim,int codim,int range,int rank>
auto
SpaceElement<dim,codim,range,rank>::
operator=(const self_t &element) -> self_t &
{
    this->shallow_copy_from(element);
    return (*this);
}



template<int dim,int codim,int range,int rank>
void
SpaceElement<dim,codim,range,rank>::
print_info(LogStream &out) const
{
    base_t::print_info(out);
}

template<int dim,int codim,int range,int rank>
void
SpaceElement<dim,codim,range,rank>::
print_cache_info(LogStream &out) const
{
    base_t::print_cache_info(out);

    Assert(all_sub_elems_cache_ != nullptr,ExcNullPtr());
    all_sub_elems_cache_->print_info(out);
}


#ifdef SERIALIZATION
template<int dim,int codim,int range,int rank>
template<class Archive>
void
SpaceElement<dim,codim,range,rank>::
serialize(Archive &ar, const unsigned int version)
{
    ar &boost::serialization::make_nvp("SpaceElement_base_t",
                                       boost::serialization::base_object<SpaceElementBase<dim>>(*this));

    ar &boost::serialization::make_nvp("all_sub_elems_cache_",all_sub_elems_cache_);
}
#endif // SERIALIZATION



IGA_NAMESPACE_CLOSE


#include <igatools/basis_functions/space_element.inst>


