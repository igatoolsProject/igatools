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


#include <igatools/functions/function_element.h>
#include <igatools/functions/function.h>


IGA_NAMESPACE_OPEN


template<int dim, int codim, int range, int rank>
FunctionElement<dim, codim, range, rank>::
FunctionElement(const std::shared_ptr<const Func> func,
                const Index elem_index)
    :
    CartesianGridElement<dim>(func->get_grid(),elem_index),
    func_(std::const_pointer_cast<Func>(func))
{
    Assert(func_ != nullptr ,ExcNullPtr());
}


template<int dim, int codim, int range, int rank>
FunctionElement<dim, codim, range, rank>::
FunctionElement(const FunctionElement<dim,codim,range,rank> &elem,
                const CopyPolicy &copy_policy)
    :
    CartesianGridElement<dim>(elem,copy_policy),
    func_(elem.func_)
{
    if (copy_policy == CopyPolicy::shallow)
        all_sub_elems_cache_ = elem.all_sub_elems_cache_;
    else
        all_sub_elems_cache_ = std::make_shared<AllSubElementsCache<Cache>>(*elem.all_sub_elems_cache_);
}


template<int dim, int codim, int range, int rank>
FunctionElement<dim,codim,range,rank> &
FunctionElement<dim, codim, range, rank>::
operator=(const FunctionElement<dim,codim,range,rank> &element)
{
    shallow_copy_from(element);
    return *this;
}


template<int dim, int codim, int range, int rank>
std::shared_ptr<FunctionElement<dim,codim,range,rank> >
FunctionElement<dim, codim, range, rank>::
clone() const
{
    auto elem = std::make_shared<FunctionElement<dim,codim,range,rank> >(*this,CopyPolicy::deep);
    Assert(elem != nullptr, ExcNullPtr());
    return elem;
}


template<int dim, int codim, int range, int rank>
ValueFlags
FunctionElement<dim, codim, range, rank>::
get_valid_flags()
{
    return cacheutils::get_valid_flags_from_cache_type(CType());
}


#ifdef SERIALIZATION
template<int dim, int codim, int range, int rank>
template<class Archive>
void
FunctionElement<dim, codim, range, rank>::
serialize(Archive &ar, const unsigned int version)
{
    ar &boost::serialization::make_nvp("FunctionElement_base_t",
                                       boost::serialization::base_object<CartesianGridElement<dim>>(*this));

    ar &boost::serialization::make_nvp("all_sub_elems_cache_",all_sub_elems_cache_);

    ar &boost::serialization::make_nvp("func_",func_);
}
#endif // SERIALIZATION


IGA_NAMESPACE_CLOSE

#include <igatools/functions/function_element.inst>

