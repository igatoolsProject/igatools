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

#include <igatools/geometry/domain_element.h>
#include <igatools/functions/function.h>
#include <igatools/functions/function_element.h>



IGA_NAMESPACE_OPEN


template<int dim, int codim, int range, int rank,  class ContainerType_>
FunctionElementBase<dim, codim, range, rank, ContainerType_>::
FunctionElementBase(const std::shared_ptr<ContainerType_> func,
                    const ListIt &index,
                    const PropId &prop)
  :
  func_(func),
  domain_elem_(func->get_domain()->create_element(index,prop))
{}


#if 0
template<int dim, int codim, int range, int rank,  class ContainerType_>
FunctionElementBase<dim, codim, range, rank, ContainerType_>::
FunctionElementBase(const self_t &elem,
                    const CopyPolicy &copy_policy)
  :
  func_(elem.func_)
{
  if (copy_policy == CopyPolicy::shallow)
  {
    local_cache_ = elem.local_cache_;
    domain_elem_ = elem.domain_elem_;
  }
  else
  {
    local_cache_ = std::make_shared<AllSubElementsCache<Cache>>(*elem.local_cache_);
    domain_elem_ = std::make_shared<DomainElem>(*elem.domain_elem_,CopyPolicy::deep);
  }
}
#endif


template<int dim, int codim, int range, int rank,  class ContainerType_>
auto
FunctionElementBase<dim, codim, range, rank, ContainerType_>::
get_domain_element() const -> const DomainElem &
{
  return *domain_elem_;
}

template<int dim, int codim, int range, int rank,  class ContainerType_>
auto
FunctionElementBase<dim, codim, range, rank, ContainerType_>::
get_domain_element() -> DomainElem &
{
  return *domain_elem_;
}



template<int dim, int codim, int range, int rank,  class ContainerType_>
bool
FunctionElementBase<dim, codim, range, rank, ContainerType_>::
operator==(const self_t &a) const
{
  Assert(func_ == a.func_,
         ExcMessage("The elements cannot be compared because defined with different functions."));
  return (*domain_elem_ == *(a.domain_elem_));
}


template<int dim, int codim, int range, int rank,  class ContainerType_>
bool
FunctionElementBase<dim, codim, range, rank, ContainerType_>::
operator!=(const self_t &a) const
{
  Assert(func_ == a.func_,
         ExcMessage("The elements cannot be compared because defined with different functions."));
  return (*domain_elem_ != *(a.domain_elem_));
}

template<int dim, int codim, int range, int rank,  class ContainerType_>
bool
FunctionElementBase<dim, codim, range, rank, ContainerType_>::
operator<(const self_t &a) const
{
  Assert(func_ == a.func_,
         ExcMessage("The elements cannot be compared because defined with different functions."));
  return (*domain_elem_ < *(a.domain_elem_));
}


template<int dim, int codim, int range, int rank,  class ContainerType_>
bool
FunctionElementBase<dim, codim, range, rank, ContainerType_>::
operator>(const self_t &a) const
{
  Assert(func_ == a.func_,
         ExcMessage("The elements cannot be compared because defined with different functions."));
  return (*domain_elem_ > *(a.domain_elem_));
}



template<int dim, int codim, int range, int rank,  class ContainerType_>
void
FunctionElementBase<dim, codim, range, rank, ContainerType_>::
print_info(LogStream &out) const
{
  Assert(false, ExcNotImplemented());
}



template<int dim, int codim, int range, int rank,  class ContainerType_>
void
FunctionElementBase<dim, codim, range, rank, ContainerType_>::
print_cache_info(LogStream &out) const
{
  local_cache_->print_info(out);
}

#if 0
#ifdef SERIALIZATION
template<int dim, int codim, int range, int rank,  class ContainerType_>
template<class Archive>
void
FunctionElementBase<dim, codim, range, rank, ContainerType_>::
serialize(Archive &ar, const unsigned int version)
{
  AssertThrow(false,ExcNotImplemented());

  ar &boost::serialization::make_nvp("FunctionElement_base_t",
                                     boost::serialization::base_object<GridElement<dim>>(*this));

  ar &boost::serialization::make_nvp("all_sub_elems_cache_",local_cache_);

  ar &boost::serialization::make_nvp("func_",func_);
  ar &boost::serialization::make_nvp("grid_elem_",grid_elem_);

  ar &boost::serialization::make_nvp("phys_domain_elem_",domain_elem_);
}
#endif // SERIALIZATION
#endif

IGA_NAMESPACE_CLOSE

#include <igatools/functions/function_element.inst>

