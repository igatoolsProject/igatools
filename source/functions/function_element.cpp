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

#include <igatools/geometry/mapping_element.h>

IGA_NAMESPACE_OPEN


template<int dim, int codim, int range, int rank>
FunctionElement<dim, codim, range, rank>::
FunctionElement(const std::shared_ptr<const Func> func,
                const Index elem_index)
    :
    func_(std::const_pointer_cast<Func>(func)),
	grid_elem_(func->get_grid()->create_element(elem_index))
{
    Assert(func_ != nullptr, ExcNullPtr());
    Assert(grid_elem_ != nullptr, ExcNullPtr());

	auto phys_domain = PhysDomain::create(func);
	phys_domain_elem_ = phys_domain->create_element(elem_index);
    Assert(phys_domain_elem_ != nullptr, ExcNullPtr());
}


template<int dim, int codim, int range, int rank>
FunctionElement<dim, codim, range, rank>::
FunctionElement(const FunctionElement<dim,codim,range,rank> &elem,
                const CopyPolicy &copy_policy)
    :
    func_(elem.func_)
{
    if (copy_policy == CopyPolicy::shallow)
    {
        all_sub_elems_cache_ = elem.all_sub_elems_cache_;
        grid_elem_ = elem.grid_elem_;
        phys_domain_elem_ = elem.phys_domain_elem_;
    }
    else
    {
        all_sub_elems_cache_ = std::make_shared<AllSubElementsCache<Cache>>(*elem.all_sub_elems_cache_);
        grid_elem_ = std::make_shared<GridElem>(*elem.grid_elem_,CopyPolicy::deep);
        phys_domain_elem_ = std::make_shared<PhysDomainElem>(*elem.phys_domain_elem_,CopyPolicy::deep);
    }
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


template<int dim, int codim, int range, int rank>
auto
FunctionElement<dim, codim, range, rank>::
get_grid_element() const -> const GridElem &
{
	return *grid_elem_;
}




template<int dim, int codim, int range, int rank>
bool
FunctionElement<dim, codim, range, rank>::
operator==(const self_t &a) const
{
	Assert (this->get_grid() == a.get_grid(),
			ExcMessage("The elements cannot be compared because defined on different grids."));
	Assert (func_ == a.func_,
			ExcMessage("The elements cannot be compared because defined with different functions."));
	return (this->get_flat_index() == a.get_flat_index());
}


template<int dim, int codim, int range, int rank>
bool
FunctionElement<dim, codim, range, rank>::
operator!=(const self_t &a) const
{
	Assert (this->get_grid() == a.get_grid(),
			ExcMessage("The elements cannot be compared because defined on different grids."));
	Assert (func_ == a.func_,
			ExcMessage("The elements cannot be compared because defined with different functions."));
	return (this->get_flat_index() != a.get_flat_index());
}

template<int dim, int codim, int range, int rank>
bool
FunctionElement<dim, codim, range, rank>::
operator<(const self_t &a) const
{
	Assert (this->get_grid() == a.get_grid(),
			ExcMessage("The elements cannot be compared because defined on different grids."));
	Assert (func_ == a.func_,
			ExcMessage("The elements cannot be compared because defined with different functions."));
	return (this->get_flat_index() < a.get_flat_index());
}


template<int dim, int codim, int range, int rank>
bool
FunctionElement<dim, codim, range, rank>::
operator>(const self_t &a) const
{
	Assert (this->get_grid() == a.get_grid(),
			ExcMessage("The elements cannot be compared because defined on different grids."));
	Assert (func_ == a.func_,
			ExcMessage("The elements cannot be compared because defined with different functions."));
	return (this->get_flat_index() > a.get_flat_index());
}


template<int dim, int codim, int range, int rank>
void
FunctionElement<dim, codim, range, rank>::
move_to(const Index flat_index)
{
	grid_elem_->move_to(flat_index);
	phys_domain_elem_->move_to(flat_index);
}


template<int dim, int codim, int range, int rank>
Index
FunctionElement<dim, codim, range, rank>::
get_flat_index() const
{
	return grid_elem_->get_flat_index();
}

template<int dim, int codim, int range, int rank>
TensorIndex<dim>
FunctionElement<dim, codim, range, rank>::
get_tensor_index() const
{
	return grid_elem_->get_tensor_index();
}

template<int dim, int codim, int range, int rank>
void
FunctionElement<dim, codim, range, rank>::
print_info(LogStream &out) const
{
	grid_elem_->print_info(out);
}

template<int dim, int codim, int range, int rank>
void
FunctionElement<dim, codim, range, rank>::
print_cache_info(LogStream &out) const
{
	grid_elem_->print_cache_info(out);
}


template<int dim, int codim, int range, int rank>
std::shared_ptr<const CartesianGrid<dim>>
FunctionElement<dim, codim, range, rank>::
get_grid() const
{
	return grid_elem_->get_grid();
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
    ar &boost::serialization::make_nvp("grid_elem_",grid_elem_);
    ar &boost::serialization::make_nvp("phys_domain_elem_",phys_domain_elem_);
}
#endif // SERIALIZATION


IGA_NAMESPACE_CLOSE

#include <igatools/functions/function_element.inst>

