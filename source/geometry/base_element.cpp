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

#include <igatools/geometry/base_element.h>
#include <igatools/geometry/unit_element.h>
#include <algorithm>

using std::endl;

IGA_NAMESPACE_OPEN

template <int dim>
BaseElement<dim>::
BaseElement(const Index fi, const TI &ti)
:
flat_index_(fi),
tensor_index_(ti)
{}



template <int dim>
inline
auto
BaseElement<dim>::
get_flat_index() const -> Index
{
    return flat_index_ ;
}



template <int dim>
inline
auto
BaseElement<dim>::
get_tensor_index() const -> TensorIndex<dim>
{
    return tensor_index_ ;
}



template <int dim>
bool
BaseElement<dim>::
operator ==(const BaseElement<dim> &elem) const
{
    return (flat_index_ == elem.flat_index_);
}



template <int dim>
bool
BaseElement<dim>::
operator !=(const BaseElement<dim> &elem) const
{
    return (flat_index_ != elem.flat_index_);
}

template <int dim>
bool
BaseElement<dim>::
operator <(const BaseElement<dim> &elem) const
{
    return (flat_index_ < elem.flat_index_);
}

template <int dim>
bool
BaseElement<dim>::
operator >(const BaseElement<dim> &elem) const
{
    return (flat_index_ > elem.flat_index_);
}


template <int dim>
void
BaseElement<dim>::
print_info(LogStream &out) const
{
    out << "Flat id = "   << flat_index_ << "    ";
    out << "Tensor id = " << tensor_index_ << endl;
}



#ifdef SERIALIZATION
template <int dim>
template<class Archive>
void
BaseElement<dim>::
serialize(Archive &ar, const unsigned int version)
{
    auto non_const_grid = std::const_pointer_cast<CartesianGrid<dim>>(grid_);
    ar &boost::serialization::make_nvp("grid_",non_const_grid);
    grid_ = non_const_grid;
    Assert(grid_ != nullptr,ExcNullPtr());

    ar &boost::serialization::make_nvp("flat_index_",flat_index_);

    ar &boost::serialization::make_nvp("tensor_index_",tensor_index_);
}
#endif // SERIALIZATION


IGA_NAMESPACE_CLOSE

#include <igatools/geometry/base_element.inst>
