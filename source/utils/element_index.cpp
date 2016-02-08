//-+--------------------------------------------------------------------
// Igatools a general purpose Isogeometric analysis library.
// Copyright (C) 2012-2016  by the igatools authors (see authors.txt).
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


#include <igatools/utils/element_index.h>


IGA_NAMESPACE_OPEN

template <int dim>
ElementIndex<dim>::
ElementIndex()
  :
  flat_id_(0),
  tensor_id_(TensorIndex<dim>())
{}

template <int dim>
ElementIndex<dim>::
ElementIndex(const int flat_id, const TensorIndex<dim> &tensor_id)
  :
  flat_id_(flat_id),
  tensor_id_(tensor_id)
{}


template <int dim>
const TensorIndex<dim> &
ElementIndex<dim>::
get_tensor_index() const
{
  return tensor_id_;
}

template <int dim>
int
ElementIndex<dim>::
get_flat_index() const
{
  return flat_id_;
}


template <int dim>
bool
ElementIndex<dim>::
operator==(const ElementIndex<dim> &elem) const
{
  const bool same_tid = (this->tensor_id_ == elem.tensor_id_);
  const bool same_fid = (this->flat_id_ == elem.flat_id_);
  return (same_tid && same_fid);
}

template <>
bool
ElementIndex<0>::
operator==(const ElementIndex<0> &elem) const
{
  // If dim==0 the element index is always the same
  return true;
}

template <int dim>
bool
ElementIndex<dim>::
operator!=(const ElementIndex<dim> &elem) const
{
  const bool different_tid = (this->tensor_id_ != elem.tensor_id_);
  const bool different_fid = (this->flat_id_ != elem.flat_id_);
  return (different_tid || different_fid);
}

template <>
bool
ElementIndex<0>::
operator!=(const ElementIndex<0> &elem) const
{
  // If dim==0 the element index is always the same
  return false;
}

template <int dim>
bool
ElementIndex<dim>::
operator<(const ElementIndex<dim> &elem) const
{
  return (this->tensor_id_ < elem.tensor_id_);
}

template <int dim>
bool
ElementIndex<dim>::
operator>(const ElementIndex<dim> &elem) const
{
  return (this->tensor_id_ > elem.tensor_id_);
}

template <int dim>
void
ElementIndex<dim>::
print_info(LogStream &out) const
{
  out << "Flat ID: " << flat_id_ << "    Tensor ID: " << tensor_id_ << std::endl;
//    tensor_id_.print_info(out);
}


#ifdef SERIALIZATION
template <int dim>
template<class Archive>
void
ElementIndex<dim>::
serialize(Archive &ar)
{
  ar &make_nvp("flat_id_",flat_id_);
  ar &make_nvp("tensor_id_",tensor_id_);
}
#endif

IGA_NAMESPACE_CLOSE

#include <igatools/utils/element_index.inst>
