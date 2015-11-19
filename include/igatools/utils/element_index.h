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

#ifndef __ELEMENT_INDEX_H_
#define __ELEMENT_INDEX_H_

#include <igatools/utils/tensor_index.h>

IGA_NAMESPACE_OPEN

template <int dim>
class ElementIndex
{
public:
  ElementIndex() = delete;

  ElementIndex(const int flat_id, const TensorIndex<dim> &tensor_id)
    :
    flat_id_(flat_id),
    tensor_id_(tensor_id)
  {}


  const TensorIndex<dim> &
  get_tensor_index() const
  {
    return tensor_id_;
  }

  const int
  get_flat_index() const
  {
    return flat_id_;
  }


  bool operator==(const ElementIndex<dim> &elem) const
  {
    const bool same_tid = (this->tensor_id_ == elem.tensor_id_);
    const bool same_fid = (this->flat_id_ == elem.flat_id_);
    return (same_tid && same_fid);
  }

  bool operator!=(const ElementIndex<dim> &elem) const
  {
    const bool different_tid = (this->tensor_id_ != elem.tensor_id_);
    const bool different_fid = (this->flat_id_ != elem.flat_id_);
    return (different_tid || different_fid);
  }

  bool operator<(const ElementIndex<dim> &elem) const
  {
    return (this->tensor_id_ < elem.tensor_id_);
  }

  bool operator>(const ElementIndex<dim> &elem) const
  {
    return (this->tensor_id_ > elem.tensor_id_);
  }

  void print_info(LogStream &out) const
  {
    out << "Flat ID: " << flat_id_ << "    Tensor ID: " << tensor_id_ << std::endl;
//    tensor_id_.print_info(out);
  }
private:
  int flat_id_ = 0;
  TensorIndex<dim> tensor_id_;


private:
#ifdef SERIALIZATION
  /**
   * @name Functions needed for serialization
   * @see <a href="http://uscilab.github.io/cereal/serialization_functions.html">Cereal serialization</a>
   */
  ///@{
  friend class cereal::access;

  template<class Archive>
  void serialize(Archive &ar)
  {
    ar &make_nvp("flat_id_",flat_id_);
    ar &make_nvp("tensor_id_",tensor_id_);
  }
  ///@}
#endif // SERIALIZATION

};

template <int dim>
LogStream &
operator<<(LogStream &out, const ElementIndex<dim> &elem_id)
{
  out << "Flat ID: " << elem_id.get_flat_index()
      << "    Tensor ID: " << elem_id.get_tensor_index();
//  out << elem_id.get_tensor_index();
  return out;
}


IGA_NAMESPACE_CLOSE





#ifdef SERIALIZATION
using ElementIndexAlias0 = iga::ElementIndex<0>;
CEREAL_SPECIALIZE_FOR_ALL_ARCHIVES(ElementIndexAlias0,cereal::specialization::member_serialize);
using ElementIndexAlias1 = iga::ElementIndex<1>;
CEREAL_SPECIALIZE_FOR_ALL_ARCHIVES(ElementIndexAlias1,cereal::specialization::member_serialize);
using ElementIndexAlias2 = iga::ElementIndex<2>;
CEREAL_SPECIALIZE_FOR_ALL_ARCHIVES(ElementIndexAlias2,cereal::specialization::member_serialize);
using ElementIndexAlias3 = iga::ElementIndex<3>;
CEREAL_SPECIALIZE_FOR_ALL_ARCHIVES(ElementIndexAlias3,cereal::specialization::member_serialize);

//#include <igatools/utils/element_index.serialization>
#endif // SERIALIZATION



#endif // __ELEMENT_INDEX_H_



