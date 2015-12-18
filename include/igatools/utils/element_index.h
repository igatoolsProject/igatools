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
  /**
   * Default constructor. Sets the element index to 0.
   */
  ElementIndex();

  ElementIndex(const int flat_id, const TensorIndex<dim> &tensor_id);

  /**
   * Copy constructor.
   */
  ElementIndex(const ElementIndex<dim> &elem_id) = default;

  /**
   * Move constructor.
   */
  ElementIndex(ElementIndex<dim> &&elem_id) = default;

  /**
   * Destructor.
   */
  ~ElementIndex() = default;

  /**
   * Copy assignment operator.
   */
  ElementIndex<dim> &operator=(const ElementIndex<dim> &elem_id) = delete;

  /**
   * Move assignment operator.
   */
  ElementIndex<dim> &operator=(ElementIndex<dim> &&elem_id) = default;


  const TensorIndex<dim> &
  get_tensor_index() const;

  int get_flat_index() const;


  bool operator==(const ElementIndex<dim> &elem) const;

  bool operator!=(const ElementIndex<dim> &elem) const;

  bool operator<(const ElementIndex<dim> &elem) const;

  bool operator>(const ElementIndex<dim> &elem) const;

  void print_info(LogStream &out) const;

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
#include <igatools/utils/element_index.serial>
#endif // SERIALIZATION


#endif // __ELEMENT_INDEX_H_



