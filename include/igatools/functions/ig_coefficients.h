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

#ifndef __IG_COEFFICIENTS_H
#define __IG_COEFFICIENTS_H


#include <igatools/base/config.h>
#include <igatools/base/logstream.h>

#include <map>
#include <set>


IGA_NAMESPACE_OPEN

/**
 * @brief Coefficients for the IgFunction and IgGridFunction classes.
 *
 * Basically it is a <tt>std::map<Index,Real></tt> in which its <tt>key</tt> is the global dof
 * and the associated <tt>value</tt> is the coefficient associated to the dof.
 *
 * @note We do not use the EpetraTools::Vector because we want the IgCoefficient class to be
 * <em>serializable</em> and to make EpetraTools::Vector serializable is not an easy task
 * (it requires to make serializable all the attributes of EpetraTools::Vector).
 *
 * @ingroup serializable
 */
class IgCoefficients
  : public std::map<Index,Real>
{
public:

  IgCoefficients() = default;

  /**
   * /brief Initialize to zero the coefficients associated with the <tt>global_dofs</tt>
   * used in the input argument.
   */
  IgCoefficients(const std::set<Index> &global_dofs);


  /**
   * Access operator.
   * Returns a reference to the mapped value of the element identified with key @p global_dof.
   * If @p global_dof does not match the key of any element in the container,
   * the function throws an <tt>out_of_range</tt> exception.
   */
  const Real &operator[](const Index global_dof) const;

  /**
   * Access operator.
   * If @p global_dof matches the key of an element in the container, the function returns a reference to its mapped value.
   * If @p global_dof does not match the key of any element in the container,
   * the function inserts a new element with that key and returns a reference to its mapped value.
   */
  Real &operator[](const Index global_dof);

  /** Return the number of coefficients stored in the container. */
  Index size() const;

  void print_info(LogStream &out) const;

private:


#ifdef SERIALIZATION
  /**
   * @name Functions needed for the serialization
   */
  ///@{
  friend class cereal::access;

  template<class Archive>
  void serialize(Archive &ar);
  ///@}
#endif // SERIALIZATION

};


IGA_NAMESPACE_CLOSE


#ifdef SERIALIZATION
CEREAL_SPECIALIZE_FOR_ARCHIVE(IArchive,iga::IgCoefficients,cereal::specialization::member_serialize)
CEREAL_SPECIALIZE_FOR_ARCHIVE(OArchive,iga::IgCoefficients,cereal::specialization::member_serialize)
#endif // SERIALIZATION


#endif // __IG_COEFFICIENTS_H

