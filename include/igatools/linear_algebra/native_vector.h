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

#ifndef __NATIVE_VECTOR_H_
#define __NATIVE_VECTOR_H_

#include <igatools/base/config.h>

IGA_NAMESPACE_OPEN

#if 0
namespace NativeTools
{
class Vector : public std::map<Index,Real>
{
public:
  using std::map<Index,Real>::map;

  IgCoefficients(const std::set<Index> &global_dofs, const SafeSTLVector<Real> &coeffs)
  {
    Assert(Index(global_dofs.size()) == coeffs.size(),
           ExcDimensionMismatch(global_dofs.size(),coeffs.size()));

    int i = 0;
    for (const auto dof : global_dofs)
      (*this)[dof] = coeffs[i++];
  }

  template <class Space>
  IgCoefficients(
    const Space &space,
    const std::string &dofs_property,
    const SafeSTLVector<Real> &coeffs)
    :
    IgCoefficients(space.get_ptr_const_dof_distribution()->get_dofs_id_same_property(dofs_property),coeffs)
  {}

  template <class Space>
  IgCoefficients(
    const Space &space,
    const std::string &dofs_property)
    :
    IgCoefficients(
     space.get_ptr_const_dof_distribution()->get_dofs_id_same_property(dofs_property),
     SafeSTLVector<Real>(space.get_ptr_const_dof_distribution()->
                         get_dofs_id_same_property(dofs_property).size(),0.0))
  {}


  Real &operator()(const Index &global_dof)
  {
#ifdef NDEBUG
    return (*this)[global_dof];
#else
    return (*this).at(global_dof);
#endif
  }

  const Real &operator()(const Index &global_dof) const
  {
#ifdef NDEBUG
    return (*this)[global_dof];
#else
    return (*this).at(global_dof);
#endif
  }


  Size size() const
  {
    return std::map<Index,Real>::size();
  }

  IgCoefficients &operator+=(const IgCoefficients &coeffs)
  {
    Assert(this->size() == coeffs.size(),ExcDimensionMismatch(this->size(),coeffs.size()));
#ifdef NDEBUG
    for (const auto &c : coeffs)
      (*this)[c.first] += c.second;
#else
    for (const auto &c : coeffs)
      (*this).at(c.first) += c.second;
#endif

    return *this;
  }

  SafeSTLVector<Real> get_local_coeffs(const SafeSTLVector<Index> &elem_dofs) const
  {

    SafeSTLVector<Real> loc_coeff;
    for (const auto &dof : elem_dofs)
      loc_coeff.emplace_back((*this)(dof));
    return  loc_coeff;
  }

  void print_info(LogStream &out) const
  {
    int i = 0;
    for (const auto &c : (*this))
    {
      out << "coeff[local_id=" << i << ", global_id=" << c.first << "] = " << c.second << std::endl;
      ++i;
    }
  }
};
};
#endif

IGA_NAMESPACE_CLOSE

#endif
