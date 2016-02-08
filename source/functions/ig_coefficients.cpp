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

#include <igatools/functions/ig_coefficients.h>


IGA_NAMESPACE_OPEN

IgCoefficients::
IgCoefficients(const std::set<Index> &global_dofs)
{
  for (const auto &dof : global_dofs)
    (*this)[dof] = 0.0;
}

IgCoefficients::
IgCoefficients(const EpetraTools::Vector &epetra_vector)
{
  Assert(false,ExcNotImplemented());
}


const Real &
IgCoefficients::
operator[](const Index global_dof) const
{
  return std::map<Index,Real>::at(global_dof);
}

Real &
IgCoefficients::
operator[](const Index global_dof)
{
  return std::map<Index,Real>::operator[](global_dof);
}

Index
IgCoefficients::
size() const
{
  return std::map<Index,Real>::size();
}

void
IgCoefficients::
print_info(LogStream &out) const
{
  int loc_id = 0;
  for (const auto &dof_value : (*this))
  {
    out << "Coef[loc_id=" << loc_id
        << " , glob_id=" << dof_value.first << "] = "
        << dof_value.second << std::endl;

    ++loc_id;
  }
}


#ifdef SERIALIZATION
template<class Archive>
void
IgCoefficients::
serialize(Archive &ar)
{
  ar &make_nvp("IgCoeff_base_t",
               base_class<std::map<Index,Real>>(this));
}
#endif // SERIALIZATION



IGA_NAMESPACE_CLOSE

#ifdef SERIALIZATION

template void iga::IgCoefficients::serialize(OArchive &);
template void iga::IgCoefficients::serialize(IArchive &);

#endif //SERIALIZATION
