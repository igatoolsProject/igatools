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

#include <igatools/linear_algebra/epetra_vector.h>

IGA_NAMESPACE_OPEN

namespace EpetraTools
{

Size Vector::size() const
{
  return GlobalLength();
}



Vector &Vector::operator +=(const Vector &vec)
{
  Update(1., vec, 1.);
  return *this;
}



void Vector::add_block(const SafeSTLVector<Index> &vec_id,
                       const DenseVector &local_vector)
{
  const auto   NumEntries = vec_id.size();
  const double *Values    = &(local_vector.data()[0]);
  const int    *Indices   = vec_id.data();

  Epetra_Vector::SumIntoGlobalValues(NumEntries, Values, Indices);
}



//TODO (pauletti, Apr 3, 2015): both SafeSTLVector<Real> and std::vector<Index>
// should be replace by a typedef and a proper type for fast comuniction with LA
SafeSTLVector<Real>
Vector::get_local_coeffs(const std::vector<Index> &global_ids) const
{
  SafeSTLVector<Real> local_coefs;
  for (const auto &global_id : global_ids)
    local_coefs.emplace_back((*this)[global_id]);

  return local_coefs;
}



void Vector::print_info(LogStream &out) const
{
  using std::endl;
  out << "-----------------------------" << endl;

  const Index n_entries = GlobalLength();
  const auto &map = Map();
  out << "Global_ID        Value" << endl;

  for (Index i = 0 ; i < n_entries ; ++i)
    out << map.GID(i) << "        " << (*this)[i] << std::endl ;

  out << "-----------------------------" << endl;
}



VectorPtr
create_vector(const Map &map)
{
  return std::make_shared<Vector>(map);
}

}

IGA_NAMESPACE_CLOSE
