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

#include <igatools/linear_algebra/dense_vector.h>
#include <igatools/base/exceptions.h>

IGA_NAMESPACE_OPEN


DenseVector &
DenseVector::
operator=(const Real value)
{
  using zero_vector = boost::numeric::ublas::zero_vector<Real>;
  Assert(value==0, ExcNonZero());
  *this = zero_vector(this->size());
  return *this;
}

int
DenseVector::
size() const
{
  return int(boost::numeric::ublas::vector<Real>::size());
}

IGA_NAMESPACE_CLOSE

