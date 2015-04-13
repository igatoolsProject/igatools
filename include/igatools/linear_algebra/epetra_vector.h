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

#ifndef __EPETRA_VECTOR_H_
#define __EPETRA_VECTOR_H_

#include <igatools/base/config.h>
#include <igatools/linear_algebra/epetra_map.h>
#include <igatools/linear_algebra/dense_vector.h>
#include <igatools/utils/vector.h>
#include <igatools/base/properties.h>

#include <Epetra_Vector.h>


IGA_NAMESPACE_OPEN

namespace EpetraTools
{

class  Vector : public Epetra_Vector
{
public:
    using Epetra_Vector::Epetra_Vector;

    Size size() const;

    Vector &operator +=(const Vector &vec);

    void add_block(const vector<Index> &vec_id,
                   const DenseVector &local_vector);

    //TODO (pauletti, Apr 3, 2015): both vector<Real> and std::vector<Index>
    // should be replace by a typedef and a proper type for fast comuniction with LA
    vector<Real>
    get_local_coeffs(const std::vector<Index> &global_ids) const;

    void print_info(LogStream &out) const;

};

using VectorPtr = std::shared_ptr<Vector>;

VectorPtr create_vector(MapPtr map);

template <class SpacePtr>
VectorPtr
create_vector(SpacePtr space, const std::string &prop = DofProperties::active)
{
    Epetra_SerialComm comm;
    auto map = EpetraTools::create_map(space, prop, comm);
    return EpetraTools::create_vector(map);
}

};

IGA_NAMESPACE_CLOSE

#endif
