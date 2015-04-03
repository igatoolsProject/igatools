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

#ifndef __EPETRA_MAP_H_
#define __EPETRA_MAP_H_

#include <igatools/base/config.h>

#include <igatools/linear_algebra/dense_vector.h>
#include <igatools/linear_algebra/dense_matrix.h>
#include <igatools/utils/vector.h>

#include <Epetra_SerialComm.h>
#include <Epetra_Map.h>


IGA_NAMESPACE_OPEN

namespace EpetraTools
{
    using Comm = Epetra_Comm;
    using CommPtr = std::shared_ptr<const Comm>;

    using Map = Epetra_Map;
    using MapPtr = std::shared_ptr<Map>;

    template<class SpacePtr>
    MapPtr create_map(const SpacePtr space,
                      const std::string &property,
                      Comm &comm)
    {
        const auto dof_dist = space->get_dof_distribution();
        const auto dofs = dof_dist->get_dofs_id_same_property(property);
        //TODO (pauletti, Mar 28, 2015): this is double copy of data
        const vector<Index> dofs_vec(dofs.begin(), dofs.end());
        auto map = std::make_shared<Map>(-1, dofs_vec.size(), dofs_vec.data(), 0, comm);
        return map;
    }

};

IGA_NAMESPACE_CLOSE

#endif
