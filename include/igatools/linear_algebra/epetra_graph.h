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

#ifndef __EPETRA_GRAPH_H_
#define __EPETRA_GRAPH_H_

#include <igatools/base/config.h>
#include <igatools/linear_algebra/epetra_map.h>

#include <Epetra_CrsGraph.h>

IGA_NAMESPACE_OPEN

namespace EpetraTools
{
using Graph = Epetra_CrsGraph;
using GraphPtr = std::shared_ptr<Graph>;


template<class RowSpacePtr, class ColSpacePtr>
GraphPtr
create_graph(const RowSpacePtr row_space, const std::string &row_property,
             const ColSpacePtr col_space, const std::string &col_property,
             MapPtr row_map_, MapPtr col_map_)
{
    /*
    LogStream out;
    out.begin_item("row space");
    row_space->print_info(out);
    out.end_item();
    //*/
    const auto n_rows = row_map_->NumMyElements();
    const auto dof_distribution_row_space = row_space->get_dof_distribution();
    /*
    out.begin_item("Dof dof_distribution_row_space");
    dof_distribution_row_space->print_info(out);
    out.end_item();
    out << "Dof distribution n.dofs = " << dof_distribution_row_space->get_num_dofs(row_property) << std::endl;
    out << "n_rows = " << n_rows << std::endl;
    //*/
    Assert(dof_distribution_row_space->get_num_dofs(row_property) == n_rows,
           ExcDimensionMismatch(dof_distribution_row_space->get_num_dofs(row_property),n_rows));

    SafeSTLVector<SafeSTLVector<Index>> loc_dofs(n_rows);
    auto r_elem = row_space->begin();
    auto c_elem = col_space->begin();
    const auto end = row_space->end();
    for (; r_elem != end; ++r_elem, ++c_elem)
    {
        auto r_dofs = r_elem->get_local_to_global(row_property);
        auto c_dofs = c_elem->get_local_to_global(col_property);
        for (auto &r_dof : r_dofs)
        {
//          const int loc_r_dof = row_map_->LID(r_dof);
//            auto &dof_vec = loc_dofs[loc_r_dof];
            auto &dof_vec = loc_dofs[r_dof];
            dof_vec.insert(dof_vec.begin(), c_dofs.begin(), c_dofs.end());
        }
    }

    SafeSTLVector<Size> n_dofs_per_row(n_rows);
    {
        Index j=0;
        for (auto &dofs : loc_dofs)
        {
            std::sort(dofs.begin(), dofs.end());
            auto it = std::unique(dofs.begin(), dofs.end());
            dofs.resize(std::distance(dofs.begin(),it));
            n_dofs_per_row[j] = dofs.size();
            ++j;
        }
    }

    const bool is_static_profile = true;
    auto graph_ = std::make_shared<Graph>(Epetra_DataAccess::Copy,
                                          *row_map_, *col_map_,
                                          n_dofs_per_row.data(),
                                          is_static_profile);

    Index j=0;
    for (auto &dofs : loc_dofs)
    {
        const Index row_id = row_map_->GID(j);
        graph_->InsertGlobalIndices(row_id, dofs.size(), dofs.data());
        ++j;
    }

    int res = graph_->FillComplete(*col_map_,*row_map_);
    AssertThrow(res==0, ExcMessage(" "));

    return graph_;
}

};

IGA_NAMESPACE_CLOSE

#endif
