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

#include <igatools/linear_algebra/epetra_graph.h>

IGA_NAMESPACE_OPEN

#ifdef IGATOOLS_USES_TRILINOS

namespace EpetraTools
{

GraphPtr
create_graph(const std::map<Index,std::set<Index>> &dofs_connectivity,
             const Comm &comm)
{
  const int n_rows = dofs_connectivity.size();
  SafeSTLVector<Size> n_dofs_per_row(n_rows);

  std::set<Index> row_all_dofs;
  std::set<Index> col_all_dofs;

  Index j = 0;
  for (const auto &row_id_and_dofs : dofs_connectivity)
  {
    const auto &col_dofs = row_id_and_dofs.second;
    n_dofs_per_row[j++] = col_dofs.size();

    row_all_dofs.insert(row_id_and_dofs.first);
    col_all_dofs.insert(col_dofs.begin(),col_dofs.end());
  }
  const auto row_map = create_map(row_all_dofs,comm);
  const auto col_map = create_map(col_all_dofs,comm);


  const bool is_static_profile = true;
  auto graph = std::make_shared<Graph>(Epetra_DataAccess::Copy,
                                       *row_map, *col_map,
                                       n_dofs_per_row.data(),
                                       is_static_profile);

  for (const auto &row_id_and_dofs : dofs_connectivity)
  {
    const Index row_id = row_id_and_dofs.first;
    SafeSTLVector<Index> cols_id(row_id_and_dofs.second.begin(),row_id_and_dofs.second.end());
    graph->InsertGlobalIndices(row_id, cols_id.size(), cols_id.data());
  }

  int res = graph->FillComplete(*col_map,*row_map);
  AssertThrow(res == 0, ExcMessage("Error raised by Epetra_CrsGraph::FillComplete()"));

  return graph;
}
}

#endif //IGATOOLS_USES_TRILINOS

IGA_NAMESPACE_CLOSE
