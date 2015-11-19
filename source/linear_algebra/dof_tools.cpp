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


#include <igatools/linear_algebra/dof_tools.h>
#include <igatools/base/exceptions.h>

using std::map;
using std::set;
using std::pair;
using std::shared_ptr;

IGA_NAMESPACE_OPEN

namespace dof_tools
{

void apply_boundary_values(const std::map<Index,Real> &boundary_values,
                           Matrix &matrix,
                           Vector &rhs,
                           Vector &solution)
{
  auto dof = boundary_values.begin();
  const auto dof_end = boundary_values.end();
  const auto &graph = matrix.Graph();
  const auto &map = graph.RowMap();

  Assert(graph.IndicesAreLocal(),ExcMessage("Indices in the graph are not local."));

  SafeSTLVector<int> global_dofs;
  SafeSTLVector<Real> rhs_values;
  SafeSTLVector<Real> solution_values;
  int n_dofs = 0;
  for (; dof != dof_end; ++dof, ++n_dofs)
  {
    const Index row_id = dof->first;
    const Index loc_id = map.LID(dof->first);

    global_dofs.emplace_back(row_id);

    const Real bc_value = dof->second;

    int NumIndices;
    int *Indices;
    graph.ExtractMyRowView(loc_id, NumIndices, Indices);

    int NumEntries;
    double *Values;
    matrix.ExtractGlobalRowView(row_id, NumEntries, Values);


    Real mat_value = NAN;
    for (int i=0; i<NumEntries; ++i)
      if (map.GID(Indices[i]) != row_id)
        Values[i] = 0.;
      else mat_value = Values[i];

    Assert(mat_value != NAN, ExcInternalError());
    rhs_values.emplace_back(bc_value * mat_value);
    solution_values.emplace_back(bc_value);
  }

  rhs.ReplaceGlobalValues(n_dofs,rhs_values.data(),global_dofs.data());
  solution.ReplaceGlobalValues(n_dofs,solution_values.data(),global_dofs.data());
}

}

IGA_NAMESPACE_CLOSE

#include <igatools/linear_algebra/dof_tools.inst>
