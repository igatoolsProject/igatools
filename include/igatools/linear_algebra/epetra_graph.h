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


/**
 * Create an Epetra_CrsGraph object (wrapped by a shared pointer) from the global @p dofs_connectivity.
 *
 * The @p dofs_connectivity is a std:map in which the key is the global row id and the associated value
 * is a std::set containing the global columns id associated to the row.
 *
 */
GraphPtr
create_graph(const std::map<Index,std::set<Index>> &dofs_connectivity,
             const Comm &comm);




template<class RowSpace, class ColSpace>
GraphPtr
create_graph(const RowSpace &row_space, const std::string &row_property,
             const ColSpace &col_space, const std::string &col_property,
             const Comm &comm)
{
  std::map<Index,std::set<Index>> dofs_connectivity;

  Assert(row_space.get_grid() == col_space.get_grid(),
         ExcMessage("Row and column basis built on different grids."));

  auto r_elem = row_space.begin();
  auto c_elem = col_space.begin();
  const auto r_end = row_space.end();

  LogStream myout;
  for (; r_elem != r_end ;)
  {
#if 0
    myout.begin_item("Row elem");
    r_elem->get_index().print_info(myout);
    myout.end_item();

    myout.begin_item("Col elem");
    c_elem->get_index().print_info(myout);
    myout.end_item();
#endif
    const auto r_dofs = r_elem->get_local_to_global(row_property);
    const auto c_dofs = c_elem->get_local_to_global(col_property);
    for (auto &r_dof : r_dofs)
      dofs_connectivity[r_dof].insert(c_dofs.begin(),c_dofs.end());

    ++r_elem;

    ++c_elem;
  }

  return create_graph(dofs_connectivity,comm);
}

}

IGA_NAMESPACE_CLOSE

#endif
