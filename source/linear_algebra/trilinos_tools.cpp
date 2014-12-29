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


#include <igatools/linear_algebra/trilinos_tools.h>

IGA_NAMESPACE_OPEN

#ifdef USE_TRILINOS


namespace trilinos_tools
{
DofsMapPtr build_row_map(const SpaceManager &space_manager,CommPtr comm)
{
    const auto dofs_set = space_manager.get_row_dofs();
    const vector<GO> dofs_vec(dofs_set.begin(),dofs_set.end());

    return Tpetra::createNonContigMap<LO,GO>(dofs_vec,comm);
}

DofsMapPtr build_col_map(const SpaceManager &space_manager,const CommPtr comm)
{
    const auto dofs_set = space_manager.get_col_dofs();
    const vector<GO> dofs_vec(dofs_set.begin(),dofs_set.end());

    return Tpetra::createNonContigMap<LO,GO>(dofs_vec,comm);
}

GraphPtr build_graph(const SpaceManager &space_manager,const DofsMapPtr row_map,const DofsMapPtr col_map)
{
    auto sparsity_pattern_ptr = space_manager.get_sparsity_pattern();
    Assert(sparsity_pattern_ptr!=nullptr,ExcNullPtr());

    const auto &sparsity_pattern = *sparsity_pattern_ptr;
    using LongUInt = long unsigned int;
    Teuchos::ArrayRCP<LongUInt> n_dofs_per_row(sparsity_pattern.get_num_rows());

    Index i = 0;
    for (const auto &map_entry : sparsity_pattern)
        n_dofs_per_row[i++] = map_entry.second.size();



    GraphPtr graph = Teuchos::rcp(new Graph(row_map,col_map,n_dofs_per_row,Tpetra::StaticProfile));
    for (const auto &row : sparsity_pattern)
    {
        const Index row_id = row.first ;
        const auto &cols_id = row.second;

        auto cols_id_vec = vector<Index>(cols_id.begin(),cols_id.end());

        auto cols_id_view = Teuchos::ArrayView<const GO>(std::move(cols_id_vec));

        graph->insertGlobalIndices(row_id,cols_id_view);
    }
    graph->fillComplete(col_map,row_map);

    /*
    Teuchos::RCP<Teuchos::FancyOStream>
          tout = Teuchos::VerboseObjectBase::getDefaultOStream();
    graph->describe(*tout, Teuchos::EVerbosityLevel::VERB_EXTREME);
    graph->print(std::cout);
    //*/
    return graph;
}


MatrixImplPtr build_matrix(GraphPtr graph)
{
    auto matrix = Teuchos::rcp(new MatrixImpl(graph));
    matrix->setAllToScalar(0.0);

    return matrix;
}

};


#endif // #ifdef USE_TRILINOS

IGA_NAMESPACE_CLOSE


