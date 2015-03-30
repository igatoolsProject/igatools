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

#ifndef __MATRIX_MAP_H_
#define __MATRIX_MAP_H_

#include <igatools/base/config.h>
#include <igatools/linear_algebra/dense_vector.h>

#include <Epetra_SerialComm.h>
#include <Epetra_Map.h>
#include <Epetra_CrsMatrix.h>

#include <Epetra_Vector.h>



IGA_NAMESPACE_OPEN

namespace EpetraTools
{

class  Vector : public Epetra_Vector
{
public:
    using Epetra_Vector::Epetra_Vector;

    void add_block(const vector<Index> &vec_id,
                   const DenseVector &local_vector)
    {
        const auto   NumEntries = vec_id.size();
        const double *Values    = &(local_vector.data()[0]);
        const int    *Indices   = vec_id.data();

        Epetra_Vector::SumIntoGlobalValues(NumEntries, Values, Indices);

    }
};


class  Matrix : public Epetra_CrsMatrix
{
public:
    using Epetra_CrsMatrix::Epetra_CrsMatrix;

    void add_block(const vector<Index> &rows_id,
                   const vector<Index> &cols_id,
                   const DenseMatrix &loc_matrix)
    {
        const auto n_rows = rows_id.size();
        const auto n_cols = cols_id.size();

        for (int i = 0 ; i < n_rows ; ++i)
        {
            const double *i_row_data =  &(loc_matrix.data()[i*n_cols]);
            SumIntoGlobalValues(rows_id[i], n_cols, i_row_data, cols_id.data());
        }
    }
};

    using Comm = Epetra_Comm;
    using CommPtr = std::shared_ptr<const Comm>;

    using Map = Epetra_Map;
    using MapPtr = std::shared_ptr<Map>;

    using Graph = Epetra_CrsGraph;
    using GraphPtr = std::shared_ptr<Graph>;

    //using Matrix = Epetra_CrsMatrix;
    using MatrixPtr = std::shared_ptr<Matrix>;

    //using Vector = Epetra_MultiVector;
    using VectorPtr = std::shared_ptr<Vector>;


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

    template<class RowSpacePtr, class ColSpacePtr>
    GraphPtr
    create_graph(const RowSpacePtr row_space, const std::string &row_property,
                 const ColSpacePtr col_space, const std::string &col_property,
                 MapPtr row_map_, MapPtr col_map_)
    {
        const auto n_rows = row_map_->NumMyElements();
        vector<vector<Index>> loc_dofs(n_rows);
        auto r_elem = row_space->begin();
        auto c_elem = col_space->begin();
        const auto end = row_space->end();
        for (;r_elem != end; ++r_elem, ++c_elem)
        {
            auto r_dofs = r_elem->get_local_to_global(row_property);
            auto c_dofs = c_elem->get_local_to_global(col_property);
            for(auto &r_dof : r_dofs)
            {
                auto &dof_vec = loc_dofs[row_map_->LID(r_dof)];
                dof_vec.insert(dof_vec.begin(), c_dofs.begin(), c_dofs.end());
            }
        }

        vector<Size> n_dofs_per_row(n_rows);
        {
            Index j=0;
            for (auto &dofs : loc_dofs)
            {
                std::sort(dofs.begin(), dofs.end());
                auto it = std::unique (dofs.begin(), dofs.end());
                dofs.resize(std::distance(dofs.begin(),it) );
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
        Assert(res==0, ExcMessage(" "));

        return graph_;
    }


    MatrixPtr
    create_matrix(GraphPtr graph)
    {
        return std::make_shared<Matrix>(Epetra_DataAccess::Copy, *graph);
    }


    VectorPtr
    create_vector(MapPtr map)
    {
        return std::make_shared<Vector>(*map);
    }
};




IGA_NAMESPACE_CLOSE

#endif
