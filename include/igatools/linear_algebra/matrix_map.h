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

#include <Teuchos_RCP.hpp>
#include <Epetra_SerialComm.h>
#include <Epetra_Map.h>

#include <Epetra_CrsMatrix.h>

//#include <Epetra_Vector.h>



IGA_NAMESPACE_OPEN

template<LAPack lapack> class MatrixGraph;

template <>
class MatrixGraph<LAPack::trilinos_epetra>
{
	using Comm = Epetra_Comm;
	using CommPtr = Teuchos::RCP<const Comm>;

	using Map = Epetra_Map;
	using MapPtr = Teuchos::RCP<Map>;

	using Graph = Epetra_CrsGraph;
	using GraphPtr = Teuchos::RCP<Graph>;

	using Matrix = Epetra_CrsMatrix;
	using MatrixPtr = Teuchos::RCP<Matrix>;

	using Vector = Epetra_MultiVector;
	using VectorPtr = Teuchos::RCP<Vector>;

public:
	template<class RowSpacePtr, class ColSpacePtr>
	MatrixGraph(const RowSpacePtr row_space, const std::string &row_property,
			    const ColSpacePtr col_space, const std::string &col_property,
				CommPtr comm = Teuchos::rcp(new Epetra_SerialComm()))
	{
		const auto row_dof_dist = row_space->get_dof_distribution();
		const auto row_dofs = row_dof_dist->get_dofs_id_same_property(row_property);
		const auto col_dof_dist = col_space->get_dof_distribution();
		const auto col_dofs = col_dof_dist->get_dofs_id_same_property(col_property);
		//TODO (pauletti, Mar 28, 2015): this is double copy of data
		const vector<Index> row_dofs_vec(row_dofs.begin(), row_dofs.end());
		row_map_ = Teuchos::rcp(
				new Epetra_Map(-1, row_dofs_vec.size(), row_dofs_vec.data(),0,*comm));
		const vector<Index> col_dofs_vec(col_dofs.begin(), col_dofs.end());
		col_map_ = Teuchos::rcp(
				new Epetra_Map(-1, col_dofs_vec.size(), col_dofs_vec.data(),0,*comm));

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

		LogStream out;
		n_dofs_per_row.print_info(out);
		loc_dofs.print_info(out);

		const bool is_static_profile = true;
		auto graph_ = Teuchos::rcp(
				new Graph(Epetra_DataAccess::Copy,
						*row_map_, *col_map_,
						n_dofs_per_row.data(),
						is_static_profile));

		Index j=0;
		for (auto &dofs : loc_dofs)
		{
			const Index row_id = row_map_->GID(j);
			graph_->InsertGlobalIndices(row_id, dofs.size(), dofs.data());
			++j;
		}

		int res = graph_->FillComplete(*col_map_,*row_map_);
		Assert(res==0, ExcMessage(" "));
		out << graph_->NumMyCols () << std::endl;
		graph_->Print(out.get_console());

	}


	MatrixPtr create_matrix() const
	{
		auto matrix = Teuchos::rcp(new Matrix(Epetra_DataAccess::Copy,*graph_));
		matrix->PutScalar(0.0);


	    return matrix;
	}

	void print_info(LogStream &out) const
	{
		out.begin_item("Row map:");
		row_map_->Print(out.get_file_stream());
		out.end_item();

		out.begin_item("Col map:");
		col_map_->Print(out.get_file_stream());
		out.end_item();
	//	out << graph_->NumMyCols () << std::endl;
//		out.begin_item("Graph:");
//		graph_->Print(out.get_file_stream());
//		out.end_item();
	}



private:
	GraphPtr graph_;
	MapPtr row_map_;
	MapPtr col_map_;
};


IGA_NAMESPACE_CLOSE

#endif
