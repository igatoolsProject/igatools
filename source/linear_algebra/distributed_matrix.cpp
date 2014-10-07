//-+--------------------------------------------------------------------
// Igatools a general purpose Isogeometric analysis library.
// Copyright (C) 2012-2014  by the igatools authors (see authors.txt).
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

#include <igatools/linear_algebra/distributed_matrix.h>
#include <igatools/base/exceptions.h>

#include <boost/numeric/ublas/matrix_proxy.hpp>



using std::shared_ptr;

IGA_NAMESPACE_OPEN

namespace
{
DeclException2(ExcInvalidIndex,
               int, int,
               << "The entry with index <" << arg1 << ',' << arg2
               << "> does not exist.");

DeclException4(ExcAccessToNonLocalElement,
               int, int, int, int,
               << "You tried to access element (" << arg1
               << "/" << arg2 << ")"
               << " of a distributed matrix, but only rows "
               << arg3 << " through " << arg4
               << " are stored locally and can be accessed.");

DeclException0(ExcNotQuadratic);
};

#ifdef USE_TRILINOS


Matrix<LAPack::trilinos>::
Matrix(const SpaceManager &space_manager,CommPtr comm)
    :
    matrix_(trilinos_tools::build_matrix(
                trilinos_tools::build_graph(
                    space_manager,
                    trilinos_tools::build_row_map(space_manager,comm),
                    trilinos_tools::build_col_map(space_manager,comm))))
{
    /*
    matrix_.reset(new MatrixImpl(
                trilinos_tools::build_graph(
                        space_manager,
                        trilinos_tools::build_row_map(space_manager,comm),
                        trilinos_tools::build_col_map(space_manager,comm))));
    matrix_->setAllToScalar(0.0);
    //*/
}

#if 0
Matrix<LAPack::trilinos>::
Matrix(const SparsityPattern &sparsity_pattern,Teuchos::RCP<const Teuchos::Comm<int>> comm)
    :
    comm_(comm)
{
    //-------------------------------------------------------------------------------------
    const auto row_dofs = sparsity_pattern.get_row_dofs();
//  std::sort(row_dofs.begin(),row_dofs.end());

    const auto col_dofs = sparsity_pattern.get_col_dofs();
//  std::sort(col_dofs.begin(),col_dofs.end());

    //-------------------------------------------------------------------------------------


    column_space_map_.reset(new DofsMap(col_dofs.size(),col_dofs,0,comm_));

    row_space_map_.reset(new DofsMap(row_dofs.size(),row_dofs,0,comm_));

    using LongUInt = long unsigned int;
    Teuchos::ArrayRCP<const LongUInt> n_dofs_per_row =
        Teuchos::arcp(
            Teuchos::RCP<const std::vector<LongUInt> >(
                new vector<LongUInt>(sparsity_pattern.get_num_dof_connections()))) ;
    //-------------------------------------------------------------------------------------



    //-------------------------------------------------------------------------------------
    // allocating the entries in the graph corresponding to the sparsitiy pattern
    // (this step is required by the Tpetra matrix, in order to use sumIntoGlobalValues()

    graph_.reset(new Graph(row_space_map_,column_space_map_,n_dofs_per_row,Tpetra::StaticProfile));
    for (const auto &row : sparsity_pattern)
    {
        const Index row_id = row.first ;
        const auto &cols_id = row.second;

        auto cols_id_vec = vector<Index>(cols_id.begin(),cols_id.end());

        auto cols_id_view = Teuchos::ArrayView<const GO>(std::move(cols_id_vec));

        graph_->insertGlobalIndices(row_id,cols_id_view);
    }
    graph_->fillComplete(column_space_map_,row_space_map_);
    //-------------------------------------------------------------------------------------

    matrix_.reset(new MatrixImpl(graph_));
    matrix_->setAllToScalar(0.0);
    //-------------------------------------------------------------------------------------
};



shared_ptr<Matrix<LAPack::trilinos> >
Matrix<LAPack::trilinos>::
create(const SparsityPattern &sparsity_pattern)
{
    return std::make_shared<Matrix>(Matrix(sparsity_pattern));
}
#endif


shared_ptr<Matrix<LAPack::trilinos> >
Matrix<LAPack::trilinos>::
create(const SpaceManager &space_manager)
{
    return std::make_shared<Matrix>(Matrix(space_manager));
}



void
Matrix<LAPack::trilinos>::
clear()
{
    matrix_->resumeFill();
    matrix_->setAllToScalar(0.);
}



void
Matrix<LAPack::trilinos>::
add_entry(const Index row_id, const Index column_id, const Real value)
{
    Teuchos::Array<Index> columns_id(1,column_id) ;
    Teuchos::Array<Real> values(1,value) ;

    matrix_->sumIntoGlobalValues(row_id,columns_id,values);
};



void
Matrix<LAPack::trilinos>::
add_block(
    const vector<Index> &rows_id,
    const vector<Index> &cols_id,
    const DenseMatrix &local_matrix)
{
    Assert(!rows_id.empty(), ExcEmptyObject()) ;
    Assert(!cols_id.empty(), ExcEmptyObject()) ;

    const Index n_rows = rows_id.size();
    const Index n_cols = cols_id.size();

    Assert(n_rows == Index(local_matrix.size1()),
           ExcDimensionMismatch(n_rows,Index(local_matrix.size1()))) ;
    Assert(n_cols == Index(local_matrix.size2()),
           ExcDimensionMismatch(n_cols,Index(local_matrix.size2()))) ;

    vector<Real> row_values(n_cols) ;
    for (int i = 0 ; i < n_rows ; ++i)
    {
        const auto row_local_matrix = local_matrix.get_row(i) ;
        for (int j = 0 ; j < n_cols ; ++j)
            row_values[j] = row_local_matrix(j) ;

        matrix_->sumIntoGlobalValues(rows_id[i],cols_id,row_values);
    }
};



void
Matrix<LAPack::trilinos>::
fill_complete()
{
    const auto graph = matrix_->getGraph();
    matrix_->fillComplete(graph->getDomainMap(),graph->getRangeMap());
};

void
Matrix<LAPack::trilinos>::
resume_fill()
{
    matrix_->resumeFill();
};


auto
Matrix<LAPack::trilinos>::
get_trilinos_matrix() -> Teuchos::RCP<MatrixImpl>
{
    return matrix_;
};

auto
Matrix<LAPack::trilinos>::
get_trilinos_matrix() const -> Teuchos::RCP<const MatrixImpl>
{
    return matrix_;
};


Real
Matrix<LAPack::trilinos>::
operator()(const Index row, const Index col) const
{
    const auto graph = matrix_->getGraph();

    const auto row_map = graph->getRowMap();
    const auto col_map = graph->getColMap();

    const auto local_row = row_map->getLocalElement(row);
    const auto local_col = col_map->getLocalElement(col);
//    std::cout << "Global=("<<row<<","<<col<<")   Local=("<<local_row<<","<<local_col<<")"<<std::endl;

    // If the data is not on the present processor, we throw an exception.
    // This is one of the two tiny differences to the el(i,j) call,
    //which does not throw any assertions, but returns 0.0.
    Assert(local_row != Teuchos::OrdinalTraits<Index>::invalid(),
           ExcAccessToNonLocalElement(row,col,
                                      row_map->getMinGlobalIndex(),
                                      row_map->getMaxGlobalIndex()+1));

    Teuchos::ArrayView<const Index> local_col_ids ;
    Teuchos::ArrayView<const Real> values ;
    matrix_->getLocalRowView(local_row,local_col_ids,values) ;

//    std::cout << local_col_ids <<std::endl ;
//    std::cout << "local_col=" << local_col << std::endl;
    // Search the index where we look for the value, and then finally get it.
    const Index *col_find = std::find(local_col_ids.begin(),local_col_ids.end(),local_col);


    // This is actually the only difference to the el(i,j) function,
    // which means that we throw an exception in this case instead of just
    // returning zero for an element that is not present in the sparsity pattern.
    Assert(col_find != local_col_ids.end(), ExcInvalidIndex(row,col));

    const Index id_find = static_cast<Index>(col_find-local_col_ids.begin());
//    std::cout << "id_find="<<id_find<<std::endl;

    return values[id_find];
}


void
Matrix<LAPack::trilinos>::
clear_row(const Index row)
{
    const auto graph = matrix_->getGraph();

    const auto row_map = graph->getRowMap();
    const auto local_row = row_map->getLocalElement(row);

    // Only do this on the rows owned locally on this processor.
    if (local_row != Teuchos::OrdinalTraits<Index>::invalid())
    {
        auto n_entries_row = graph->getNumEntriesInLocalRow(local_row);

        vector<Index> local_col_ids(n_entries_row);

        graph->getLocalRowCopy(local_row,local_col_ids,n_entries_row);
        Assert(n_entries_row == graph->getNumEntriesInLocalRow(local_row),
               ExcDimensionMismatch(n_entries_row,graph->getNumEntriesInLocalRow(local_row)));

        vector<Real> zeros(n_entries_row,0.0);
        matrix_->replaceLocalValues(local_row,local_col_ids,zeros);
    }
}



void
Matrix<LAPack::trilinos>::
print(LogStream &out) const
{
    using std::endl;

    const auto n_rows = this->get_num_rows();


    out << "-----------------------------" << endl;
    // Commented as different trilinos version show different information here
    //   out << *matrix_ ;
    out << "Num. rows    = " << n_rows << endl;
    out << "Num. cols    = " << this->get_num_columns() << endl;
    out << "Num. entries = " << this->get_num_entries() << endl;
    out << endl;
    out << "Row Index        Col Index        Value" << endl;
    for (Index row_id = 0 ; row_id < n_rows ; ++row_id)
    {
        auto n_entries_row = matrix_->getNumEntriesInGlobalRow(row_id);

        vector<Index> columns_id(n_entries_row);

        vector<Real> values(n_entries_row) ;

        matrix_->getGlobalRowCopy(row_id,columns_id,values,n_entries_row);

        for (const auto col_id : columns_id)
            out << row_id << "       "
                << col_id << "        "
                << (*this)(row_id,col_id) << endl;
    }
    out << "-----------------------------" << endl;

}





auto
Matrix<LAPack::trilinos>::
get_num_rows() const -> Index
{
    return matrix_->getGlobalNumRows() ;
}

auto
Matrix<LAPack::trilinos>::
get_num_columns() const -> Index
{
    return matrix_->getGlobalNumCols() ;
}

auto
Matrix<LAPack::trilinos>::
get_num_entries() const -> Index
{
    return matrix_->getGlobalNumEntries() ;
}


void
Matrix<LAPack::trilinos>::
multiply_by_right_vector(const vector_t &x,vector_t &y,const Real alpha,const Real beta) const
{
    matrix_->apply(*x.get_trilinos_vector(),
                   *y.get_trilinos_vector(),
                   Teuchos::NO_TRANS,alpha,beta);
}

void
Matrix<LAPack::trilinos>::
multiply_by_left_vector(const vector_t &x,vector_t &y,const Real alpha,const Real beta) const
{
    matrix_->apply(*x.get_trilinos_vector(),
                   *y.get_trilinos_vector(),
                   Teuchos::TRANS,alpha,beta);
}

#endif //#ifdef USE_TRILINOS








#ifdef USE_PETSC


Matrix<LAPack::petsc>::
Matrix(const SparsityPattern &sparsity_pattern)
{
    init(sparsity_pattern) ;
};


void
Matrix<LAPack::petsc>::
init(const SparsityPattern &sparsity_pattern)
{
    PetscErrorCode ierr;
    comm_ = PETSC_COMM_WORLD;

    int n_rows = sparsity_pattern.get_num_row_dofs();
    int n_cols = sparsity_pattern.get_num_col_dofs();
    vector<PetscInt> nnz;

    ierr = MatCreate(comm_, &matrix_);  // CHKERRQ(ierr);
    ierr = MatCreate(comm_, &matrix_);  // CHKERRQ(ierr);
    ierr = MatSetSizes(matrix_,
                       PETSC_DECIDE,PETSC_DECIDE,
                       n_rows, n_cols);
    ierr = MatSetType(matrix_,MATAIJ);

    auto row     = sparsity_pattern.cbegin() ;
    auto row_end = sparsity_pattern.cend() ;
    for (; row != row_end ; ++row)
    {
        const vector<Index> columns_id(row->second.begin(),row->second.end()) ;
        nnz.push_back(columns_id.size());
    }

    ierr = MatSeqAIJSetPreallocation(matrix_, 0, nnz.data());
    ierr = MatSetOption(matrix_, MAT_KEEP_NONZERO_PATTERN, PETSC_TRUE); //CHKERRQ(ierr);
    ierr = MatSetUp(matrix_); //CHKERRQ(ierr);
    ierr = MatZeroEntries(matrix_); //CHKERRQ(ierr);

};


shared_ptr<Matrix<LAPack::petsc> >
Matrix<LAPack::petsc>::
create(const SparsityPattern &sparsity_pattern)
{
    return std::make_shared<self_t>(self_t(sparsity_pattern));
}

void
Matrix<LAPack::petsc>::
add_entry(const Index row_id, const Index column_id, const Real value)
{
    PetscErrorCode ierr;
    ierr = MatSetValue(matrix_,row_id,column_id,value,ADD_VALUES); //CHKERRQ(ierr);
};



void
Matrix<LAPack::petsc>::
add_block(
    const vector<Index> &rows_id,
    const vector<Index> &cols_id,
    const DenseMatrix &local_matrix)
{
    Assert(!rows_id.empty(), ExcEmptyObject()) ;
    Assert(!cols_id.empty(), ExcEmptyObject()) ;

    const Index n_rows = rows_id.size();
    const Index n_cols = cols_id.size();

    Assert(n_rows == Index(local_matrix.size1()),
           ExcDimensionMismatch(n_rows,Index(local_matrix.size1()))) ;
    Assert(n_cols == Index(local_matrix.size2()),
           ExcDimensionMismatch(n_cols,Index(local_matrix.size2()))) ;

    vector<PetscScalar> row_values ;
    for (int i = 0 ; i < n_rows ; ++i)
    {
        const auto row_local_matrix = local_matrix.get_row(i) ;
        for (int j = 0 ; j < n_cols ; ++j)
//            row_values.push_back(row_local_matrix(j)) ;
            row_values.push_back(PetscScalar(local_matrix(i,j))) ;
    }

    PetscErrorCode ierr;
    ierr = MatSetValues(matrix_,n_rows,rows_id.data(),n_cols,cols_id.data(),
                        row_values.data(),ADD_VALUES); //CHKERRQ(ierr);

};



void
Matrix<LAPack::petsc>::
fill_complete()
{
    PetscErrorCode ierr;
    ierr = MatAssemblyBegin(matrix_, MAT_FINAL_ASSEMBLY); // CHKERRQ(ierr);
    ierr = MatAssemblyEnd(matrix_, MAT_FINAL_ASSEMBLY); // CHKERRQ(ierr);
};

void
Matrix<LAPack::petsc>::
resume_fill()
{
    /* todo:I think this is not necessary in PETSc. Must check in parallel computations */
};


auto
Matrix<LAPack::petsc>::
get_petsc_matrix() -> Mat
{
    return matrix_;
};

auto
Matrix<LAPack::petsc>::
get_petsc_matrix() const -> Mat
{
    return matrix_;
};


Real
Matrix<LAPack::petsc>::
operator()(const Index row, const Index col) const
{
    PetscErrorCode ierr;
    PetscScalar value;

    /* todo: the assert should work with a distributed matrix */
    AssertIndexRange(row, this->get_num_rows());
    AssertIndexRange(col, this->get_num_columns());

    /* todo: the function returns zero if the entry is not in the sparsity
     pattern, instead it should give an error. I have not found a function in
     Petsc to check whether the entry is in the spartsity pattern or not. */

    ierr = MatGetValues(matrix_,1,&row,1,&col,&value);

    return value;
}


void
Matrix<LAPack::petsc>::
clear_row(const Index row)
{
    PetscErrorCode ierr;

    Assert(false,ExcNotImplemented());
    AssertThrow(false,ExcNotImplemented());
    /*
        const auto graph = matrix_->getGraph();

        const auto row_map = graph->getRowMap();
        const auto local_row = row_map->getLocalElement(row);

        // Only do this on the rows owned locally on this processor.
        if (local_row != Teuchos::OrdinalTraits<Index>::invalid())
        {
            auto n_entries_row = graph->getNumEntriesInLocalRow(local_row);

            vector<Index> local_col_ids(n_entries_row);

            graph->getLocalRowCopy(local_row,local_col_ids,n_entries_row);
            Assert(n_entries_row == graph->getNumEntriesInLocalRow(local_row),
                   ExcDimensionMismatch(n_entries_row,graph->getNumEntriesInLocalRow(local_row)));

            vector<Real> zeros(n_entries_row,0.0);
            matrix_->replaceLocalValues(local_row,local_col_ids,zeros);
        }
     //*/
}



void
Matrix<LAPack::petsc>::
print(LogStream &out) const
{
    Assert(false,ExcNotImplemented());
    AssertThrow(false,ExcNotImplemented());
    /*
        using std::endl;

        const auto n_global_rows = this->get_num_rows();

        out << "-----------------------------" << endl;
        // Commented as different petsc version show different information here
        //   out << *matrix_ ;

        out << "Row Index        Col Index        Value" << endl;
        for (Index row_id = 0 ; row_id < n_global_rows ; ++row_id)
        {
            auto n_entries_row = matrix_->getNumEntriesInGlobalRow(row_id);

            vector<Index> columns_id(n_entries_row);

            vector<Real> values(n_entries_row) ;

            matrix_->getGlobalRowCopy(row_id,columns_id,values,n_entries_row);

            for (uint row_entry = 0 ; row_entry < n_entries_row ; ++row_entry)
                out << row_id << "       " <<
                    columns_id[row_entry] << "        " <<
                    values[row_entry] << endl;
        }
        out << "-----------------------------" << endl;
    //*/
}





auto
Matrix<LAPack::petsc>::
get_num_rows() const -> Index
{
    PetscErrorCode ierr;
    PetscInt n_rows, n_cols;

    ierr = MatGetSize(matrix_, &n_rows, &n_cols);

    return n_rows;
}

auto
Matrix<LAPack::petsc>::
get_num_columns() const -> Index
{
    PetscErrorCode ierr;
    PetscInt n_rows, n_cols;

    ierr = MatGetSize(matrix_, &n_rows, &n_cols);

    return n_cols;
}


void
Matrix<LAPack::petsc>::
multiply_by_right_vector(const vector_t &x,vector_t &y,const Real alpha,const Real beta) const
{
    Assert(false,ExcNotImplemented());
    AssertThrow(false,ExcNotImplemented());
    /*
        matrix_->apply(*x.get_petsc_vector(),
                       *y.get_petsc_vector(),
                       Teuchos::NO_TRANS,alpha,beta);
                       //*/
}

void
Matrix<LAPack::petsc>::
multiply_by_left_vector(const vector_t &x,vector_t &y,const Real alpha,const Real beta) const
{
    Assert(false,ExcNotImplemented());
    AssertThrow(false,ExcNotImplemented());
    /*
        matrix_->apply(*x.get_petsc_vector(),
                       *y.get_petsc_vector(),
                       Teuchos::TRANS,alpha,beta);
                       //*/
}

#endif //#ifdef USE_PETSC



IGA_NAMESPACE_CLOSE



