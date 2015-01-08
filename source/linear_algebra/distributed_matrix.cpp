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

template<TrilinosImpl trilinos_impl>
MatrixTrilinos<trilinos_impl>::
MatrixTrilinos(const SpaceManager &space_manager,CommPtr comm)
    :
    matrix_(Tools::build_matrix(
                Tools::build_graph(
                    space_manager,
                    Tools::build_row_map(space_manager,comm),
                    Tools::build_col_map(space_manager,comm))))
{}


template<TrilinosImpl trilinos_impl>
auto
MatrixTrilinos<trilinos_impl>::
get_trilinos_matrix() -> Teuchos::RCP<WrappedMatrix>
{
    return matrix_;
};

template<TrilinosImpl trilinos_impl>
auto
MatrixTrilinos<trilinos_impl>::
get_trilinos_matrix() const -> Teuchos::RCP<const WrappedMatrix>
{
    return matrix_;
};


Matrix<LAPack::trilinos_tpetra>::
Matrix(const SpaceManager &space_manager,CommPtr comm)
    :
    MatrixTrilinos<TrilinosImpl::tpetra>(space_manager,comm)
{}



shared_ptr<Matrix<LAPack::trilinos_tpetra> >
Matrix<LAPack::trilinos_tpetra>::
create(const SpaceManager &space_manager)
{
    using Mat = Matrix<LAPack::trilinos_tpetra>;
    return std::make_shared<Mat>(Mat(space_manager));
}



void
Matrix<LAPack::trilinos_tpetra>::
clear()
{
    matrix_->resumeFill();
    matrix_->setAllToScalar(0.);
}



void
Matrix<LAPack::trilinos_tpetra>::
add_entry(const Index row_id, const Index column_id, const Real value)
{
    Teuchos::Array<Index> columns_id(1,column_id) ;
    Teuchos::Array<Real> values(1,value) ;

    matrix_->sumIntoGlobalValues(row_id,columns_id,values);
};



void
Matrix<LAPack::trilinos_tpetra>::
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

//    vector<Real> row_values(n_cols) ;
    for (int i = 0 ; i < n_rows ; ++i)
    {
//        const auto row_local_matrix = local_matrix.get_row(i) ;
//        for (int j = 0 ; j < n_cols ; ++j)
//            row_values[j] = row_local_matrix(j) ;

        auto row_values = Teuchos::ArrayView<const Real>(&local_matrix(i,0),n_cols);

        matrix_->sumIntoGlobalValues(rows_id[i],cols_id,row_values);
    }
};



void
Matrix<LAPack::trilinos_tpetra>::
fill_complete()
{
    const auto graph = matrix_->getGraph();
    matrix_->fillComplete(graph->getDomainMap(),graph->getRangeMap());
};

void
Matrix<LAPack::trilinos_tpetra>::
resume_fill()
{
    matrix_->resumeFill();
};


Real
Matrix<LAPack::trilinos_tpetra>::
operator()(const Index row, const Index col) const
{
    const auto graph = matrix_->getGraph();

    const auto row_map = graph->getRowMap();
    const auto col_map = graph->getColMap();

    const auto local_row = row_map->getLocalElement(row);
    const auto local_col = col_map->getLocalElement(col);

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

    // Search the index where we look for the value, and then finally get it.
    const auto col_find = std::find(local_col_ids.begin(),local_col_ids.end(),local_col);


    // This is actually the only difference to the el(i,j) function,
    // which means that we throw an exception in this case instead of just
    // returning zero for an element that is not present in the sparsity pattern.
    Assert(col_find != local_col_ids.end(), ExcInvalidIndex(row,col));

    const Index id_find = static_cast<Index>(col_find-local_col_ids.begin());

    return values[id_find];
}


void
Matrix<LAPack::trilinos_tpetra>::
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
Matrix<LAPack::trilinos_tpetra>::
print_info(LogStream &out) const
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
    for (Index local_row = 0 ; local_row < n_rows ; ++local_row)
    {
        const auto graph = matrix_->getGraph();

        const auto row_map = graph->getRowMap();

        const auto global_row = row_map->getGlobalElement(local_row);

        auto n_entries_row = matrix_->getNumEntriesInGlobalRow(global_row);

        vector<Index> columns_id(n_entries_row);

        vector<Real> values(n_entries_row) ;

        matrix_->getGlobalRowCopy(global_row,columns_id,values,n_entries_row);

        for (const auto global_col : columns_id)
            out << global_row << "       "
                << global_col << "        "
                << (*this)(global_row,global_col) << endl;
    }
    out << "-----------------------------" << endl;

}





auto
Matrix<LAPack::trilinos_tpetra>::
get_num_rows() const -> Index
{
    return matrix_->getGlobalNumRows() ;
}

auto
Matrix<LAPack::trilinos_tpetra>::
get_num_columns() const -> Index
{
    return matrix_->getGlobalNumCols() ;
}

auto
Matrix<LAPack::trilinos_tpetra>::
get_num_entries() const -> Index
{
    return matrix_->getGlobalNumEntries() ;
}


void
Matrix<LAPack::trilinos_tpetra>::
multiply_by_right_vector(const vector_t &x,vector_t &y,const Real alpha,const Real beta) const
{
    matrix_->apply(*x.get_trilinos_vector(),
                   *y.get_trilinos_vector(),
                   Teuchos::NO_TRANS,alpha,beta);
}

void
Matrix<LAPack::trilinos_tpetra>::
multiply_by_left_vector(const vector_t &x,vector_t &y,const Real alpha,const Real beta) const
{
    matrix_->apply(*x.get_trilinos_vector(),
                   *y.get_trilinos_vector(),
                   Teuchos::TRANS,alpha,beta);
}







Matrix<LAPack::trilinos_epetra>::
Matrix(const SpaceManager &space_manager,CommPtr comm)
    :
    MatrixTrilinos<TrilinosImpl::epetra>(space_manager,comm)
{}

shared_ptr<Matrix<LAPack::trilinos_epetra> >
Matrix<LAPack::trilinos_epetra>::
create(const SpaceManager &space_manager)
{
    using Mat = Matrix<LAPack::trilinos_epetra>;
    return std::make_shared<Mat>(Mat(space_manager));
}


void
Matrix<LAPack::trilinos_epetra>::
add_entry(const Index row_id, const Index column_id, const Real value)
{
    matrix_->SumIntoGlobalValues(row_id,1,&value,&column_id);
};


void
Matrix<LAPack::trilinos_epetra>::
fill_complete()
{
//    const auto &graph = matrix_->Graph();
//    matrix_->FillComplete(graph.DomainMap(),graph.RangeMap());
    matrix_->FillComplete(matrix_->ColMap(),matrix_->RowMap());
};

void
Matrix<LAPack::trilinos_epetra>::
resume_fill()
{};

void
Matrix<LAPack::trilinos_epetra>::
clear()
{
    matrix_->PutScalar(0.);
}

auto
Matrix<LAPack::trilinos_epetra>::
get_num_rows() const -> Index
{
    return matrix_->NumGlobalRows() ;
}

auto
Matrix<LAPack::trilinos_epetra>::
get_num_columns() const -> Index
{
    return matrix_->NumGlobalCols() ;
}

auto
Matrix<LAPack::trilinos_epetra>::
get_num_entries() const -> Index
{
    return matrix_->GlobalMaxNumEntries() ;
}


void
Matrix<LAPack::trilinos_epetra>::
multiply_by_right_vector(const vector_t &x,vector_t &y,const Real alpha,const Real beta) const
{
    auto x_ptr = x.get_trilinos_vector();
    using TrilinosVec = typename vector_t::WrappedVector;

    auto A_times_X = Teuchos::rcp(new TrilinosVec(matrix_->RangeMap(),x_ptr->NumVectors()));
    matrix_->Multiply(false,*x_ptr,*A_times_X);

    y.get_trilinos_vector()->Update(alpha, *A_times_X, beta);
}

void
Matrix<LAPack::trilinos_epetra>::
multiply_by_left_vector(const vector_t &x,vector_t &y,const Real alpha,const Real beta) const
{
    auto x_ptr = x.get_trilinos_vector();
    auto A_times_X = *x_ptr;
    matrix_->Multiply(true,*x_ptr,A_times_X);

    y.get_trilinos_vector()->Update(alpha, A_times_X, beta);
}


void
Matrix<LAPack::trilinos_epetra>::
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

        matrix_->SumIntoGlobalValues(rows_id[i],n_cols,row_values.data(),cols_id.data());
    }
};




Real
Matrix<LAPack::trilinos_epetra>::
operator()(const Index row, const Index col) const
{
    const auto &graph = matrix_->Graph();

    const auto &row_map = graph.RowMap();
    const auto &col_map = graph.ColMap();

    const auto local_row = row_map.LID(row);
    const auto local_col = col_map.LID(col);

    // If the data is not on the present processor, we throw an exception.
    // This is one of the two tiny differences to the el(i,j) call,
    //which does not throw any assertions, but returns 0.0.
    Assert(local_row != Teuchos::OrdinalTraits<Index>::invalid(),
           ExcAccessToNonLocalElement(row,col,
                                      row_map.MinAllGID(),
                                      row_map.MaxAllGID()+1));

    Index n_cols = matrix_->NumMyEntries(local_row);

    Real *val_ptr;
    Index *id_ptr;
    matrix_->ExtractMyRowView(local_row, n_cols, val_ptr, id_ptr);
    Teuchos::ArrayView<const Index> local_col_ids(id_ptr,n_cols) ;
    Teuchos::ArrayView<const Real> values(val_ptr,n_cols) ;

    // Search the index where we look for the value, and then finally get it.
    const auto col_find = std::find(local_col_ids.begin(),local_col_ids.end(),local_col);


    // This is actually the only difference to the el(i,j) function,
    // which means that we throw an exception in this case instead of just
    // returning zero for an element that is not present in the sparsity pattern.
    Assert(col_find != local_col_ids.end(), ExcInvalidIndex(row,col));

    const Index id_find = static_cast<Index>(col_find-local_col_ids.begin());

    return values[id_find];
}


void
Matrix<LAPack::trilinos_epetra>::
clear_row(const Index row)
{
    const auto &graph = matrix_->Graph();

    const auto &row_map = graph.RowMap();
    const auto local_row = row_map.LID(row);

    // Only do this on the rows owned locally on this processor.
    if (local_row != Teuchos::OrdinalTraits<Index>::invalid())
    {
        auto n_entries_row = graph.NumMyIndices(local_row);

        vector<Index> local_col_ids(n_entries_row);

        graph.ExtractMyRowCopy(local_row,n_entries_row,n_entries_row,local_col_ids.data());

        vector<Real> zeros(n_entries_row,0.0);
        matrix_->ReplaceMyValues(local_row,n_entries_row,zeros.data(),local_col_ids.data());
    }
}


void
Matrix<LAPack::trilinos_epetra>::
print_info(LogStream &out) const
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

    for (Index local_row = 0 ; local_row < n_rows ; ++local_row)
    {
        const auto &graph = matrix_->Graph();

        const auto &row_map = graph.RowMap();

        const auto global_row = row_map.GID(local_row);

        Index n_entries_row = matrix_->NumGlobalEntries(global_row);
        vector<Real> values(n_entries_row);
        vector<Index> columns_id(n_entries_row);

        Index nnz = 0;
        matrix_->ExtractGlobalRowCopy(global_row,n_entries_row,nnz,values.data(),columns_id.data());

        for (Index i = 0 ; i < n_entries_row ; ++i)
            out << global_row << "       "
                << columns_id[i] << "        "
                << values[i] << endl;
    }
    out << "-----------------------------" << endl;

}



template class MatrixTrilinos<TrilinosImpl::tpetra>;
template class MatrixTrilinos<TrilinosImpl::epetra>;

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
print_info(LogStream &out) const
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



