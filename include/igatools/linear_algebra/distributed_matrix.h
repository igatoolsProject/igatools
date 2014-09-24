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

#ifndef DISTRIBUTED_MATRIX_H_
#define DISTRIBUTED_MATRIX_H_

#include <igatools/base/config.h>
#include <igatools/base/logstream.h>
#include <igatools/linear_algebra/dense_matrix.h>
#include <igatools/linear_algebra/sparsity_pattern.h>
#include <igatools/linear_algebra/distributed_vector.h>

#ifdef USE_TRILINOS
#include <Tpetra_CrsMatrix.hpp>
#endif

#ifdef USE_PETSC
#include <petscmat.h>
#endif

IGA_NAMESPACE_OPEN

template <LAPack la_pack>
class Matrix;

#ifdef USE_TRILINOS
/**
 * @todo Missing documentation
 *
 *
 */
template <>
class Matrix<LAPack::trilinos>
{
private:
    /** Type alias for the local ordinal types (i.e. the types for the local indices). */
    using LO = Index;

    /** Type alias for the global ordinal types (i.e. the types for the global indices). */
    using GO = Index;
#if 0
    /**
     The Kokkos "Node" type describes the type of shared-memory
     parallelism that Tpetra will use _within_ an MPI process.  The
     available Node types depend on Trilinos' build options and the
     availability of certain third-party libraries.  Here are a few
     examples:

     Kokkos::SerialNode: No parallelism

     Kokkos::TPINode: Uses a custom Pthreads wrapper

     Kokkos::TBBNode: Uses Intel's Threading Building Blocks

     Kokkos::ThrustNode: Uses Thrust, a C++ CUDA wrapper,
     for GPU parallelism.

     Using a GPU-oriented Node means that Tpetra objects that store a
     lot of data (vectors and sparse matrices, for example) will store
     that data on the GPU, and operate on it there whenever possible.

     Kokkos::DefaultNode gives you a default Node type.  It may be
     different, depending on Trilinos' build options.  Currently, for
     example, building Trilinos with Pthreads enabled gives you
     Kokkos::TPINode by default.  That means your default Node is a
     parallel node!
    */
    using Node = typename Kokkos::SerialNode;
#endif

    /** Type alias for the dofs map across the processors */
    using DofsMap = typename Tpetra::Map<LO,GO>;


    /** Type alias for the connecitivty graph */
    using Graph = typename Tpetra::CrsGraph<LO,GO>;



public:
    /** Typedef for the matrix type */
    using WrappedMatrixType = Tpetra::CrsMatrix<Real,LO,GO> ;

    using self_t = Matrix<LAPack::trilinos>;

    using vector_t = Vector<LAPack::trilinos>;

public:
    /**@name Constructor and destructor */
    ///@{
    /** Default constructor */
    Matrix() = delete;

    /**
     * Construct a distributed matrix with the dof distribution for its rows and column
     * specified by the SparsityPattern @p sparsity_pattern.
     */
    Matrix(const SparsityPattern &sparsity_pattern,
           Teuchos::RCP<const Teuchos::Comm<int>> comm = Teuchos::createSerialComm<int>());

    /**
     * Copy constructor. Not allowed to be used.
     */
    Matrix(const self_t &matrix) = delete;

    /**
     * Move constructor.
     */
    Matrix(self_t &&matrix) = default;

    /** Destructor */
    ~Matrix() = default;
    ///@}

    /**
     * @name Function for creating Matrix objects
     * @note These methods implement the "create idiom" and return
     * a Matrix object wrapped by a std::shared_ptr
     */
    ///@{
    /**
     * Create a distributed matrix with the dof dostribution for its rows and column
     * specified by the SparsityPattern @p sparsity_pattern.
     */
    static std::shared_ptr<self_t> create(const SparsityPattern &sparsity_pattern);
    ///@}

    /** @name Assignment operators */
    ///@{

    /**
     * Copy assignment operator. Not allowed to be used.
     */
    Matrix &operator=(const self_t &matrix) = delete;

    /**
     * Move assignment operator. Not allowed to be used.
     */
    Matrix &operator=(self_t &&matrix) = delete;

    ///@}

    /** @name Methods for getting and/or modifying the matrix entries */
    ///@{

    /**
     * Add the value @p input to the matrix entry (i,j).
     * @note @p i and @p j are the global indices of the entry (i,j).
     */
    void add_entry(const Index i,const  Index j,const Real input);

    /**
     * \brief This function add the local matrix values to the global matrix, the local-to-global ids
     * are passed as input argument.
     * \param[in] row_glob_ids The vector containing the row global ids
     * associated to the local matrix entries
     * \param[in] col_glob_ids The vector containing the colume global ids
     * associated to the local matrix entries
     * \param[in] local_matrix The local matrix that must be added to the global matrix.
     * \note The number of rows and columns of the local matrix must be equal to the dimension of the
     * local_to_global vector, otherwise an exception will be raised.
     */
    void add_block(
        const vector<Index> &row_glob_ids,
        const vector<Index> &cln_glob_ids,
        const DenseMatrix &local_matrix) ;

    /**
     * Given two global indices @p i and @p j, return the (i,j) entry of the matrix.
     */
    Real operator()(const Index i, const Index j) const;


    /** Set the row specified by the global id @p row to zero. */
    void clear_row(const Index row);
    ///@}

    /** @name Methods for retrieving the Trilinos objects wrapped by this class */
    ///@{
    /**
     * Return the Trilinos RCP (Reference-Counted-Pointer, i.e. a smart pointer) wrapping the
     * concrete Trilinos distributed matrix. Const version.
     */
    Teuchos::RCP<const WrappedMatrixType> get_trilinos_matrix() const ;

    /**
     * Return the Trilinos RCP (Reference-Counted-Pointer, i.e. a smart pointer) wrapping the
     * concrete Trilinos distributed matrix. Non-const version.
     */
    Teuchos::RCP<WrappedMatrixType> get_trilinos_matrix() ;
    ///@}

    /** @name Methods for retrieving or printing the matrix informations */
    ///@{
    void print(LogStream &out) const ;

    /** Return the number of global rows of the matrix */
    Index get_num_rows() const ;

    /** Return the number of global columns of the matrix */
    Index get_num_columns() const ;
    ///@}

    /** @name Matrix-by-vector multiplication */
    ///@{
    /**
     * Compute a sparse matrix-MultiVector multiply.
     *
     * This method computes y := beta * y + alpha * A * x.
     * If beta == 0, this operation will enjoy overwrite semantics: y's entries will be ignored,
     * and y will be overwritten with the result of the multiplication, even if it contains NaN
     * (not-a-number) floating-point entries.
     */
    void multiply_by_right_vector(const vector_t &x,vector_t &y,const Real alpha,const Real beta) const;

    /**
     * Compute a sparse matrix-MultiVector multiply.
     *
     * This method computes y^T := beta * y^T + alpha * x^T * A.
     * If beta == 0, this operation will enjoy overwrite semantics: y's entries will be ignored,
     * and y will be overwritten with the result of the multiplication, even if it contains NaN
     * (not-a-number) floating-point entries.
     */
    void multiply_by_left_vector(const vector_t &x,vector_t &y,const Real alpha,const Real beta) const;
    ///@}

    /** @name Methods for changing the filling state. */
    ///@{
    /**
     * Signal that data entry is complete, specifying domain and range maps.
     *
     * Off-node indices are distributed, indices are sorted,
     * redundant indices are eliminated, and global indices are transformed to local indices.
     */
    void fill_complete();

    /**
     * Resume operations that may change the values or structure of the matrix.
     *
     * This method must be called as a collective operation.
     * Calling fill_complete() "freezes" both the values and the structure of the matrix.
     * If you want to modify the matrix again, you must first call resume_fill().
     * You then may not call resume_fill() again on that matrix until you first call fill_complete().
     * You may make sequences of fill_complete(), resume_fill() calls as many times as you wish.
     */
    void resume_fill();
    ///@}

private:
    /** The real Trilinos::TPetra matrix */
    Teuchos::RCP<WrappedMatrixType> matrix_ ;


//    Teuchos::RCP<DofsMap> all_dofs_map_;
    Teuchos::RCP<DofsMap> row_space_map_;
    Teuchos::RCP<DofsMap> column_space_map_;

    Teuchos::RCP<Graph> graph_;

    Teuchos::RCP<const Teuchos::Comm<int>> comm_;

    void init(const SparsityPattern &sparsity_pattern);
};
#endif // #ifdef USE_TRILINOS



#ifdef USE_PETSC
/**
 * @todo Missing documentation
 *
 *
 */
template <>
class Matrix<LAPack::petsc>
{
public:
    /** Typedef for the matrix type */
    using self_t = Matrix<LAPack::petsc>;
    using vector_t = Vector<LAPack::petsc>;

public:
    /**@name Constructor and destructor */
    ///@{
    /** Default constructor */
    Matrix() = delete;

    /**
     * Construct a distributed matrix with the dof distribution for its rows and column
     * specified by the SparsityPattern @p sparsity_pattern.
     */
    Matrix(const SparsityPattern &sparsity_pattern);

    /**
     * Copy constructor. Not allowed to be used.
     */
    Matrix(const self_t &matrix) = delete;

    /**
     * Move constructor.
     */
    Matrix(self_t &&matrix) = default;

    /** Destructor */
    ~Matrix() = default;
    ///@}

    /**
     * @name Function for creating Matrix objects
     * @note These methods implement the "create idiom" and return
     * a Matrix object wrapped by a std::shared_ptr
     */
    ///@{
    /**
     * Create a distributed matrix with the dof dostribution for its rows and column
     * specified by the SparsityPattern @p sparsity_pattern.
     */
    static std::shared_ptr<self_t> create(const SparsityPattern &sparsity_pattern);
    ///@}

    /** @name Assignment operators */
    ///@{

    /**
     * Copy assignment operator. Not allowed to be used.
     */
    Matrix &operator=(const self_t &matrix) = delete;

    /**
     * Move assignment operator. Not allowed to be used.
     */
    Matrix &operator=(self_t &&matrix) = delete;
    ///@}

    /** @name Methods for getting and/or modifying the matrix entries */
    ///@{

    /**
     * Add the value @p input to the matrix entry (i,j).
     * @note @p i and @p j are the global indices of the entry (i,j).
     */
    void add_entry(const Index i,const  Index j,const Real input);

    /**
     * \brief This function add the local matrix values to the global matrix, the local-to-global ids
     * are passed as input argument.
     * \param[in] row_glob_ids The vector containing the row global ids
     * associated to the local matrix entries
     * \param[in] col_glob_ids The vector containing the colume global ids
     * associated to the local matrix entries
     * \param[in] local_matrix The local matrix that must be added to the global matrix.
     * \note The number of rows and columns of the local matrix must be equal to the dimension of the
     * local_to_global vector, otherwise an exception will be raised.
     */
    void add_block(
        const vector<Index> &row_glob_ids,
        const vector<Index> &cln_glob_ids,
        const DenseMatrix &local_matrix) ;

    /**
     * Given two global indices @p i and @p j, return the (i,j) entry of the matrix.
     */
    Real operator()(const Index i, const Index j) const;


    /** Set the row specified by the global id @p row to zero. */
    void clear_row(const Index row);
    ///@}

    /** @name Methods for retrieving the PETSc objects wrapped by this class */
    ///@{
    /**
     * todo: add documentation
     */
    Mat get_petsc_matrix() const ;

    /**
     * todo: add documentation
     */
    Mat get_petsc_matrix() ;
    ///@}

    /** @name Methods for retrieving or printing the matrix informations */
    ///@{
    void print(LogStream &out) const ;


    /** Return the number of global rows of the matrix */
    Index get_num_rows() const ;

    /** Return the number of global columns of the matrix */
    Index get_num_columns() const ;
    ///@}

    /** @name Matrix-by-vector multiplication */
    ///@{
    /**
     * Compute a sparse matrix-MultiVector multiply.
     *
     * This method computes y := beta * y + alpha * A * x.
     * If beta == 0, this operation will enjoy overwrite semantics: y's entries will be ignored,
     * and y will be overwritten with the result of the multiplication, even if it contains NaN
     * (not-a-number) floating-point entries.
     */
    void multiply_by_right_vector(const vector_t &x,vector_t &y,const Real alpha,const Real beta) const;

    /**
     * Compute a sparse matrix-MultiVector multiply.
     *
     * This method computes y^T := beta * y^T + alpha * x^T * A.
     * If beta == 0, this operation will enjoy overwrite semantics: y's entries will be ignored,
     * and y will be overwritten with the result of the multiplication, even if it contains NaN
     * (not-a-number) floating-point entries.
     */
    void multiply_by_left_vector(const vector_t &x,vector_t &y,const Real alpha,const Real beta) const;
    ///@}


    /** @name Methods for changing the filling state. */
    ///@{
    /**
     * Signal that data entry is complete, specifying domain and range maps.
     *
     * Off-node indices are distributed, indices are sorted,
     * redundant indices are eliminated, and global indices are transformed to local indices.
     */
    void fill_complete();

    /**
     * Resume operations that may change the values or structure of the matrix.
     *
     * This method must be called as a collective operation.
     * Calling fill_complete() "freezes" both the values and the structure of the matrix.
     * If you want to modify the matrix again, you must first call resume_fill().
     * You then may not call resume_fill() again on that matrix until you first call fill_complete().
     * You may make sequences of fill_complete(), resume_fill() calls as many times as you wish.
     */
    void resume_fill();
    ///@}

private:
    MPI_Comm comm_;
    Mat matrix_;

    void init(const SparsityPattern &sparsity_pattern);
};
#endif // #ifdef USE_PETSC

IGA_NAMESPACE_CLOSE

#endif /* DISTRIBUTED_MATRIX_H_ */
