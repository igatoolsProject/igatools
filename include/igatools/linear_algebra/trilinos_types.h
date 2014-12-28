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

#ifndef TRILINOS_TYPES_H_
#define TRILINOS_TYPES_H_

#include <igatools/base/config.h>


#ifdef USE_TRILINOS

#include <Tpetra_CrsMatrix.hpp>
#include <Epetra_CrsMatrix.h>
#include <Epetra_Vector.h>
#include <Epetra_SerialComm.h>

IGA_NAMESPACE_OPEN






/**
 * @brief Type for specifying which kind of Trilinos implementation to use (Tpetra-based or Epetra-based).
 */
enum class TrilinosImpl : int
{
    /** Use the Trilinos Tpetra linear algebra implementation.*/
    tpetra = 1,

    /** Use the Trilinos Epetra linear algebra implementation.*/
    epetra = 2
};



/**
 * @brief Struct containing the type aliases for the Trilinos objects.
 */
template <TrilinosImpl trilinos_impl>
struct TrilinosTypes;


/**
 * @brief Struct containing the type aliases for the Trilinos/Tpetra implementation.
 */
template <>
struct TrilinosTypes<TrilinosImpl::tpetra>
{
    static const LAPack la_pack = LAPack::trilinos_tpetra;

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
//    using Node = typename Kokkos::SerialNode;

    /** Type alias for the local indices. */
    using LO = Index;

    /** Type alias for the global indices. */
    using GO = Index;

    /** Type alias for the communicator. */
    using Comm = Teuchos::Comm<int>;
    using CommPtr = Teuchos::RCP<const Comm>;

    /** Type alias for the dofs map across the processors. */
    using Map = Tpetra::Map<LO,GO>;
    using MapPtr = Teuchos::RCP<const Map>;

    /** Type alias for the connectivty graph. */
    using Graph = Tpetra::CrsGraph<LO,GO>;
    using GraphPtr = Teuchos::RCP<Graph>;

    /** Type alias for the Trilinos matrix.*/
    using Matrix = Tpetra::CrsMatrix<Real,LO,GO>;
    using MatrixPtr = Teuchos::RCP<Matrix>;

    /** Type alias for the Trilinos (multi) vector.*/
    using Vector = Tpetra::MultiVector<Real,LO,GO>;
    using VectorPtr = Teuchos::RCP<Vector>;


    /** Type alias for the linear operator.*/
    using Op = Tpetra::Operator<Real,LO,GO> ;
    using OpPtr = Teuchos::RCP<Op>;

    /** Available direct solvers.*/
    enum class DirectSolver
    {
        SUPERLU = 0,
        LAPACK = 1,
        PARDISO_MKL = 2,
        UMFPACK = 3,
        ENUM_SIZE = 4,
    };
};

/**
 * @brief Struct containing the type aliases for the Trilinos/Epetra implementation.
 */
template <>
struct TrilinosTypes<TrilinosImpl::epetra>
{
    static const LAPack la_pack = LAPack::trilinos_epetra;

    /** Type alias for the communicator. */
    using Comm = Epetra_Comm;
    using CommPtr = Teuchos::RCP<const Comm>;

    /** Type alias for the dofs map across the processors. */
    using Map = Epetra_Map;
    using MapPtr = Teuchos::RCP<const Map>;

    /** Type alias for the connectivty graph. */
    using Graph = Epetra_CrsGraph;
    using GraphPtr = Teuchos::RCP<Graph>;

    /** Type alias for the Trilinos matrix.*/
    using Matrix = Epetra_CrsMatrix;
    using MatrixPtr = Teuchos::RCP<Matrix>;

    /** Type alias for the Trilinos (multi) vector.*/
    using Vector = Epetra_MultiVector;
    using VectorPtr = Teuchos::RCP<Vector>;


    /** Type alias for the linear operator.*/
    using Op = Epetra_Operator;
    using OpPtr = Teuchos::RCP<Op>;

    /** Available direct solvers.*/
    enum class DirectSolver
    {
        SUPERLU = 0,
        LAPACK = 1,
        PARDISO_MKL = 2,
        MUMPS = 3,
        UMFPACK = 4,
        KLU = 5,
        ENUM_SIZE = 6,
    };
};

IGA_NAMESPACE_CLOSE


#endif // #ifdef USE_TRILINOS


#endif // #ifndef TRILINOS_TYPES_H_

