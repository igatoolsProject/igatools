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

#ifndef TRILINOS_TOOLS_H_
#define TRILINOS_TOOLS_H_

#include <igatools/base/config.h>

#include <igatools/basis_functions/space_manager.h>

#ifdef USE_TRILINOS
#include <Tpetra_CrsMatrix.hpp>

IGA_NAMESPACE_OPEN





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

/** Type alias for the communicator. */
using Comm = Teuchos::Comm<int>;
using CommPtr = Teuchos::RCP<const Comm>;

/** Type alias for the dofs map across the processors */
using DofsMap = Tpetra::Map<LO,GO>;
using DofsMapPtr = Teuchos::RCP<const DofsMap>;

/** Type alias for the connectivty graph */
using Graph = Tpetra::CrsGraph<LO,GO>;
using GraphPtr = Teuchos::RCP<Graph>;

/** Type alias for the Trilinos matrix. */
using MatrixImpl = Tpetra::CrsMatrix<Real,LO,GO>;
using MatrixImplPtr = Teuchos::RCP<MatrixImpl>;


namespace trilinos_tools
{
DofsMapPtr build_row_map(const SpaceManager &space_manager, const CommPtr comm);

DofsMapPtr build_col_map(const SpaceManager &space_manager, const CommPtr comm);

GraphPtr build_graph(const SpaceManager &space_manager,const DofsMapPtr row_map,const DofsMapPtr col_map);

MatrixImplPtr build_matrix(GraphPtr graph);
};


IGA_NAMESPACE_CLOSE

#endif // #ifdef USE_TRILINOS

#endif // #ifndef TRILINOS_TOOLS_H_
