//-+--------------------------------------------------------------------
// Igatools a general purpose Isogeometric analysis library.
// Copyright (C) 2012-2016  by the igatools authors (see authors.txt).
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

#ifndef __EPETRA_SOLVER_H_
#define __EPETRA_SOLVER_H_

#include <igatools/base/config.h>

#ifdef USE_TRILINOS
#include <BelosSolverFactory.hpp>
#include <BelosEpetraAdapter.hpp>
#endif // USE_TRILINOS

#include <igatools/linear_algebra/epetra_vector.h>
#include <igatools/linear_algebra/epetra_matrix.h>

IGA_NAMESPACE_OPEN

#ifdef USE_TRILINOS

namespace EpetraTools
{
using OP = Epetra_Operator;
using MV = Epetra_MultiVector;
using SolverPtr = Teuchos::RCP<Belos::SolverManager<double, MV, OP> >;
SolverPtr
create_solver(const Matrix &A, Vector &x, const Vector &b,
              const std::string &solver_type = "CG",
              const Real tolerance = 1.0e-8,
              const int max_num_iters = 400);
}

#endif // USE_TRILINOS

IGA_NAMESPACE_CLOSE

#endif
