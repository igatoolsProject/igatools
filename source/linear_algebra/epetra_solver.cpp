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

#include <igatools/linear_algebra/epetra_solver.h>
#include "ml_epetra_preconditioner.h"

IGA_NAMESPACE_OPEN

namespace EpetraTools
{

    SolverPtr create_solver(MatrixPtr A, VectorPtr x, VectorPtr b)
    {
    	using Teuchos::ParameterList;
    	using Teuchos::parameterList;
    	using Teuchos::RCP;
    	using Teuchos::rcp;
    	Belos::SolverFactory<double, MV, OP> factory;
    	RCP<ParameterList> solverParams = parameterList();
    	solverParams->set ("Num Blocks", 40);
    	solverParams->set ("Maximum Iterations", 400);
    	solverParams->set ("Convergence Tolerance", 1.0e-8);

    	SolverPtr solver =
    			factory.create ("CG", solverParams);
    	RCP<Belos::LinearProblem<double, MV, OP> > problem =
    			rcp (new Belos::LinearProblem<double, MV, OP> (
    					rcp<OP>(A.get(),false),
						rcp<MV>(x.get(),false),
						rcp<MV>(b.get(),false)));

    	RCP<ML_Epetra::MultiLevelPreconditioner> Prec =
    			rcp( new ML_Epetra::MultiLevelPreconditioner(*(A.get()), true));

    	RCP<Belos::EpetraPrecOp> belosPrec = rcp(new Belos::EpetraPrecOp(Prec));
    	problem->setLeftPrec(belosPrec);
    	problem->setProblem();

    	solver->setProblem (problem);

    	return solver;
    }

};

IGA_NAMESPACE_CLOSE

