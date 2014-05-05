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


#include <igatools/linear_algebra/linear_solver.h>


#ifdef USE_TRILINOS
#include <BelosTpetraAdapter.hpp>
#include <BelosSolverFactory.hpp>


#ifdef REAL_IS_LONG_DOUBLE
namespace Teuchos
{
template<>
iga::Real
ScalarTraits<iga::Real>::zero()
{
    return iga::Real(0.0);
}

template<>
iga::Real
ScalarTraits<iga::Real>::one()
{
    return iga::Real(1.0);
}

template<>
iga::Real
ScalarTraits<iga::Real>::eps()
{
    return std::numeric_limits<iga::Real>::epsilon();
}

template<>
string
ScalarTraits<iga::Real>::name()
{
    return "iga::Real";
}

template<>
iga::Real
ScalarTraits<iga::Real>::rmax()
{
    return std::numeric_limits<iga::Real>::max();
}

template<>
auto
ScalarTraits<iga::Real>::Real(iga::Real number) -> magnitudeType
{
    return std::Real(number);
}

template<>
auto
ScalarTraits<iga::Real>::conjugate(iga::Real number) -> magnitudeType
{
    return std::Real(number);
}

template<>
auto
ScalarTraits<iga::Real>::magnitude(iga::Real number) -> magnitudeType
{
    return std::abs(number);
}

template<>
auto
ScalarTraits<iga::Real>::squareroot(iga::Real number) -> magnitudeType
{
    return std::sqrt(number);
}


template<>
iga::Real
EnhancedNumberTraits<iga::Real>::defaultStep()
{
    return iga::Real(1.0);
}

template<>
unsigned short
EnhancedNumberTraits<iga::Real>::defaultPrecision()
{
    return 100;
}

};
#endif // #ifdef REAL_IS_LONG_DOUBLE

using namespace Teuchos ;

#endif // #ifdef USE_TRILINOS




IGA_NAMESPACE_OPEN


#ifdef USE_TRILINOS

LinearSolver<LinearAlgebraPackage::trilinos>::
LinearSolver(const Type solver_type, const Real tolerance, const int max_num_iter)
    :
    solver_params_(parameterList())
{
    // map the SolverType enum elements to the name aliases used by Belos
    solver_type_enum_to_alias_[to_integral(Type::GMRES)] = "GMRES";
    solver_type_enum_to_alias_[to_integral(Type::FlexibleGMRES)] = "Flexible GMRES";
    solver_type_enum_to_alias_[to_integral(Type::CG)] = "CG";
    solver_type_enum_to_alias_[to_integral(Type::StochasticCG)] = "Stochastic CG";
    solver_type_enum_to_alias_[to_integral(Type::RecyclingCG)] = "Recycling CG";
    solver_type_enum_to_alias_[to_integral(Type::RecyclingGMRES)] = "Recycling GMRES";
    solver_type_enum_to_alias_[to_integral(Type::PseudoBlockGMRES)] = "Pseudo Block GMRES";
    solver_type_enum_to_alias_[to_integral(Type::PseudoBlockCG)] = "Pseudo Block CG";


    const std::string solver_name = solver_type_enum_to_alias_[to_integral(solver_type)] ;

    Belos::SolverFactory<Real,vector_t,matrix_t> factory;
    //
    // "Num Blocks" = Maximum number of Krylov vectors to store.  This
    // is also the restart length.  "Block" here refers to the ability
    // of this particular solver (and many other Belos solvers) to solve
    // multiple linear systems at a time, even though we are only solving
    // one linear system in this example.
    solver_params_->set("Num Blocks", 40);
    solver_params_->set("Maximum Iterations", max_num_iter);
    solver_params_->set("Convergence Tolerance", tolerance);

    // Create the solver.
    solver_ = factory.create(solver_name,solver_params_);
}


void
LinearSolver<LinearAlgebraPackage::trilinos>::
set_solver_parameters(Teuchos::RCP<Teuchos::ParameterList> solver_params)
{
    solver_params_ = solver_params;
    solver_->setParameters(solver_params_);
}

void
LinearSolver<LinearAlgebraPackage::trilinos>::
set_max_num_iterations(const int max_num_iter)
{
    solver_params_->set("Maximum Iterations", max_num_iter);
    solver_->setParameters(solver_params_);
}

/**
 * Set the level that residual norms must reach to decide convergence.
 */
void
LinearSolver<LinearAlgebraPackage::trilinos>::
set_tolerance(const Real tolerance)
{
    solver_params_->set("Convergence Tolerance", tolerance);
    solver_->setParameters(solver_params_);
}


void
LinearSolver<LinearAlgebraPackage::trilinos>::
solve(Matrix<LinearAlgebraPackage::trilinos> &A,
      Vector<LinearAlgebraPackage::trilinos> &b,
      Vector<LinearAlgebraPackage::trilinos> &x)
{
    // Create a LinearProblem struct with the problem to solve.
    // A, X, B, and M are passed by (smart) pointer, not copied.
    auto linear_system = rcp(
                             new Belos::LinearProblem<Real,vector_t,matrix_t>(
                                 A.get_trilinos_matrix(),
                                 x.get_trilinos_vector(),
                                 b.get_trilinos_vector()));
    linear_system->setProblem() ;

    solver_->setProblem(linear_system);

    // Attempt to solve the linear system.  result == Belos::Converged
    // means that it was solved to the desired tolerance.  This call
    // overwrites X with the computed approximate solution.
    Belos::ReturnType result = solver_->solve();

    AssertThrow(result == Belos::ReturnType::Converged, ExcMessage("No convergence."));
}


Real
LinearSolver<LinearAlgebraPackage::trilinos>::
get_achieved_tolerance() const
{
    return solver_->achievedTol() ;
}

int
LinearSolver<LinearAlgebraPackage::trilinos>::
get_num_iterations() const
{
    return solver_->getNumIters() ;
}

#endif // #ifdef USE_TRILINOS







#ifdef USE_PETSC

LinearSolver<LinearAlgebraPackage::petsc>::
LinearSolver(const Type solver_type, const Real tolerance, const int max_num_iter)
{
    Assert(false,ExcNotImplemented());
    AssertThrow(false,ExcNotImplemented());
    /*
    // map the SolverType enum elements to the name aliases used by PETSC::????????
    solver_type_enum_to_alias_[to_integral(Type::GMRES)] = "GMRES";
    solver_type_enum_to_alias_[to_integral(Type::CG)] = "CG";


    const std::string solver_name = solver_type_enum_to_alias_[to_integral(solver_type)] ;

    Belos::SolverFactory<Real,vector_t,matrix_t> factory;
    //
    // "Num Blocks" = Maximum number of Krylov vectors to store.  This
    // is also the restart length.  "Block" here refers to the ability
    // of this particular solver (and many other Belos solvers) to solve
    // multiple linear systems at a time, even though we are only solving
    // one linear system in this example.
    solver_params_->set("Num Blocks", 40);
    solver_params_->set("Maximum Iterations", max_num_iter);
    solver_params_->set("Convergence Tolerance", tolerance);

    // Create the solver.
    solver_ = factory.create(solver_name,solver_params_);
       //*/
}

/*
void
LinearSolver<LinearAlgebraPackage::petsc>::
set_solver_parameters(Teuchos::RCP<Teuchos::ParameterList> solver_params)
{
    solver_params_ = solver_params;
    solver_->setParameters(solver_params_);
}
//*/

void
LinearSolver<LinearAlgebraPackage::petsc>::
set_max_num_iterations(const int max_num_iter)
{
    Assert(false,ExcNotImplemented());
    AssertThrow(false,ExcNotImplemented());
    /*
    solver_params_->set("Maximum Iterations", max_num_iter);
    solver_->setParameters(solver_params_);
    //*/
}

/**
 * Set the level that residual norms must reach to decide convergence.
 */
void
LinearSolver<LinearAlgebraPackage::petsc>::
set_tolerance(const Real tolerance)
{
    Assert(false,ExcNotImplemented());
    AssertThrow(false,ExcNotImplemented());
    /*
    solver_params_->set("Convergence Tolerance", tolerance);
    solver_->setParameters(solver_params_);
    //*/
}


void
LinearSolver<LinearAlgebraPackage::petsc>::
solve(Matrix<LinearAlgebraPackage::petsc> &A,
      Vector<LinearAlgebraPackage::petsc> &b,
      Vector<LinearAlgebraPackage::petsc> &x)
{
    Assert(false,ExcNotImplemented());
    AssertThrow(false,ExcNotImplemented());
    /*
    // Create a LinearProblem struct with the problem to solve.
    // A, X, B, and M are passed by (smart) pointer, not copied.
    auto linear_system = rcp(
                             new Belos::LinearProblem<Real,vector_t,matrix_t>(
                                 A.get_petsc_matrix(),
                                 x.get_petsc_vector(),
                                 b.get_petsc_vector()));
    linear_system->setProblem() ;

    solver_->setProblem(linear_system);

    // Attempt to solve the linear system.  result == Belos::Converged
    // means that it was solved to the desired tolerance.  This call
    // overwrites X with the computed approximate solution.
    Belos::ReturnType result = solver_->solve();

    AssertThrow(result == Belos::ReturnType::Converged, ExcMessage("No convergence."));
    //*/
}


Real
LinearSolver<LinearAlgebraPackage::petsc>::
get_achieved_tolerance() const
{
    Assert(false,ExcNotImplemented());
    AssertThrow(false,ExcNotImplemented());
    /*
    return solver_->achievedTol() ;
    //*/

    return 0.0;
}

int
LinearSolver<LinearAlgebraPackage::petsc>::
get_num_iterations() const
{
    Assert(false,ExcNotImplemented());
    AssertThrow(false,ExcNotImplemented());
    /*
    return solver_->getNumIters() ;
    //*/

    return 0;
}

#endif // #ifdef USE_PETSC


IGA_NAMESPACE_CLOSE


