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

#if 0
#include <igatools/linear_algebra/linear_solver.h>

using std::shared_ptr;


#ifdef USE_TRILINOS
#include <BelosTpetraAdapter.hpp>
#include <BelosEpetraAdapter.hpp>
#include <BelosSolverFactory.hpp>


#include "Amesos2.hpp"
#include "Amesos2_Version.hpp"
#include "Amesos.h"



#ifdef REAL_IS_LONG_DOUBLE
namespace Teuchos
{
template<>
Real
ScalarTraits<Real>::zero()
{
    return Real(0.0);
}

template<>
Real
ScalarTraits<Real>::one()
{
    return Real(1.0);
}

template<>
Real
ScalarTraits<Real>::eps()
{
    return std::numeric_limits<Real>::epsilon();
}

template<>
string
ScalarTraits<Real>::name()
{
    return "Real";
}

template<>
Real
ScalarTraits<Real>::rmax()
{
    return std::numeric_limits<Real>::max();
}

template<>
auto
ScalarTraits<Real>::Real(Real number) -> magnitudeType
{
    return std::Real(number);
}

template<>
auto
ScalarTraits<Real>::conjugate(Real number) -> magnitudeType
{
    return std::Real(number);
}

template<>
auto
ScalarTraits<Real>::magnitude(Real number) -> magnitudeType
{
    return std::abs(number);
}

template<>
auto
ScalarTraits<Real>::squareroot(Real number) -> magnitudeType
{
    return std::sqrt(number);
}


template<>
Real
EnhancedNumberTraits<Real>::defaultStep()
{
    return Real(1.0);
}

template<>
unsigned short
EnhancedNumberTraits<Real>::defaultPrecision()
{
    return 100;
}

};
#endif // #ifdef REAL_IS_LONG_DOUBLE

using namespace Teuchos ;

#endif // #ifdef USE_TRILINOS




IGA_NAMESPACE_OPEN


#ifdef USE_TRILINOS






template<TrilinosImpl trilinos_impl>
LinearSolverIterativeTrilinos<trilinos_impl>::
LinearSolverIterativeTrilinos(
    const SolverType solver_type,
    const PreconditionerType prec_type,
    const Real tolerance,
    const int max_num_iter)
    :
    solver_params_(parameterList())
{
    // map the SolverType enum elements to the name aliases used by Belos
    solver_type_enum_to_alias_[to_integral(SolverType::GMRES)] = "GMRES";
    solver_type_enum_to_alias_[to_integral(SolverType::FlexibleGMRES)] = "Flexible GMRES";
    solver_type_enum_to_alias_[to_integral(SolverType::CG)] = "CG";
    solver_type_enum_to_alias_[to_integral(SolverType::StochasticCG)] = "Stochastic CG";
    solver_type_enum_to_alias_[to_integral(SolverType::RecyclingCG)] = "Recycling CG";
    solver_type_enum_to_alias_[to_integral(SolverType::RecyclingGMRES)] = "Recycling GMRES";
    solver_type_enum_to_alias_[to_integral(SolverType::PseudoBlockGMRES)] = "Pseudo Block GMRES";
    solver_type_enum_to_alias_[to_integral(SolverType::PseudoBlockCG)] = "Pseudo Block CG";
    solver_type_enum_to_alias_[to_integral(SolverType::TFQMR)] = "TFQMR";
    solver_type_enum_to_alias_[to_integral(SolverType::GCRODR)] = "GCRODR";
    solver_type_enum_to_alias_[to_integral(SolverType::MINRES)] = "MINRES";
    solver_type_enum_to_alias_[to_integral(SolverType::PCPG)] = "PCPG";


    solver_name_ = solver_type_enum_to_alias_[to_integral(solver_type)] ;


    preconditioner_type_ = prec_type;



    Belos::SolverFactory<Real,Vec,Op> factory;
    //
    // "Num Blocks" = Maximum number of Krylov vectors to store.  This
    // is also the restart length.  "Block" here refers to the ability
    // of this particular solver (and many other Belos solvers) to solve
    // multiple linear systems at a time, even though we are only solving
    // one linear system.
    solver_params_->set("Num Blocks", 1000);
    solver_params_->set("Maximum Iterations", max_num_iter);
    solver_params_->set("Convergence Tolerance", tolerance);


//    int verb = Belos::Warnings + Belos::Errors + Belos::FinalSummary + Belos::TimingDetails;
    int verb = Belos::Warnings + Belos::Errors;
    solver_params_->set("Verbosity", verb);

    // Create the solver.
    solver_ = factory.create(solver_name_, solver_params_);
}


template<TrilinosImpl trilinos_impl>
shared_ptr<LinearSolverIterativeTrilinos<trilinos_impl> >
LinearSolverIterativeTrilinos<trilinos_impl>::
create(const SolverType solver_type,
       const PreconditionerType prec_type,
       const Real tolerance,
       const int max_num_iter)
{
    return shared_ptr<LinearSolverIterativeTrilinos<trilinos_impl>>(
               new LinearSolverIterativeTrilinos<trilinos_impl>(
                   solver_type, prec_type, tolerance, max_num_iter));
}



template<TrilinosImpl trilinos_impl>
void
LinearSolverIterativeTrilinos<trilinos_impl>::
set_solver_parameters(Teuchos::RCP<Teuchos::ParameterList> solver_params)
{
    solver_params_ = solver_params;
    solver_->setParameters(solver_params_);
}

template<TrilinosImpl trilinos_impl>
void
LinearSolverIterativeTrilinos<trilinos_impl>::
set_max_num_iterations(const int max_num_iter)
{
    solver_params_->set("Maximum Iterations", max_num_iter);
    solver_->setParameters(solver_params_);
}

template<TrilinosImpl trilinos_impl>
void
LinearSolverIterativeTrilinos<trilinos_impl>::
set_tolerance(const Real tolerance)
{
    solver_params_->set("Convergence Tolerance", tolerance);
    solver_->setParameters(solver_params_);
}


template<TrilinosImpl trilinos_impl>
void
LinearSolverIterativeTrilinos<trilinos_impl>::
solve(const Matrix<la_pack> &A, const Vector<la_pack> &b, Vector<la_pack> &x)
{
    Matrix<la_pack> &A_nonconst = const_cast<Matrix<la_pack> &>(A);

    if (preconditioner_type_ == PreconditionerType::ILU0)
        preconditioner_ = TrilinosTools<trilinos_impl>::compute_preconditioner_ILU(
                              A_nonconst.get_trilinos_matrix(),0);
    else if (preconditioner_type_ == PreconditionerType::ILU1)
        preconditioner_ = TrilinosTools<trilinos_impl>::compute_preconditioner_ILU(
                              A_nonconst.get_trilinos_matrix(),1);
    else if (preconditioner_type_ == PreconditionerType::ILU2)
        preconditioner_ = TrilinosTools<trilinos_impl>::compute_preconditioner_ILU(
                              A_nonconst.get_trilinos_matrix(),2);



    // Create a LinearProblem struct with the problem to solve.
    // A, X, B, and M are passed by (smart) pointer, not copied.
    auto linear_system = rcp(
                             new Belos::LinearProblem<Real,Vec,Op>(
                                 A_nonconst.get_trilinos_matrix(),
                                 x.get_trilinos_vector(),
                                 const_cast<Vector<la_pack> &>(b).get_trilinos_vector()));

    if (!preconditioner_ .is_null())
        linear_system->setLeftPrec(preconditioner_);

    linear_system->setProblem() ;

    solver_->setProblem(linear_system);

    // Attempt to solve the linear system.  result == Belos::Converged
    // means that it was solved to the desired tolerance.  This call
    // overwrites X with the computed approximate solution.
    Belos::ReturnType result = solver_->solve();

    AssertThrow(result == Belos::ReturnType::Converged, ExcMessage("No convergence."));
}


template<TrilinosImpl trilinos_impl>
Real
LinearSolverIterativeTrilinos<trilinos_impl>::
get_achieved_tolerance() const
{
    return solver_->achievedTol() ;
}

template<TrilinosImpl trilinos_impl>
int
LinearSolverIterativeTrilinos<trilinos_impl>::
get_num_iterations() const
{
    return solver_->getNumIters() ;
}



template<TrilinosImpl trilinos_impl>
string
LinearSolverIterativeTrilinos<trilinos_impl>::
get_solver_name() const
{
    return solver_name_;
}



template<TrilinosImpl trilinos_impl>
auto
LinearSolverIterativeTrilinos<trilinos_impl>::
get_preconditioner_type() const -> PreconditionerType
{
    return preconditioner_type_;
}



template<TrilinosImpl trilinos_impl>
Teuchos::RCP<Teuchos::ParameterList>
LinearSolverIterativeTrilinos<trilinos_impl>::
get_solver_parameters() const
{
    return solver_params_;
}






template<>
LinearSolverDirectTrilinos<TrilinosImpl::tpetra>::
LinearSolverDirectTrilinos(const SolverType &solver_type)
{
    if (solver_type == SolverType::SUPERLU)
        solver_name_ = "superlu";
    else if (solver_type == SolverType::LAPACK)
        solver_name_ = "lapack";
    else if (solver_type == SolverType::PARDISO_MKL)
        solver_name_ = "pardiso_mkl";
    else if (solver_type == SolverType::UMFPACK)
        solver_name_ = "umfpack";
};

template<>
LinearSolverDirectTrilinos<TrilinosImpl::epetra>::
LinearSolverDirectTrilinos(const SolverType &solver_type)
{
    if (solver_type == SolverType::SUPERLU)
        solver_name_ = "Superlu";
    else if (solver_type == SolverType::LAPACK)
        solver_name_ = "Lapack";
    else if (solver_type == SolverType::PARDISO_MKL)
        solver_name_ = "Pardiso";
    else if (solver_type == SolverType::UMFPACK)
        solver_name_ = "Umfpack";
    else if (solver_type == SolverType::MUMPS)
        solver_name_ = "Mumps";
    else if (solver_type == SolverType::KLU)
        solver_name_ = "Klu";
};



template<TrilinosImpl trilinos_impl>
shared_ptr<LinearSolverDirectTrilinos<trilinos_impl> >
LinearSolverDirectTrilinos<trilinos_impl>::
create(const SolverType &solver_type)
{
    return shared_ptr<LinearSolverDirectTrilinos<trilinos_impl>>(
               new LinearSolverDirectTrilinos<trilinos_impl>(solver_type));
};


template<>
void
LinearSolverDirectTrilinos<TrilinosImpl::epetra>::
solve(const Matrix<la_pack> &matrix,
      const Vector<la_pack> &rhs,
      Vector<la_pack> &solution)
{
    Assert(matrix.get_trilinos_matrix().get() != nullptr, ExcNullPtr());
    Assert(rhs.get_trilinos_vector().get() != nullptr, ExcNullPtr());
    Assert(solution.get_trilinos_vector().get() != nullptr, ExcNullPtr());

    using Mat = Matrix<la_pack>;
    using Vec = Vector<la_pack>;
    auto m = const_cast<Mat &>(matrix).get_trilinos_matrix();
    auto r = const_cast<Vec &>(rhs).get_trilinos_vector();
    auto s = solution.get_trilinos_vector();
    auto epetra_linear_problem =
        Epetra_LinearProblem(m.get(),s.get(),r.get());

    Amesos amesos_factory;
    auto solver =  amesos_factory.Create(solver_name_,epetra_linear_problem);
    Assert(solver != nullptr,ExcNullPtr());

    solver->SymbolicFactorization();
    solver->NumericFactorization();
    solver->Solve();
}



template<>
void
LinearSolverDirectTrilinos<TrilinosImpl::tpetra>::
solve(const Matrix<la_pack> &matrix,
      const Vector<la_pack> &rhs,
      Vector<la_pack> &solution)
{
    Assert(matrix.get_trilinos_matrix().get() != nullptr, ExcNullPtr());
    Assert(rhs.get_trilinos_vector().get() != nullptr, ExcNullPtr());
    Assert(solution.get_trilinos_vector().get() != nullptr, ExcNullPtr());

    using Mat = Matrix<la_pack>;
    using Vec = Vector<la_pack>;
    auto solver = Amesos2::create<
                  typename Mat::WrappedMatrix,
                  typename Vec::WrappedVector>(
                      solver_name_,
                      const_cast<Mat &>(matrix).get_trilinos_matrix(),
                      solution.get_trilinos_vector(),
                      const_cast<Vec &>(rhs).get_trilinos_vector());

    solver->symbolicFactorization().numericFactorization().solve();
}



template<TrilinosImpl trilinos_impl>
string
LinearSolverDirectTrilinos<trilinos_impl>::
get_solver_name() const
{
    return solver_name_;
}


template class LinearSolverIterativeTrilinos<TrilinosImpl::tpetra>;
template class LinearSolverIterativeTrilinos<TrilinosImpl::epetra>;

template class LinearSolverDirectTrilinos<TrilinosImpl::tpetra>;
template class LinearSolverDirectTrilinos<TrilinosImpl::epetra>;



#endif // #ifdef USE_TRILINOS







#ifdef USE_PETSC

LinearSolverIterative<LAPack::petsc>::
LinearSolverIterative(const SolverType solver_type, const Real tolerance, const int max_num_iter)
{
    PetscErrorCode ierr;
    comm_ = PETSC_COMM_WORLD;
    std::string prec_name;

    // map the SolverType enum elements to the name aliases used by PETSc
    solver_type_enum_to_alias_[to_integral(SolverType::GMRES)] = "gmres";
    solver_type_enum_to_alias_[to_integral(SolverType::CG)] = "cg";
    solver_type_enum_to_alias_[to_integral(SolverType::LU)] = "preonly";

    solver_name_ = solver_type_enum_to_alias_[to_integral(solver_type)] ;


    if (solver_type == SolverType::LU)
    {
        preconditioner_type_ = PreconditionerType::ILU;
        prec_name = "lu" ;
    }
    else
    {
        preconditioner_type_ = PreconditionerType::NONE;
        prec_name = "none" ;
    }

    ierr = KSPCreate(comm_, &ksp_);
    ierr = KSPGetPC(ksp_,&pc_);

    ierr = PCSetType(pc_,prec_name.c_str());
    ierr = KSPSetType(ksp_,solver_name_.c_str());

    ierr = KSPSetTolerances(ksp_,tolerance,PETSC_DEFAULT,PETSC_DEFAULT,max_num_iter);
}

LinearSolverIterative<LAPack::petsc>::
LinearSolverIterative(const SolverType solver_type, const PreconditionerType prec_type,
                      const Real tolerance, const int max_num_iter)
{
    PetscErrorCode ierr;
    comm_ = PETSC_COMM_WORLD;

    // map the SolverType enum elements to the name aliases used by PETSc
    solver_type_enum_to_alias_[to_integral(SolverType::GMRES)] = "gmres";
    solver_type_enum_to_alias_[to_integral(SolverType::CG)] = "cg";
    solver_type_enum_to_alias_[to_integral(SolverType::LU)] = "preonly";

    // map the PreconditionerType enum elements to the name aliases used by PETSc
    prec_type_enum_to_alias_[to_integral(PreconditionerType::NONE)] = "none";
    prec_type_enum_to_alias_[to_integral(PreconditionerType::ILU)] = "ilu";
    prec_type_enum_to_alias_[to_integral(PreconditionerType::JACOBI)] = "jacobi";


    preconditioner_type_ = prec_type;

    solver_name_ = solver_type_enum_to_alias_[to_integral(solver_type)] ;
    std::string prec_name = prec_type_enum_to_alias_[to_integral(prec_type)] ;

    if (solver_type == SolverType::LU)
        prec_name = "lu" ;

    ierr = KSPCreate(comm_, &ksp_);
    ierr = KSPGetPC(ksp_,&pc_);

    ierr = PCSetType(pc_,prec_name.c_str());
    ierr = KSPSetType(ksp_, solver_name_.c_str());

    ierr = KSPSetTolerances(ksp_,tolerance,PETSC_DEFAULT,PETSC_DEFAULT,max_num_iter);
}

/*
void
LinearSolverIterative<LAPack::petsc>::
set_solver_parameters(Teuchos::RCP<Teuchos::ParameterList> solver_params)
{
    solver_params_ = solver_params;
    solver_->setParameters(solver_params_);
}
//*/

void
LinearSolverIterative<LAPack::petsc>::
set_max_num_iterations(const int max_num_iter)
{
    PetscErrorCode ierr;
    PetscReal rtol, abstol, dtol;
    PetscInt maxits;
    ierr = KSPGetTolerances(ksp_, &rtol, &abstol, &dtol, &maxits);
    ierr = KSPSetTolerances(ksp_, rtol, abstol, dtol, max_num_iter);

}

/**
 * Set the level that residual norms must reach to decide convergence.
 */
void
LinearSolverIterative<LAPack::petsc>::
set_tolerance(const Real tolerance)
{
    PetscErrorCode ierr;
    PetscReal rtol, abstol, dtol;
    PetscInt maxits;
    ierr = KSPGetTolerances(ksp_, &rtol, &abstol, &dtol, &maxits);
    ierr = KSPSetTolerances(ksp_, tolerance, abstol, dtol, maxits);
}


void
LinearSolverIterative<LAPack::petsc>::
solve(Matrix<LAPack::petsc> &A,
      Vector<LAPack::petsc> &b,
      Vector<LAPack::petsc> &x)
{
    PetscErrorCode ierr;
    ierr = KSPSetOperators(ksp_,A.get_petsc_matrix(),A.get_petsc_matrix(),SAME_NONZERO_PATTERN);
    ierr = KSPSolve(ksp_,b.get_petsc_vector(),x.get_petsc_vector());
}


Real
LinearSolverIterative<LAPack::petsc>::
get_achieved_tolerance() const
{
    PetscErrorCode ierr;
    PetscReal achieved_tol;
    ierr = KSPGetResidualNorm(ksp_, &achieved_tol);

    return achieved_tol;
}

int
LinearSolverIterative<LAPack::petsc>::
get_num_iterations() const
{
    PetscErrorCode ierr;
    int num_iters;
    ierr = KSPGetIterationNumber(ksp_, &num_iters);

    return num_iters;
}


string
LinearSolverIterative<LAPack::petsc>::
get_solver_name() const
{
    return solver_name_;
}



auto
LinearSolverIterative<LAPack::petsc>::
get_preconditioner_type() const -> PreconditionerType
{
    return preconditioner_type_;
}

#endif // #ifdef USE_PETSC



IGA_NAMESPACE_CLOSE
#endif

