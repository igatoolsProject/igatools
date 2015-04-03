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
#include <igatools/linear_algebra/trilinos_tools.h>


#include "Ifpack.h"
#include <Ifpack2_Factory.hpp>
#if 0
#include <Teuchos_oblackholestream.hpp>
#include <Teuchos_TimeMonitor.hpp>
#endif

IGA_NAMESPACE_OPEN

#ifdef USE_TRILINOS




template <>
Teuchos::RCP<Teuchos::ParameterList>
TrilinosTools<TrilinosImpl::tpetra>::
ILUKPreconditionerFactory::
get_parameter_list_for_preconditioner_creation(const int fill_level) const
{
    // The name of the type of preconditioner to use.
    const std::string precondType("RILUK");

    // Ifpack2 expects double-precision arguments here.
    const double drop_tol = 0.0;
    const double abs_threshold = 0.0;

    auto pl = Teuchos::parameterList("Preconditioner");
    pl->set("Ifpack2::Preconditioner", precondType);

    Teuchos::ParameterList prec_params("Ifpack2");
    prec_params.set("fact: iluk level-of-fill", fill_level);
    prec_params.set("fact: drop tolerance", drop_tol);
    prec_params.set("fact: absolute threshold", abs_threshold);

    pl->set("Ifpack2", prec_params);
    return pl;
}

template <>
auto
TrilinosTools<TrilinosImpl::tpetra>::
ILUKPreconditionerFactory::
create(const Teuchos::RCP<const Matrix> &A,
       const Teuchos::RCP<Teuchos::ParameterList> &plist) const -> OpPtr
{
    using Teuchos::Time;
    using Teuchos::TimeMonitor;
    using std::cout;
    using std::cerr;
    using std::endl;

    /*
    cerr << "Creating RILU preconditioner" << endl
        << "-- Configuring" << endl;
    //*/

    // An Ifpack2::Preconditioner is-a Tpetra::Operator.  Ifpack2
    // creates a Preconditioner object, but users of iterative methods
    // want a Tpetra::Operator.  That's why create() returns an
    // Operator instead of a Preconditioner.
    using prec_type = Ifpack2::Preconditioner<Real,
          typename Types::LO,
          typename Types::GO,
          typename Kokkos::SerialNode>;

    Teuchos::RCP<prec_type> prec;
    {
        Ifpack2::Factory factory;

        // Get the preconditioner type.
        const std::string precName =
            plist->get<std::string>("Ifpack2::Preconditioner");

        // Set up the preconditioner of that type.
        prec = factory.create(precName, A);

        Teuchos::ParameterList ifpack2Params;
        if (plist->isSublist("Ifpack2"))
            ifpack2Params = plist->sublist("Ifpack2");
        else
            ifpack2Params.setName("Ifpack2");

        prec->setParameters(ifpack2Params);
    }

    // Initializing
    prec->initialize();

    // Computing
    prec->compute();

    /*
    cerr << "-- Estimating condition number" << endl;
    magnitude_type condest = STM::one();
    {
        TimeMonitor mon(*condestTimer);
        condest = prec->computeCondEst(Ifpack2::Cheap);
        cerr << "---- Estimated condition number = " << condest << endl;
    }
    //*/
    return prec;
}



template <>
Teuchos::RCP<Teuchos::ParameterList>
TrilinosTools<TrilinosImpl::epetra>::
ILUKPreconditionerFactory::
get_parameter_list_for_preconditioner_creation(const int fill_level) const
{
    const double drop_tol = 0.0;
    const double abs_threshold = 0.0;

    auto prec_params = Teuchos::parameterList("Ifpack");
    prec_params->set("fact: level-of-fill", fill_level);
    prec_params->set("fact: drop tolerance", drop_tol);
    prec_params->set("fact: absolute threshold", abs_threshold);

    // IFPACK's implementation of overlapping Schwarz domain
    // decomposition offers different modes for combining overlapping
    // values on different MPI processes:
    //
    // "Add", "Zero", "Insert", "InsertAdd", "Average", "AbsMax"
    //
    // See Epetra_CombineMode.h for the meaning of each of these modes.
    prec_params->set("schwarz: combine mode", "Add");

    return prec_params;
}

template <>
auto
TrilinosTools<TrilinosImpl::epetra>::
ILUKPreconditionerFactory::
create(const Teuchos::RCP<const Matrix> &A,
       const Teuchos::RCP<Teuchos::ParameterList> &plist) const -> OpPtr
{
    // Allocate an IFPACK factory.  No data is associated with this
    // object.  It only exposes the Create() method.
    Ifpack Factory;

    Matrix &A_nonconst = const_cast<Matrix &>(*A);
    //
    // Create the preconditioner using the Factory.  Please check the
    // IFPACK documentation for valid values of PrecType.
    //
    std::string PrecType = "ILU"; // Incomplete LU factorization
    int OverlapLevel = 1;    // Must be >= 0. Ignored if Comm.NumProc() == 1.
    auto prec = Teuchos::rcp(Factory.Create(PrecType, &A_nonconst, OverlapLevel));


    // Set the parameters.
    AssertThrow(prec->SetParameters(*plist) == 0,ExcInvalidState());

    // Initialize the preconditioner. At this point the matrix must have
    // been FillComplete()'d, but the values of its entries are ignored.
    AssertThrow(prec->Initialize() == 0,ExcInvalidState());

    // Build the preconditioner by examining the values in the matrix.
    AssertThrow(prec->Compute() == 0,ExcInvalidState());

    return prec;
}



template <>
auto
TrilinosTools<TrilinosImpl::tpetra>::
build_row_map(const SpaceManager &space_manager,const CommPtr comm) -> MapPtr
{
    using LO = typename Types::LO;
    using GO = typename Types::GO;

    const auto dofs_set = space_manager.get_row_dofs();
    const vector<GO> dofs_vec(dofs_set.begin(),dofs_set.end());

    return Tpetra::createNonContigMap<LO,GO>(dofs_vec,comm);
}

template <>
auto
TrilinosTools<TrilinosImpl::tpetra>::
build_col_map(const SpaceManager &space_manager,const CommPtr comm) -> MapPtr
{
    using LO = typename Types::LO;
    using GO = typename Types::GO;

    const auto dofs_set = space_manager.get_col_dofs();
    const vector<GO> dofs_vec(dofs_set.begin(),dofs_set.end());

    return Tpetra::createNonContigMap<LO,GO>(dofs_vec,comm);
}

template <>
auto
TrilinosTools<TrilinosImpl::tpetra>::
build_graph(const SpaceManager &space_manager,const MapPtr row_map,const MapPtr col_map) -> GraphPtr
{
    auto sparsity_pattern_ptr = space_manager.get_sparsity_pattern();
    Assert(sparsity_pattern_ptr!=nullptr,ExcNullPtr());

    const auto &sparsity_pattern = *sparsity_pattern_ptr;
    using LongUInt = long unsigned int;
    Teuchos::ArrayRCP<LongUInt> n_dofs_per_row(sparsity_pattern.get_num_rows());

    Index i = 0;
    for (const auto &map_entry : sparsity_pattern)
        n_dofs_per_row[i++] = map_entry.second.size();

    auto graph = Teuchos::rcp(new Graph(row_map,col_map,n_dofs_per_row,Tpetra::StaticProfile));
    for (const auto &row : sparsity_pattern)
    {
        const Index row_id = row.first ;
        const auto &cols_id = row.second;

        auto cols_id_vec = vector<Index>(cols_id.begin(),cols_id.end());

        auto cols_id_view = Teuchos::ArrayView<const typename Types::GO>(std::move(cols_id_vec));

        graph->insertGlobalIndices(row_id,cols_id_view);
    }
    graph->fillComplete(col_map,row_map);

    return graph;
}



template <>
auto
TrilinosTools<TrilinosImpl::tpetra>::
build_matrix(GraphPtr graph) -> MatrixPtr
{
    auto matrix = Teuchos::rcp(new Matrix(graph));
    matrix->setAllToScalar(0.0);

    return matrix;
}

template <>
auto
TrilinosTools<TrilinosImpl::tpetra>::
build_vector(const MapPtr map) -> VectorPtr
{
    using LO = typename Types::LO;
    using GO = typename Types::GO;
    return Tpetra::createMultiVector<Real,LO,GO>(map,1);
}

template <TrilinosImpl trilinos_impl>
auto
//TrilinosTools<TrilinosImpl::tpetra>::
TrilinosTools<trilinos_impl>::
compute_preconditioner_ILU(const MatrixPtr matrix, const int fill_level) -> OpPtr
{
    // Create a factory for creating a preconditioner.
    ILUKPreconditionerFactory factory;

    // ParameterList for creating an Ifpack2 preconditioner.
    auto param_list = factory.get_parameter_list_for_preconditioner_creation(fill_level);

    // Compute the preconditioner using the input matrix.
    // The matrix itself is not modified.
    auto preconditioner = factory.create(matrix, param_list);

    return preconditioner;
}





template <>
auto
TrilinosTools<TrilinosImpl::epetra>::
build_row_map(const SpaceManager &space_manager,const CommPtr comm) -> MapPtr
{
    const auto dofs_set = space_manager.get_row_dofs();
    const vector<Index> dofs_vec(dofs_set.begin(),dofs_set.end());

    return Teuchos::rcp(new Epetra_Map(-1,dofs_vec.size(),dofs_vec.data(),0,*comm));
}


template <>
auto
TrilinosTools<TrilinosImpl::epetra>::
build_col_map(const SpaceManager &space_manager,const CommPtr comm) -> MapPtr
{
    const auto dofs_set = space_manager.get_col_dofs();
    const vector<Index> dofs_vec(dofs_set.begin(),dofs_set.end());

    return Teuchos::rcp(new Epetra_Map(-1,dofs_vec.size(),dofs_vec.data(),0,*comm));
}


template <>
auto
TrilinosTools<TrilinosImpl::epetra>::
build_graph(const SpaceManager &space_manager,const MapPtr row_map,const MapPtr col_map) -> GraphPtr
{
    auto sparsity_pattern_ptr = space_manager.get_sparsity_pattern();
    Assert(sparsity_pattern_ptr!=nullptr,ExcNullPtr());

    const auto &sparsity_pattern = *sparsity_pattern_ptr;
    vector<Index> n_dofs_per_row(sparsity_pattern.get_num_rows());

    Index i = 0;
    for (const auto &map_entry : sparsity_pattern)
        n_dofs_per_row[i++] = map_entry.second.size();


    const bool is_static_profile = true;
    auto graph = Teuchos::rcp(
        new Graph(Epetra_DataAccess::Copy, *row_map, *col_map,n_dofs_per_row.data(),is_static_profile));
    for (const auto &row : sparsity_pattern)
    {
        const Index row_id = row.first ;
        const auto &cols_id = row.second;

        auto cols_id_vec = vector<Index>(cols_id.begin(),cols_id.end());

//        auto cols_id_view = Teuchos::ArrayView<const GO>(std::move(cols_id_vec));

        graph->InsertGlobalIndices(row_id,cols_id_vec.size(),cols_id_vec.data());
    }
    graph->FillComplete(*col_map,*row_map);

    return graph;
}



template <>
auto
TrilinosTools<TrilinosImpl::epetra>::
build_matrix(GraphPtr graph) -> MatrixPtr
{
    auto matrix = Teuchos::rcp(new Matrix(Epetra_DataAccess::Copy,*graph));
    matrix->PutScalar(0.0);

    return matrix;
}

template <>
auto
TrilinosTools<TrilinosImpl::epetra>::
build_vector(const MapPtr map) -> VectorPtr
{
    return Teuchos::rcp(new Epetra_MultiVector(*map,1));
}

#if 0
template <>
auto
TrilinosTools<TrilinosImpl::epetra>::
compute_preconditioner_ILU(const MatrixPtr matrix, const int fill_level) -> OpPtr
{
    Assert(false,ExcNotImplemented());
    AssertThrow(false,ExcNotImplemented());

    OpPtr preconditioner;
    return preconditioner;
}
#endif


template class TrilinosTools<TrilinosImpl::tpetra>;
template class TrilinosTools<TrilinosImpl::epetra>;

#endif // #ifdef USE_TRILINOS

IGA_NAMESPACE_CLOSE


#endif
