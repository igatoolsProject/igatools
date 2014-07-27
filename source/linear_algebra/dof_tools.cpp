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


#include <igatools/linear_algebra/dof_tools.h>
#include <igatools/base/exceptions.h>
#include <igatools/linear_algebra/distributed_matrix.h>


using std::vector;
using std::map;
using std::set;
using std::pair;
using std::shared_ptr;

IGA_NAMESPACE_OPEN

namespace dof_tools
{
SparsityPattern
get_sparsity_pattern(const DofsManager &dofs_manager)
{
    Assert(dofs_manager.is_dofs_view_open() == false,ExcInvalidState());

    // build the dofs graph
    const auto & dofs_view = dofs_manager.get_dofs_view();

    vector<Index> dofs_copy;
    for (const auto &dof : dofs_view)
    	dofs_copy.push_back(dof);

    Assert(!dofs_copy.empty(),ExcEmptyObject());

    SparsityPattern sparsity_pattern(dofs_copy, dofs_copy);

    using DofsInRow = set<Index>;
    DofsInRow empty_set;

    // adding the global dof keys to the map representing the dof connectivity
    for (const auto &dof : dofs_copy)
        sparsity_pattern.insert(pair<Index,DofsInRow>(dof,empty_set));

    Assert(dofs_manager.are_elements_dofs_view_open() == false,ExcInvalidState());
    Assert(!dofs_manager.get_elements_dofs_view().empty(),
    		ExcEmptyObject());
    for (const auto element_dofs : dofs_manager.get_elements_dofs_view())
    	for (const auto &dof : element_dofs)
    		sparsity_pattern[dof].insert(element_dofs.begin(),element_dofs.end());

    return (sparsity_pattern);
}
SparsityPattern
get_sparsity_pattern(const DofsManager &dofs_manager_rows,const DofsManager &dofs_manager_cols)
{
    Assert(dofs_manager_rows.is_dofs_view_open() == false,ExcInvalidState());
    Assert(dofs_manager_cols.is_dofs_view_open() == false,ExcInvalidState());

    // build the dofs graph
    const auto & dofs_view_rows = dofs_manager_rows.get_dofs_view();
    const auto & dofs_view_cols = dofs_manager_cols.get_dofs_view();

    vector<Index> dofs_copy_rows;
    for (const auto &dof : dofs_view_rows)
    	dofs_copy_rows.push_back(dof);
    Assert(!dofs_copy_rows.empty(),ExcEmptyObject());

    vector<Index> dofs_copy_cols;
    for (const auto &dof : dofs_view_cols)
    	dofs_copy_cols.push_back(dof);
    Assert(!dofs_copy_cols.empty(),ExcEmptyObject());

    SparsityPattern sparsity_pattern(dofs_copy_rows, dofs_copy_cols);

    using DofsInRow = set<Index>;
    DofsInRow empty_set;

    // adding the global dof keys to the map representing the dof connectivity
    for (const auto &dof : dofs_copy_rows)
        sparsity_pattern.insert(pair<Index,DofsInRow>(dof,empty_set));


    Assert(dofs_manager_rows.are_elements_dofs_view_open() == false,ExcInvalidState());
    const auto & elements_dofs_rows = dofs_manager_rows.get_elements_dofs_view();
    Assert(!elements_dofs_rows.empty(),ExcEmptyObject());

    Assert(dofs_manager_cols.are_elements_dofs_view_open() == false,ExcInvalidState());
    const auto & elements_dofs_cols = dofs_manager_cols.get_elements_dofs_view();
    Assert(!elements_dofs_cols.empty(),ExcEmptyObject());

    Assert(elements_dofs_rows.size() == elements_dofs_cols.size(),
    		ExcDimensionMismatch(elements_dofs_rows.size(),elements_dofs_cols.size()));

    const Index n_elements = elements_dofs_rows.size();
    for (Index ielem = 0 ; ielem < n_elements ; ++ielem)
    {
    	const auto &dofs_rows = elements_dofs_rows[ielem];
    	const auto &dofs_cols = elements_dofs_cols[ielem];

    	for (const auto & dof_row : dofs_rows)
    		sparsity_pattern[dof_row].insert(dofs_cols.begin(),dofs_cols.end());
    }
    return (sparsity_pattern);
}


#if 0
template < class SpaceType >
vector<Index>
get_dofs(shared_ptr<const SpaceType> space, EnableIf<is_function_space<SpaceType>()> *)
{
//  LogStream out;
//  space->print_info(out);
    auto element = space->begin();
    const auto element_end = space->end();

    set<Index> dofs_set;

    for (; element != element_end; ++element)
    {
        const vector< Index > element_dofs = element->get_local_to_global();
        Assert(!element_dofs.empty(),ExcEmptyObject());

        for (const Index &dof : element_dofs)
            dofs_set.insert(dof);
    }

    vector<Index> space_dofs(dofs_set.begin(), dofs_set.end());

    return space_dofs;
}
#endif


template <LAPack la_pack>
void apply_boundary_values(const std::map<Index,Real> &boundary_values,
                           Matrix<la_pack> &matrix,
                           Vector<la_pack> &rhs,
                           Vector<la_pack> &solution)
{
    std::vector<Index> constrained_rows;

    auto dof = boundary_values.begin();
    const auto dof_end = boundary_values.end();
    for (; dof != dof_end; ++dof)
    {
        Index row_id = dof->first;
        const Real bc_value  = dof->second;
        const Real mat_value = matrix(row_id,row_id);


        // set the matrix in write mode
        matrix.resume_fill();

        // set the selected row to 0.0
        matrix.clear_row(row_id);

        // set the diagonal element corresponding to the entry
        // (row_id,row_id) to mat_value
        matrix.add_entry(row_id, row_id, mat_value);

        // communicate the matrix values to the different processors
        matrix.fill_complete();

        rhs(row_id) = bc_value * mat_value;
        solution(row_id) = bc_value;
    }
}



#ifdef USE_PETSC

template <>
void apply_boundary_values(const std::map<Index,Real> &boundary_values,
                           Matrix<LAPack::petsc> &matrix,
                           Vector<LAPack::petsc> &rhs,
                           Vector<LAPack::petsc> &solution)
{
    PetscErrorCode ierr;

    vector<Index> rows;
    vector<PetscScalar> values;

    for (const auto &bv : boundary_values)
    {
        rows.push_back(bv.first);
        values.push_back(bv.second);
    }

    // Set the matrix in write mode.
    matrix.resume_fill();

    // Set the boundary value in the solution vector.
    ierr = VecSetValues(solution.get_petsc_vector(), rows.size(), rows.data(),
                        values.data(), INSERT_VALUES);

    // Getting the first diagonal value (of the constrained degrees of freedom),
    // this value is going to be written in the diagonal for all the removed
    // rows/columns.
    // Note: It would be desirable to keep the correspoding values for
    // every row/column, but the petsc function MatZeroRowsColumns only allows
    // to specify a single value for the diagonal terms.
    // This could be changed with multiple calls to this function (one for every
    // degree of freedom constrained). Probably it would be expensive, but the
    // condition number would be possibly improved.
    const Real diagonal = matrix(rows[0], rows[0]);

    // Setting to zero the rows and columns corresponding to rows.
    // The diagonal terms will be set with the value of diagonal.
    // The corresponding values or rhs will be set to be equal to
    // diagonal * values (for every vector element).
    ierr = MatZeroRowsColumns(matrix.get_petsc_matrix(), rows.size(),
                              rows.data(), diagonal,
                              solution.get_petsc_vector(),
                              rhs.get_petsc_vector());

    // Communicate the matrix values to the different processors.
    matrix.fill_complete();

}
#endif //#ifdef USE_PETSC

};

IGA_NAMESPACE_CLOSE

#include <igatools/linear_algebra/dof_tools.inst>
