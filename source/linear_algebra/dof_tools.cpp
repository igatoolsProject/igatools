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


using std::map;
using std::set;
using std::pair;
using std::shared_ptr;

IGA_NAMESPACE_OPEN

namespace dof_tools
{

template <LAPack la_pack>
void apply_boundary_values(const std::map<Index,Real> &boundary_values,
                           Matrix<la_pack> &matrix,
                           Vector<la_pack> &rhs,
                           Vector<la_pack> &solution)
{
//    vector<Index> constrained_rows;
;
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

}

IGA_NAMESPACE_CLOSE

#include <igatools/linear_algebra/dof_tools.inst>
