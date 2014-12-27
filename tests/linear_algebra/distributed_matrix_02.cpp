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
/*
 *  Test for definition of non-square Epetra Matrices using sparsity pattern
 *  created by using two different spaces.
 *  this test i going to be an adaptation of:
 *  matrix_definition01.cpp
 *  author: antolin
 *  date: 2013-04-10
 *
 */

// TODO (pauletti, Jun 3, 2014): the comment is not consistent with the test
// TODO (pauletti, Nov 24, 2014): the filename is not good
// TODO (pauletti, Nov 24, 2014):  this test needs to be rewritten

#include "../tests.h"

#include <igatools/basis_functions/bspline_space.h>
#include <igatools/basis_functions/bspline_element.h>
#include <igatools/geometry/unit_element.h>

#include <igatools/linear_algebra/sparsity_pattern.h>
#include <igatools/linear_algebra/linear_solver.h>
#include <igatools/linear_algebra/dof_tools.h>

using std::set;

int main(int argc, char *argv[])
{

    const int dim_domain = 1;
    const int dim_range  = 1;
    const int rank=1;
    const int p_r = 3;
    const int p_c = 2;

    out << " Domain dim: " << dim_domain;
    out << " Range dim: " << dim_range <<endl;
    out << " Degree of rows space: " << p_r <<endl;
    out << " Degree of columns space: " << p_c <<endl;


    int n_knots = 2;
    CartesianProductArray<Real,dim_domain> coord ;
    for (int i = 0; i < dim_domain; ++i)
    {
        vector<Real> tmp_coord;
        for (int j = 0; j < n_knots; ++j)
            tmp_coord.push_back(j);
        coord.copy_data_direction(i,tmp_coord);
    }



    auto knots = CartesianGrid<dim_domain>::create(coord);

    auto bspline_space_rows = BSplineSpace< dim_domain, dim_range, rank  >::create(p_r, knots) ;
    auto bspline_space_cols = BSplineSpace< dim_domain, dim_range, rank  >::create(p_c, knots) ;

    const auto dofs_view_row = bspline_space_rows->get_dof_distribution_global().get_dofs_view();
    const auto dofs_view_col = bspline_space_cols->get_dof_distribution_global().get_dofs_view();

    const auto n_basis_sp_rows = bspline_space_rows->get_num_basis();
    const auto n_basis_sp_cols = bspline_space_cols->get_num_basis();
    out << endl;
    out << "Number of dofs of rows space: " << n_basis_sp_rows << std::endl;
    out << "Number of dofs of columns space: " << n_basis_sp_cols << std::endl;
    out << endl;

    SpaceManager space_manager;
    space_manager.spaces_insertion_open();
    space_manager.add_space(bspline_space_rows);
    space_manager.add_space(bspline_space_cols);
    space_manager.spaces_insertion_close();

    space_manager.spaces_connectivity_open();
    space_manager.add_spaces_connection(bspline_space_rows,bspline_space_cols);
    space_manager.spaces_connectivity_close();

    auto &sp_conn = space_manager.get_spaces_connection(bspline_space_rows,bspline_space_cols);
    sp_conn.add_dofs_connectivity(
        dof_tools::build_dofs_connectvity_all_to_all(dofs_view_row,dofs_view_col));

#if defined(USE_TRILINOS)
    const auto la_pack = LAPack::trilinos;
#elif defined(USE_PETSC)
    const auto la_pack = LAPack::petsc;
#endif
    using VectorType = Vector<la_pack>;
    using MatrixType = Matrix<la_pack>;


//    SparsityPattern sparsity_pattern(space_manager);

    MatrixType A(space_manager);


    const auto num_rows = n_basis_sp_rows ;
    const auto num_cols = n_basis_sp_cols ;

    for (Index i = 0; i < num_cols ; i++)
        A.add_entry(dofs_view_row[i], dofs_view_col[i], 2.0);

    for (Index i = 0; i < num_rows - num_cols ; i++)
        A.add_entry(dofs_view_row[num_cols + i], dofs_view_col[i], 1.0);

    A.fill_complete();


    out << "A matrix" << endl;
    A.print_info(out);
    out << endl;


    VectorType b(vector<Index>(dofs_view_col.begin(),dofs_view_col.end()));
    Real val = 1.0;
    for (const Index id : dofs_view_col)
    {
        b.add_entry(id,val);
        val += 1;
    }

    out << "b vector" << endl;
    b.print_info(out);
    out << endl;

    VectorType c(vector<Index>(dofs_view_row.begin(),dofs_view_row.end()));


    // c = A . b
    A.multiply_by_right_vector(b,c,1.0,0.0);

    out << "c = A . b" << endl;
    c.print_info(out);

    out << endl;



    return 0;

}
