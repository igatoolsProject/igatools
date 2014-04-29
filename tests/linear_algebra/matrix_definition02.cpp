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
/*
 *  Test for definition of non-square Epetra Matrices using sparsity pattern
 *  created by using two different spaces.
 *  this test i going to be an adaptation of:
 *  matrix_definition01.cpp
 *  author: antolin
 *  date: 2013-04-10
 *
 */

#include "../tests.h"

#include <igatools/basis_functions/bspline_space.h>
#include <igatools/basis_functions/bspline_element_accessor.h>
#include <igatools/geometry/unit_element.h>

#include <igatools/linear_algebra/dof_tools.h>
#include <igatools/linear_algebra/linear_solver.h>

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
        vector<iga::Real> tmp_coord;
        for (int j = 0; j < n_knots; ++j)
            tmp_coord.push_back(j);
        coord.copy_data_direction(i,tmp_coord);
    }



    auto knots = CartesianGrid<dim_domain>::create(coord);

    auto bspline_space_rows = BSplineSpace< dim_domain, dim_range, rank  >::create (knots, p_r) ;
    auto bspline_space_cols = BSplineSpace< dim_domain, dim_range, rank  >::create (knots, p_c) ;

    const auto n_basis_sp_rows = bspline_space_rows->get_num_basis();
    const auto n_basis_sp_cols = bspline_space_cols->get_num_basis();
    out << endl;
    out << "Rows space:" << endl;
    bspline_space_rows->print_info(out);
    out << endl;
    out << "Columns space:" << endl;
    bspline_space_cols->print_info(out);
    out << endl;
    out << "Number of dofs of rows space: " << n_basis_sp_rows << std::endl;
    out << "Number of dofs of columns space: " << n_basis_sp_cols << std::endl;
    out << endl;


    Matrix A(dof_tools::get_sparsity_pattern<BSplineSpace< dim_domain, dim_range, rank  >,BSplineSpace< dim_domain, dim_range, rank  >>(bspline_space_rows, bspline_space_cols));


    const auto num_rows = n_basis_sp_rows ;
    const auto num_cols = n_basis_sp_cols ;

    for (Index i = 0; i < num_cols ; i++)
        A.add_entry(i, i, 2.0);

    for (Index i = 0; i < num_rows - num_cols ; i++)
        A.add_entry(num_cols + i, i, 1.0);

    A.fill_complete();


    out << "A matrix" << endl;
    A.print(out);
    out << endl;


    Vector b(bspline_space_cols->get_num_basis());
    for (Index i = 0; i < num_cols ; i++)
        b.add_entry(i,i + 1.0);

    out << "b vector" << endl;
    b.print(out);
    out << endl;

    Vector c(bspline_space_rows->get_num_basis());


    // c = A . b
    A.multiply_by_right_vector(b,c,1.0,0.0);

    out << "c = A . b" << endl;
    c.print(out);

    out << endl;



    return 0;

}
