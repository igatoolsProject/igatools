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
 *  Test for the SphericalFunction class as a mapping
 *
 *  author: pauletti
 *  date: 2014-10-24
 *
 */

#include "../tests.h"

#include <igatools/linear_algebra/dense_matrix.h>


template <int dim>
void comp_determinant()
{
    OUTSTART

    DenseMatrix A(dim, dim);
    for (int i=0; i<dim; ++i)
        for (int j=0; j<dim; ++j)
            A(i,j) = pow(i,j);

    A.print_info(out);
    out << endl << "Determinant: " << A.determinant() << endl;

    OUTEND
}


int main()
{
    comp_determinant<2>();
    comp_determinant<3>();

    return 0;
}
