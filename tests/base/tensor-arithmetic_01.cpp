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
 *  Test for arithmetic tensor operators: +, -, *, /.
 *
 *  author: pauletti
 *  date: Feb 25, 2014
 *
 */

#include "../tests.h"
#include <igatools/base/tensor.h>

template <int dim, int range, int rank, int order>
void run_test()
{
    using Tensor = Derivatives<dim, range, rank, order>;
    Tensor A, B;
    Real alpha = 1.;

    out << "Testing Derivatives<" << dim << ",";
    out << range << "," << rank << "," << order << ">" << endl;

    out<< "A =" << endl;
    out << A << endl;

    out<< "B =" << endl;
    out << B << endl;

    out<< "alpha = " << alpha << endl;

    out<< "A + B =" << endl;
    out << A + B << endl;

    out<< "A - B =" << endl;
    out << A - B << endl;

    out<< "alpha * A =" << endl;
    out << alpha *A  << endl;

    out<< " A * alpha " << endl;
    out << A *alpha  << endl;

    out<< " A / alpha " << endl;
    out << A / alpha  << endl;
}


int main()
{
    out.depth_console(10);

    out << "Test for tensor arithmetic operators" << endl;
    run_test<2, 2, 1, 1>();
    run_test<3, 3, 1, 1>();

    return  0;
}
