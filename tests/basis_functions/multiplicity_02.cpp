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
 *  Test for the multiplicity class
 *  author: pauletti
 *  date: Feb 19, 2013
 *
 */

#include "../tests.h"
#include <igatools/basis_functions/multiplicity.h>


//Test the different constructors
template<int dim, int range = 1, int rank = 1>
void do_test()
{
    using MultTable = Multiplicity<dim, range, rank>;
    using DegreeTable = typename MultTable::DegreeTable;
    auto grid = CartesianGrid<dim>::create();

    MultTable mul(data2, DegreeTable(TensorIndex<dim>(deg)));

//    Multiplicity<dim, range, rank> mult1;
//    mult1.print_info(out);
//    out << endl;

    const int n_knots = 3;
    const int deg     = 2;
    CartesianProductArray<Size, dim> cp(n_knots);
    typename MultTable::parent_t data2(cp);

    MultTable mult2(data2, DegreeTable(TensorIndex<dim>(deg)));

//    mult2.print_info(out);
    out << endl;
}
#if 0
    array<vector<int>,dim> mult_vector;
    for (int i = 0; i < dim; ++i)
    {
        mult_vector[i].resize(n_knots);
        for (int j = 0; j < n_knots; ++j)
            mult_vector[i][j] = i+j;
    }
    Multiplicity<dim> mult3(mult_vector);
    mult3.print_info(out);
    out << endl;


    TensorSize<dim> n_knots1;
    int_array<dim> deg1;
    for (int i = 0; i < dim; ++i)
    {
        n_knots1(i) = i+4;
        deg1[i] = i;
    }
    Multiplicity<dim> mult4(n_knots1);
    mult4.fill_max_regularity(deg1);
    mult4.print_info(out);
    out << endl;
}



template<int dim>
void do_test1()
{
    const int n_knots = 3;
    const int deg     = 2;
    Multiplicity<dim> mult2(n_knots);
    mult2.fill_max_regularity(deg);
    mult2.print_info(out);
    out << endl;

    auto cumul = mult2.accumulate();

    cumul.print_info(out);
    out << endl;
}
#endif

int main()
{
    out.depth_console(10);
    do_test<1>();
    do_test<2>();
    do_test<3>();


//    do_test1<1>();
//    do_test1<2>();
//    do_test1<3>();
    return 0;
}
