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
 * Test for StaticMultiArray: construction, index accessing
 * author: martinelli
 * data:   24 Feb 2014
 *
 */

#include "../tests.h"

#include <igatools/utils/static_multi_array.h>
// TODO (pauletti, Aug 28, 2014): this test is too long, split in several

//Test the different constructors
template <int dim>
void do_test()
{
    out << "========== BEGIN do_test<" << dim << "> ==========" << endl;

//    out << "Default constructor "<< endl;
    StaticMultiArray<Index,dim,2> data1;
//    out << endl;


    out << "Constant value constructor (value=3)"<< endl;
    StaticMultiArray<Index,dim,2> data2(3);
    data2.print_info(out);
    out << endl;


    out << "Variable value constructor (value=dim)"<< endl;
    array<Index,dim> val;
    for (int i = 0; i < dim; ++i)
        val[i] = i;
    StaticMultiArray<Index,dim,2> data3(val);
    data3.print_info(out);
    out << endl;
    out << "========== END do_test<" << dim << "> ==========" << endl;
    out << endl;
}


//Test the different constructors
template <int dim>
void do_test_1()
{
    out << "========== BEGIN do_test_1<" << dim << "> ==========" << endl;
    StaticMultiArray<Index,dim,2> data1;

    out << "Size: " << data1.flat_size() << endl;


    for (int i = 0; i < data1.flat_size(); ++i)
        data1(i) = i;

    data1.print_info(out);
    out << endl;


    out << "Size: " << data1.flat_size() << endl;
    for (int i = 0; i < data1.flat_size(); ++i)
        data1(i) = i;

    data1.print_info(out);
    out << endl;
    out << "========== END do_test_1<" << dim << "> ==========" << endl;
    out << endl;
}


template <int dim>
void do_test_2()
{
    out << "========== BEGIN do_test_2<" << dim << "> ==========" << endl;
    out << "Fill progression from 0 "<< endl;
    StaticMultiArray<Index,dim,2> data1(3);
    data1.fill_progression();
    out << "Size: " << data1.flat_size() << endl;
    data1.print_info(out);
    out << endl;


    out << "Fill progression from 10 "<< endl;
    data1.fill_progression(10);
    out << "Size: " << data1.flat_size() << endl;
    data1.print_info(out);
    out << endl;
    out << "========== END do_test_2<" << dim << "> ==========" << endl;
    out << endl;
}


template <int dim>
void do_test_3()
{
    out << "========== BEGIN do_test_3<" << dim << "> ==========" << endl;
    out << "Constant value constructor (value=4)"<< endl;
    StaticMultiArray<Index,dim,2> data1(4);
    out << "Fill progression from 2 "<< endl;
    data1.fill_progression(2);
    out << "Size: " << data1.flat_size() << endl;
    data1.print_info(out);
    out << endl;
    out << "========== END do_test_3<" << dim << "> ==========" << endl;
    out << endl;
}


int main()
{
    out.depth_console(10);
    do_test<1>();
    do_test<2>();
    do_test<3>();

    do_test_1<1>();
    do_test_1<2>();
    do_test_1<3>();

    do_test_2<1>();
    do_test_2<2>();
    do_test_2<3>();

    do_test_3<1>();
    do_test_3<2>();
    do_test_3<3>();

    return 0;
}
