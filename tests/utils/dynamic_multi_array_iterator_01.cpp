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
 * Test for DynamicMultiArray flat view
 * author: pauletti
 * data:   2013-06-15
 *
 */

#include "../tests.h"

#include <igatools/utils/dynamic_multi_array.h>


//Test the different constructors
template <int dim>
void do_test()
{
    out << "Default constructor "<< endl;
    DynamicMultiArray<Index, dim> data1;
    out << data1 << endl;


    out << "Uniform constructor "<< endl;
    DynamicMultiArray<Index, dim> data2(3);
    out << data2 << endl;


    out << "Rectangular constructor "<< endl;
    TensorSize<dim> size;
    for (int i = 0; i < dim; ++i)
        size(i) = 2+i;
    DynamicMultiArray<Index, dim> data3(size);
    out << data3 << endl;
}


//Test the different constructors
template <int dim>
void do_test_1()
{
    out << "Default constructor "<< endl;
    DynamicMultiArray<Index, dim> data1;
    out << data1 << endl;

    data1.resize(3);

    out << "Size: " << data1.flat_size() << endl;

    auto cit1_entry = data1.cbegin();
    const auto cit1_end = data1.cend();

    for (; cit1_entry != cit1_end; ++cit1_entry)
        out << *cit1_entry << " ";

    Index id = 0;
    for (auto & d : data1)
        d = id++;

    out << data1 << endl;


    TensorSize<dim> size;
    for (int i = 0; i < dim; ++i)
        size(i) = 2+i;
    data1.resize(size);
    out << "Size: " << data1.flat_size() << endl;

    id = 0;
    for (auto & d : data1)
        d = id++;

    out << data1 << endl;
}


template <int dim>
void do_test_2()
{
    out << "Default constructor "<< endl;
    DynamicMultiArray<Index, dim> data1(3);
    data1.fill_progression();
    out << "Size: " << data1.flat_size() << endl;
    out << data1 << endl;


    TensorSize<dim> size;
    for (int i = 0; i < dim; ++i)
        size(i) = 2+i;
    data1.resize(size);
    data1.fill_progression(10);
    out << "Size: " << data1.flat_size() << endl;
    out << data1 << endl;
}


//Test the extract flat view
template <int dim>
void do_test_3()
{
    out << "Default constructor "<< endl;
    DynamicMultiArray<Index, dim> data1(4);
    data1.fill_progression();
    out << "Size: " << data1.flat_size() << endl;
    out << data1 << endl;


    TensorIndex<dim> origin;
    for (int i = 0; i < dim; ++i)
        origin[i] = 1;
    TensorIndex<dim> increment;
    for (int i = 0; i < dim; ++i)
        increment[i] = 3-i;

    out << "Flat view: " << data1.get_flat_view(origin, increment) << endl;
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
