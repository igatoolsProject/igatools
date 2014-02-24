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
 * Test for ProductArray
 * pauletti
 * 2013-01-15
 *
 */

#include "../tests.h"

#include <igatools/utils/product_array.h>
#include <igatools/utils/multi_array_utils.h>


//Test the different constructors
template <int dim>
void do_test()
{
    typedef ProductArray<Index, dim> ProdArr;

    ProdArr data1;
    out << data1 << endl;

    std::array< Index, dim> size;
    for (int i = 0; i < dim; ++i)
        size[i] = 2+i;

    ProdArr data2(size);
    out << data2 << endl;

    ProdArr data3(3);
    out << data3 << endl;
}


//Test the helper function to deal with indices
template <int dim>
void do_test1()
{
    std::array< Index, dim > size;
    for (int i = 0; i < dim; ++i)
        size[i] = 2+i;
    auto weight = MultiArrayUtils<dim>::compute_weight(size);
    out << "Size: " << endl;
    for (int i = 0; i < dim; ++i)
        out << size[i] << endl;

    out << "Tensor weight: " << endl;
    for (int i = 0; i < dim; ++i)
        out << weight[i] << endl;

    int flat_index = 10;
    auto tensor_index = MultiArrayUtils<dim>::flat_to_tensor_index(flat_index, weight);
    out <<  "Flat index: " << flat_index << "  ";
    out << "Tensor index: ( ";
    for (int i = 0; i < dim; ++i)
        out << tensor_index[i] << " ";
    out << ")" << endl;

    for (int i = 0; i < dim; ++i)
        tensor_index[i] = size[i]-2;
    out << "Tensor index: ( ";
    for (int i = 0; i < dim; ++i)
        out << tensor_index[i] << " ";
    out << ")" << "  ";
    flat_index = MultiArrayUtils<dim>::tensor_to_flat_index(tensor_index,weight);
    out <<  "Flat index: " << flat_index << endl;
}



int main(int argc, char *argv[])
{
    do_test<1>();
    do_test<2>();
    do_test<3>();

    do_test1<1>();
    do_test1<2>();
    do_test1<3>();
//  do_test1<4>();

    return 0;
}
