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
 * Test for TensorProductArray
 * martinelli
 * 24 Feb 2014
 *
 */

#include "../tests.h"

#include <igatools/utils/tensor_product_array.h>


//Test the different constructors
template <int dim>
void do_test()
{
    out << "========= BEGIN do_test<" << dim << "> =========" << endl ;

    using ClassToTest = TensorProductArray<dim>;

    out<< "Testing the default constructor" << endl;
    ClassToTest data1;
    data1.print_info(out);
    out << endl;

    TensorSize<dim> size;
    for (int i = 0; i < dim; ++i)
        size(i) = 2+i;

    out << "Testing the entry() function"<<endl;
    ClassToTest data2(size);
    Index id = 1;
    for (int i = 0; i < dim ; ++i)
        for (int j = 0; j < size(i) ; ++j, ++id)
            data2.entry(i,j) = id*1.0;

    data2.print_info(out);
    out << endl;



    out<< "Testing the copy constructor" << endl;
    ClassToTest data3 = data2;
    data3.print_info(out);
    out << endl;

    out << "Testing the get_flat_cartesian_product() function" <<endl;
    vector<Point<dim>> flat_cartesian_product = data3.get_flat_cartesian_product();
    out << flat_cartesian_product << endl << endl;

    out << "Testing the get_flat_tensor_product() function" <<endl;
    vector<Real> flat_tensor_product = data3.get_flat_tensor_product();
    out << flat_tensor_product << endl << endl;

    out<< "Testing the get_sub_product() function" << endl;
    TensorIndex<dim-1> sub_id;
    for (int i = 0 ; i < dim-1 ; ++i)
        sub_id[i] = i+1;
    TensorProductArray<dim-1> data4 = data3.get_sub_product(sub_id);
    data4.print_info(out);
    out << endl;


    out << "Testing the insert() function" << endl;
    vector<Real> new_values = {-2.0,-1.0};
    TensorProductArray<dim> data5 = insert(data4, dim-1,new_values);
    data5.print_info(out);

    out << "========= END do_test<" << dim << "> =========" << endl ;
    out << endl;
}





int main(int argc, char *argv[])
{
    do_test<1>();
    do_test<2>();
    do_test<3>();


    return 0;
}
