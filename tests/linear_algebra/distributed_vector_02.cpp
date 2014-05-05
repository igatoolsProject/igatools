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
 *  Test for Vector class
 *
 *  author: pauletti
 *  date: Apr 5, 2013
 *
 */

#include "../tests.h"

#include <igatools/linear_algebra/distributed_vector.h>


void run_test()
{
#if defined(USE_TRILINOS)
    const auto linear_algebra_package = LinearAlgebraPackage::trilinos;
#elif defined(USE_PETSC)
    const auto linear_algebra_package = LinearAlgebraPackage::petsc;
#endif
    using VectorType = Vector<linear_algebra_package>;


    /*
        VectorType v0(0);
        v0.print(out);
        out << endl;
    //*/
    VectorType v1(10);

    vector<Index> loc_dofs= {0,1,2};
    out << "local_coefs=" <<  v1.get_local_coefs(loc_dofs) << endl;
    v1.print(out);
    out << endl;

}

int main(int argc,char **argv)
{
    out.depth_console(10);

    PetscInitialize(&argc,&argv,(char *)0,"test");

    run_test();

    return  0;
}
