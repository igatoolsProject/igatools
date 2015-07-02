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
 *  Test for Vector class
 *
 *  author: pauletti
 *  date: Apr 5, 2013
 *
 */
//TODO (pauletti, Apr 4, 2015): this test is nonsense, rewrite or remove
#include "../tests.h"

#include <igatools/linear_algebra/epetra_vector.h>

using namespace EpetraTools;

void run_test()
{
    SafeSTLVector<Index> dofs_vec;
    Epetra_SerialComm comm;
    auto map = std::make_shared<Map>(-1, dofs_vec.size(), dofs_vec.data(), 0, comm);
    auto vec = create_vector(*map);

    vec->print_info(out);
    out << endl;



}

int main()
{
    out.begin_item("Testing Trilinos/Tpetra vector:");
    run_test();
    out.end_item();

    out.begin_item("Testing Trilinos/Epetra vector:");
    run_test();
    out.end_item();



    return  0;
}
