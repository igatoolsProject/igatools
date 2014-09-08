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
 * Test for iga::vector
 * author: pauletti
 * date:   2014-08-26
 *
 */

#include "../tests.h"

#include <igatools/utils/vector.h>
#include <igatools/utils/array.h>

void vector_print_info()
{
    iga::vector<Real> vec(5);
    vec.print_info(out);
    out << endl;
}


void array_print_info()
{
    special_array<Real, 5> arr;
    arr.print_info(out);
    out<< endl;
    special_array<int, 3> a{1,2,3};
    for (auto e : a)
    	out << e << endl;

}


int main()
{
    vector_print_info();
    array_print_info();
    return 0;
}
