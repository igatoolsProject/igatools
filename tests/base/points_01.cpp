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
 * Test for the Points alias
 * author: pauletti
 * date: Jun 21, 2014
 *
 */

#include "../tests.h"
#include "igatools/base/tensor.h"

template<int dim>
void default_constructor()
{
    Points<dim> p;
    out << p << endl;
}

void init_list()
{
    Points<2> p2 = {1., 2.5};
    Points<3> p3( {1., 2.5, -4.3});
    Points<3> p3a {1., 2.5, -4.3};
    out << p2 << p3 << p3a  << endl;
}

int main()
{

    default_constructor<0>();
    default_constructor<1>();
    default_constructor<2>();
    default_constructor<3>();

    init_list();

    return 0;
}
