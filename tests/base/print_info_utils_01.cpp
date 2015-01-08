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
 *  Test for
 *
 *  author: pauletti
 *  date:2014-08-26
 */

#include "../tests.h"
#include <igatools/base/print_info_utils.h>


class A
{
public:
    void print_info(LogStream &) {}
};


class B
{};


template <class Z>
EnableIf<has_print_info<Z>(0), void>
print_info(LogStream &out)
{
}

template <class A>
EnableIf<(!has_print_info<A>(0)), void>
print_info(LogStream &out)
{
}

int main()
{
    out << has_print_info<A>(0) << endl;
    out << has_print_info<B>(0) << endl;
    print_info<A>(out);
    print_info<B>(out);
    return 0;
}
