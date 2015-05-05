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
 * Test for iga::SafeSTLVector
 * author: pauletti
 * date:   2014-08-26
 *
 */

#include "../tests.h"

#include <igatools/utils/safe_stl_vector.h>
#include <igatools/utils/safe_stl_array.h>


#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>


using namespace iga;

void vector_print_info()
{
	OUTSTART

	SafeSTLVector<Real> vec(5,1);
    vec.print_info(out);
    out << endl;

	OUTEND
}


void array_print_info()
{
	OUTSTART

	SafeSTLArray<Real, 5> arr(2);
    arr.print_info(out);
    out<< endl;
    // SafeSTLArray<int, 3> a {1,2,3};
    // for (auto e : a)
    // out << e << endl;

	OUTEND
}


void vector_serialization()
{
	OUTSTART

	SafeSTLVector<Real> vec = {1.,2.,3.,4.,5.};
	vec.print_info(out);
	out << endl;


	boost::archive::text_oarchive out_archive(out.get_file_stream());

	out_archive << vec;

	OUTEND
}

int main()
{
    vector_print_info();
    array_print_info();

    vector_serialization();
    return 0;
}
