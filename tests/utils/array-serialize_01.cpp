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

/**
 *  @file
 *  @brief  SafeSTLArray serialization
 *  @author  martinelli
 *  @date 2015-05-05
 */

#include "../tests.h"
#include <igatools/utils/safe_stl_array.h>

void array_serialization()
{
    OUTSTART

    SafeSTLArray<Real,3> arr = {-1.,-2.,-3.};

    {
        // writing to an xml file
        std::ofstream xml_ostream("array.xml");
        OArchive xml_out(xml_ostream);
        xml_out << BOOST_SERIALIZATION_NVP(arr);
        xml_ostream.close();
    }

    arr.fill(0);

    {
        // reading from an xml file
        std::ifstream xml_istream("array.xml");
        IArchive xml_in(xml_istream);
        xml_in >> BOOST_SERIALIZATION_NVP(arr);
        xml_istream.close();
    }


    out << endl;

    arr.print_info(out);
    out << endl;

    OUTEND
}


int main()
{
    array_serialization();
    return 0;
}
