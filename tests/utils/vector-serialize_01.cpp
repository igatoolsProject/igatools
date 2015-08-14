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
 *  @brief  SafeSTLVector serialization
 *  @author  martinelli
 *  @date 2015-05-05
 */

#include "../tests.h"

#include <igatools/utils/safe_stl_vector.h>

using namespace iga;

void vector_print_info()
{
    OUTSTART

    SafeSTLVector<Real> vec(5,1);
    vec.print_info(out);
    out << endl;

    OUTEND
}



void vector_serialization()
{
    OUTSTART

    SafeSTLVector<Real> vec = {1.,2.,3.,4.,5.};

    {
        // writing to an xml file
        std::ofstream xml_ostream("vector.xml");
        OArchive xml_out(xml_ostream);
        xml_out << BOOST_SERIALIZATION_NVP(vec);
        xml_ostream.close();
    }

    vec.clear();

    {
        // reading from an xml file
        std::ifstream xml_istream("vector.xml");
        IArchive xml_in(xml_istream);
        xml_in >> BOOST_SERIALIZATION_NVP(vec);
        xml_istream.close();
    }


    out << endl;

    vec.print_info(out);
    out << endl;

    OUTEND
}



int main()
{
    vector_serialization();
    return 0;
}
