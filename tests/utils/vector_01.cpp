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
 * Test for SafeSTLVector, SafeSTLArray and their serialization
 * author: pauletti, martinelli
 * date:   2014-08-26
 * date:   2015-05-05 (added the tests for serialization)
 *
 */

#include "../tests.h"

#include <igatools/utils/safe_stl_vector.h>
#include <igatools/utils/safe_stl_array.h>

#include <fstream>

//#include <boost/archive/text_oarchive.hpp>
//#include <boost/archive/text_iarchive.hpp>


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

    {
        // writing to an xml file
        std::ofstream xml_ostream("vector.xml");
        boost::archive::xml_oarchive xml_out(xml_ostream);
        xml_out << BOOST_SERIALIZATION_NVP(vec);
        xml_ostream.close();
    }

    vec.clear();

    {
        // reading from an xml file
        std::ifstream xml_istream("vector.xml");
        boost::archive::xml_iarchive xml_in(xml_istream);
        xml_in >> BOOST_SERIALIZATION_NVP(vec);
        xml_istream.close();
    }

    for (auto &v : vec)
        v += 1.0;

    {
        boost::archive::xml_oarchive xml_out_1(out.get_file_stream());
        xml_out_1 << BOOST_SERIALIZATION_NVP(vec);
    }

    /*
    {
        boost::archive::text_oarchive text_out(out.get_file_stream());
        text_out << vec;
    }
    //*/

    out << endl;

    vec.print_info(out);
    out << endl;

    OUTEND
}

void array_serialization()
{
    OUTSTART

    SafeSTLArray<Real,3> arr = {-1.,-2.,-3.};

    {
        // writing to an xml file
        std::ofstream xml_ostream("array.xml");
        boost::archive::xml_oarchive xml_out(xml_ostream);
        xml_out << BOOST_SERIALIZATION_NVP(arr);
        xml_ostream.close();
    }

    arr.fill(0);

    {
        // reading from an xml file
        std::ifstream xml_istream("array.xml");
        boost::archive::xml_iarchive xml_in(xml_istream);
        xml_in >> BOOST_SERIALIZATION_NVP(arr);
        xml_istream.close();
    }

    for (auto &v : arr)
        v += 1.0;

    {
        boost::archive::xml_oarchive xml_out_1(out.get_file_stream());
        xml_out_1 << BOOST_SERIALIZATION_NVP(arr);
    }

    /*
        {
            boost::archive::text_oarchive text_out(out.get_file_stream());
            text_out << arr;
        }
    //*/

    out << endl;

    arr.print_info(out);
    out << endl;

    OUTEND
}

int main()
{
    vector_print_info();
    array_print_info();

    vector_serialization();
    array_serialization();
    return 0;
}
