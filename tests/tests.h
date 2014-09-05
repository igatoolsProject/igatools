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
 *  Common header file for all tests
 */

#include <igatools/base/logstream.h>
#include <fstream>


using std::array;
using std::shared_ptr;
using std::make_shared;
using std::const_pointer_cast;
using std::dynamic_pointer_cast;
using std::string;
using std::to_string;
using std::endl;
using std::string;
using std::ofstream;

//using std::endl;
//using std::array;
//using std::string;
//using std::ofstream;

using namespace iga;

class IgaTestOutput : public LogStream
{
public:
    IgaTestOutput(std::ostream &o)
        :
        LogStream()
    {
        LogStream::test_mode();
        LogStream::attach(o);
        LogStream::depth_console(0);
    }

    void line(string title = "")
    {
        (*this) << "========================================================================" << endl;
        (*this) << title << endl;
    }
};

std::ofstream file("output.txt");
IgaTestOutput out(file);

#define SEPARATOR "========================================================================"
#define OUTSTART out << SEPARATOR << endl << __PRETTY_FUNCTION__ << endl << SEPARATOR << endl;
#define OUTEND out << SEPARATOR << endl << endl;
