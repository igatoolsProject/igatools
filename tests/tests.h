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

using namespace std;
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
};

ofstream file("output.txt");
IgaTestOutput out(file);
