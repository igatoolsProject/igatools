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
 *  Tensor initialization
 *
 *  author: unclaimed
 *  date: 2014-11-23 (of notice)
 *
 */

//TODO: this test should be rewritten and improve to current standards
#include "../tests.h"
#include <igatools/base/tensor.h>

typedef Tensor< 1, 1, tensor::contravariant, Tdouble > ValueScalar;

int main()
{

    // testing the tensor.type for the scalar fields
    ValueScalar Value0;
    out << Value0 << endl;


    ValueScalar Value1;
    Value1 = 1.0;
    out << Value1 << endl;


    ValueScalar Value2;
    Value2[0] = 2.0;
    out << Value2 << endl;


    ValueScalar Value2a(Value2);
    out << Value2a << endl;


    ValueScalar Value2b = Value2;
    out << Value2b << endl;


    ValueScalar Value2c;
    Value2c = Value2;
    out << Value2c << endl;


    ValueScalar Value3;
    Value3 = 3.0;
    out << Value3 << endl;


    ValueScalar Value3a;
    Value3a = Value1 + Value2;
    out << Value3a << endl;

    ValueScalar Value4;
    Value4 = Value1 + Value3;
    out << Value4 << endl;


    ValueScalar Value4b;
    Value4b = Value3 + Value1;
    out << Value4b << endl;


    ValueScalar Value1a;
    Value1a = Value3 - Value2;
    out << Value1a << endl;


    ValueScalar Value4a;
    Value4a = Value3 + Value2 - Value1;
    out << Value4a << endl;


    ValueScalar Value3b;
    Value3b = Value4 - Value1;
    out << Value3b << endl;


    ValueScalar Value6;
    Value6 = Value3[0] * Value2[0];
    out << Value6 << endl;


    ValueScalar Value6a;
    Value6a = 3.0 * Value2;
    out << Value6a << endl;


    ValueScalar Value6b;
    Value6b = Value2 * 3.0;
    out << Value6b << endl;

    return 0;
}
