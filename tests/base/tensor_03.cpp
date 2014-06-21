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


#include "../tests.h"


#include <igatools/base/function.h>



int main(int argc, char *argv[])
{
    {
        typedef Function< 1, 1, 1 > Func ;

        out << "FunctionWithGradient< 1, 1, 0 >::Value" << endl ;
        Func::Value v ;
        out << v << endl ;

        v = 1.0 ;
        out << v << endl << endl ;

        out << "FunctionWithGradient< 1, 1, 0 >::Gradient" << endl ;
        Func::Gradient g ;
        out << g << endl ;

        g = 2.0 ;
        out << g << endl << endl ;

        out << "FunctionWithGradient< 1, 1, 0 >::Pt" << endl ;
        Func::PointType p ;
        out << p << endl ;

        p = 3.0 ;
        out << p << endl << endl ;
    }


    {
        typedef Function< 2, 1, 1 > Func ;

        out << "FunctionWithGradient< 2, 1, 0 >::Value" << endl ;
        Func::Value v ;
        out << v << endl ;

        v = 1.0 ;
        out << v << endl << endl ;

        out << "FunctionWithGradient< 2, 1, 0 >::Gradient" << endl ;
        Func::Gradient g ;
        out << g << endl ;

        g[0][0] = 2.0 ;
        g[1][0] = 3.0 ;
        out << g << endl << endl ;

        out << "FunctionWithGradient< 2, 1, 0 >::Pt" << endl ;
        Func::PointType p ;
        out << p << endl ;

        p[0] = 4.0 ;
        p[1] = 5.0 ;
        out << p << endl << endl ;
    }


    {
        typedef Function< 2, 2, 1 > Func ;

        out << "FunctionWithGradient< 2, 2, 1 >::Value" << endl ;
        Func::Value v ;
        out << v << endl ;

        v[0] = 1.0 ;
        v[1] = 2.0 ;
        out << v << endl << endl ;

        out << "FunctionWithGradient< 2, 2, 1 >::Gradient" << endl ;
        Func::Gradient g ;
        out << g << endl ;

        g[0][0] = 3.0 ;
        g[1][0] = 4.0 ;
        g[0][1] = 5.0 ;
        g[1][1] = 6.0 ;
        out << g << endl << endl ;

        out << "FunctionWithGradient< 2, 2, 1 >::Pt" << endl ;
        Func::PointType p ;
        out << p << endl ;

        p[0] = 7.0 ;
        p[1] = 8.0 ;
        out << p << endl << endl ;
    }

    {
        typedef Function< 2, 2, 2 > Func ;

        out << "FunctionWithGradient< 2, 2, 2 >::Value" << endl ;
        Func::Value v ;
        out << v << endl ;

        v[0][0] = 1.0 ;
        v[0][1] = 2.0 ;
        v[1][0] = 3.0 ;
        v[1][1] = 4.0 ;
        out << v << endl << endl ;

        out << "FunctionWithGradient< 2, 2, 2 >::Gradient" << endl ;
        Func::Gradient g ;
        out << g << endl ;

        g[0][0][0] = 5.0 ;
        g[1][0][0] = 6.0 ;
        g[0][0][1] = 7.0 ;
        g[1][0][1] = 8.0 ;
        g[0][1][0] = 9.0 ;
        g[1][1][0] =10.0 ;
        g[0][1][1] =11.0 ;
        g[1][1][1] =12.0 ;
        out << g << endl << endl ;

        out << "FunctionWithGradient< 2, 2, 2 >::Pt" << endl ;
        Func::PointType p ;
        out << p << endl ;

        p[0] = 13.0 ;
        p[1] = 14.0 ;
        out << p << endl << endl ;
    }

    {
        typedef Function< 3, 1, 1 > Func ;

        out << "FunctionWithGradient< 3, 1, 0 >::Value" << endl ;
        Func::Value v ;
        out << v << endl ;

        v = 1.0 ;
        out << v << endl << endl ;

        out << "FunctionWithGradient< 3, 1, 0 >::Gradient" << endl ;
        Func::Gradient g ;
        out << g << endl ;

        g[0][0] = 2.0 ;
        g[1][0] = 3.0 ;
        g[2][0] = 4.0 ;
        out << g << endl << endl ;

        out << "FunctionWithGradient< 3, 1, 0 >::Pt" << endl ;
        Func::PointType p ;
        out << p << endl ;

        p[0] = 5.0 ;
        p[1] = 6.0 ;
        p[2] = 7.0 ;
        out << p << endl << endl ;
    }


    {
        typedef Function< 3, 3, 1 > Func ;

        out << "FunctionWithGradient< 3, 3, 0 >::Value" << endl ;
        Func::Value v ;
        out << v << endl ;

        v[0] = 1.0 ;
        v[1] = 2.0 ;
        v[2] = 3.0 ;
        out << v << endl << endl ;

        out << "FunctionWithGradient< 3, 3, 0 >::Gradient" << endl ;
        Func::Gradient g ;
        out << g << endl ;

        g[0][0] = 4.0 ;
        g[1][0] = 5.0 ;
        g[2][0] = 6.0 ;
        g[0][1] = 7.0 ;
        g[1][1] = 8.0 ;
        g[2][1] = 9.0 ;
        g[0][2] =10.0 ;
        g[1][2] =11.0 ;
        g[2][2] =12.0 ;
        out << g << endl << endl ;

        out << "FunctionWithGradient< 3, 3, 0 >::Pt" << endl ;
        Func::PointType p ;
        out << p << endl ;

        p[0] = 13.0 ;
        p[1] = 14.0 ;
        p[2] = 15.0 ;
        out << p << endl << endl ;
    }


    {
        typedef Function< 3, 3, 2 > Func ;

        out << "FunctionWithGradient< 3, 3, 2 >::Value" << endl ;
        Func::Value v ;
        out << v << endl ;

        v[0][0] = 1.0 ;
        v[0][1] = 2.0 ;
        v[0][2] = 3.0 ;
        v[1][0] = 4.0 ;
        v[1][1] = 5.0 ;
        v[1][2] = 6.0 ;
        v[2][0] = 7.0 ;
        v[2][1] = 8.0 ;
        v[2][2] = 9.0 ;
        out << v << endl << endl ;

        out << "FunctionWithGradient< 3, 3, 2 >::Gradient" << endl ;
        Func::Gradient g ;
        out << g << endl ;

        g[0][0][0] = 10.0 ;
        g[1][0][0] = 11.0 ;
        g[2][0][0] = 12.0 ;
        g[0][0][1] = 13.0 ;
        g[1][0][1] = 14.0 ;
        g[2][0][1] = 15.0 ;
        g[0][0][2] = 16.0 ;
        g[1][0][2] = 17.0 ;
        g[2][0][2] = 18.0 ;
        g[0][1][0] = 19.0 ;
        g[1][1][0] = 20.0 ;
        g[2][1][0] = 21.0 ;
        g[0][1][1] = 22.0 ;
        g[1][1][1] = 23.0 ;
        g[2][1][1] = 24.0 ;
        g[0][1][2] = 25.0 ;
        g[1][1][2] = 26.0 ;
        g[2][1][2] = 27.0 ;
        g[0][2][0] = 28.0 ;
        g[1][2][0] = 29.0 ;
        g[2][2][0] = 30.0 ;
        g[0][2][1] = 31.0 ;
        g[1][2][1] = 32.0 ;
        g[2][2][1] = 33.0 ;
        g[0][2][2] = 34.0 ;
        g[1][2][2] = 35.0 ;
        g[2][2][2] = 36.0 ;
        out << g << endl << endl ;

        out << "FunctionWithGradient< 3, 3, 2 >::Pt" << endl ;
        Func::PointType p ;
        out << p << endl ;

        p[0] = 37.0 ;
        p[1] = 38.0 ;
        p[2] = 39.0 ;
        out << p << endl << endl ;
    }

    return (0) ;
}
