//-+--------------------------------------------------------------------
// Igatools a general purpose Isogeometric analysis library.
// Copyright (C) 2012-2016  by the igatools authors (see authors.txt).
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




# include <igatools/base/tensor.h>


#define dim_domain 3
#define dim_range  3

typedef Tensor< dim_domain, 1, tensor::covariant, Tensor< 1, 1, tensor::contravariant, Tdouble > >
GradientFunctionScalar ;

int main(int argc, char *argv[])
{

  GradientFunctionScalar g1 ;
  out << g1 << endl ;

  g1[0][0] = 0.0 ;
  g1[1][0] = 1.0 ;
  g1[2][0] = 2.0 ;

  out << g1 << endl ;

  GradientFunctionScalar g1a ;
  g1a[0][0] -= g1[0][0] ;
  g1a[1][0] -= g1[1][0] ;
  g1a[2][0] -= g1[2][0] ;
  out << g1a << endl ;

  GradientFunctionScalar g2 ;
  g2[0][0] = 2.0 ;
  g2[1][0] = 1.0 ;
  g2[2][0] = 0.0 ;

  out << g2 << endl ;


  GradientFunctionScalar g3 ;
  g3 = g2 ;
  out << g3 << endl ;


  GradientFunctionScalar g4 = g2 ;
  out << g4 << endl ;

  GradientFunctionScalar g5(g2) ;
  out << g5 << endl ;


  GradientFunctionScalar g6 ;
  g6 = g1 + g2 ;
  out << g6 << endl ;

  GradientFunctionScalar g7 ;
  g7 = g1 - g2 ;
  out << g7 << endl ;

//  double dotprod = DotProduct( g1, g2 ) ;
//  out << dotprod << endl ;


  return (0) ;
}
