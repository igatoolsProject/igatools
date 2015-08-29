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
 * Test for ValueVector
 * martinelli
 * 2013-06-12
 *
 */

#include "../tests.h"

#include <igatools/utils/value_vector.h>

//Test the different constructors
void do_test()
{
  out << "---------------------------------------" << endl ;
  out << "Test the different constructors" << endl << endl ;

  typedef ValueVector<Real> ClassToTest ;


  out << "Default constructor" << endl ;
  ClassToTest data1;
  out << "data1:\n";
  data1.print_info(out);
  out << endl ;


  out << "Constructor with the size of the vector" << endl ;
  const int size = 3 ;
  ClassToTest data2(size);
  out << "data2:\n";
  data2.print_info(out);
  out << endl ;


  out << "Initializer-list constructor" << endl ;
  SafeSTLVector<Real> vec = {1.0,2.0,3.0};
  ClassToTest data3(vec);
  out << "data3:\n";
  data3.print_info(out);
  out << endl ;


  out << "Copy constructor" << endl ;
  ClassToTest data4 = data3 ;
  out << "data4:\n";
  data4.print_info(out);
  out << endl;
  out << "data3:\n";
  data3.print_info(out);
  out << endl ;


  out << "Move constructor" << endl ;
  ClassToTest data5 = std::move(data3);
  out << "data5:\n";
  data5.print_info(out);
  out << endl;

  out << "---------------------------------------" << endl ;

  out << endl ;

}

//Test the copy/move assignment operators
void do_test1()
{
  out << "---------------------------------------" << endl ;
  out << "Test the copy/move assignment operators" << endl << endl ;

  typedef ValueVector<Real> ClassToTest ;

  ClassToTest data1 = {1.0,2.0,3.0};

  out << "Copy assignment" << endl ;
  ClassToTest data2 ;
  out << "before copy" << endl ;
  out.push("\t") ;
  out << "data1:\n";
  data1.print_info(out);
  out << endl ;
  out << "data2:\n";
  data2.print_info(out);
  out << endl ;
  out.pop() ;

  data2 = data1 ;
  out << "after copy" << endl ;
  out.push("\t") ;
  out << "data1:\n";
  data1.print_info(out);
  out << endl;
  out << "data2:\n";
  data2.print_info(out);
  out << endl;

  out.pop() ;
  out << endl ;

  out << "Move assignment" << endl ;
  ClassToTest data3 ;
  out << "before move" << endl ;
  out.push("\t") ;
  out << "data1:\n";
  data1.print_info(out);
  out << endl;
  out << "data3:\n";
  data3.print_info(out);
  out << endl;
  out.pop() ;

  data3 = std::move(data1) ;
  out << "after move" << endl ;
  out.push("\t") ;
//    data1.print_info(out);
  out << "data3:\n";
  data3.print_info(out);
  out << endl;
  out.pop() ;


  out << "---------------------------------------" << std::endl ;

  out << endl ;

}

//Test the scalar-by-SafeSTLVector multiplication
void do_test2()
{
  out << "---------------------------------------" << endl ;
  out << "Test the scalar-by-vector multiplication" << endl << endl ;

  typedef ValueVector<Real> ClassToTest ;

  ClassToTest data = {1.0,2.0,3.0};
  out << "Original data =" << endl ;
  data.print_info(out);
  out << endl ;

  const double a = -1.0 ;
  ClassToTest data1 = a * data ;
  out << "-1.0 * data =" << endl ;
  data1.print_info(out);
  out << endl ;

  const double b = 2.0 ;
  ClassToTest data2 = data * b ;
  out << "data * 2.0 =" << endl ;
  data2.print_info(out);
  out << endl;
  out << "---------------------------------------" << endl ;
  out << endl ;
}


int main(int argc, char *argv[])
{
  do_test();
  do_test1();
  do_test2();

  return 0;
}
