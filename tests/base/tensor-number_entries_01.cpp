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
// Test for the action of a rank 2 Tensor
// Author: Seba 2012/10/26

#include "../tests.h"

#include <igatools/base/tensor.h>


template<int rdim, int cdim>
void do_number_entries_tensor()
{
  typedef Tensor<cdim, 1, tensor::covariant, Tensor< rdim, 1, tensor::contravariant, Tdouble> > T;

  out << "The number of entries of " <<
      "Tensor<" << cdim << ", 1, tensor::covariant, Tensor< " << rdim <<
      ", 1, tensor::contravariant, Tdouble> > is: " << T::get_number_of_entries() << std::endl;
}

template<int cdim>
void do_number_entries_vector()
{
  typedef Tensor<cdim, 1, tensor::contravariant, Tdouble> T;

  out << "The number of entries of " <<
      "Tensor<" << cdim << ", 1, tensor::contravariant, Tdouble> is: " <<
      T::get_number_of_entries() << std::endl;
}

void do_number_entries_scalar()
{
  out << "The number of entries of " <<
      "Tdouble is: " << Tdouble::get_number_of_entries() << std::endl;
}

int main(int argc, char *argv[])
{
  do_number_entries_scalar();

  do_number_entries_vector<1>();
  do_number_entries_vector<2>();
  do_number_entries_vector<3>();

  do_number_entries_tensor<1,1>();
  do_number_entries_tensor<2,3>();
  do_number_entries_tensor<2,2>();
  do_number_entries_tensor<3,2>();
  do_number_entries_tensor<2,3>();
  do_number_entries_tensor<3,3>();
}

