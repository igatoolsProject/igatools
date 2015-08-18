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
 *  @brief Tensor range generator function
 *  @author pauletti
 *  @date 2015-03-07
 */


#include "../tests.h"
#include <igatools/utils/tensor_range.h>

void get_range()
{
    OUTSTART
    {
        TensorIndex<1> first{3};
        TensorIndex<1> last{7};
        out.begin_item("Tensor range 1D:");
        el_tensor_range(first, last).print_info(out);
        out.end_item();
    }

    {
        TensorIndex<2> first{3,5};
        TensorIndex<2> last{7,10};
        out.begin_item("Tensor range 2D:");
        el_tensor_range(first, last).print_info(out);
        out.end_item();
    }

    {
        TensorIndex<3> first{3,5,1};
        TensorIndex<3> last{7,10, 3};
        out.begin_item("Tensor range 3D:");
        el_tensor_range(first, last).print_info(out);
        out.end_item();
    }
    OUTEND
}


int main()
{
    get_range();
    return 0;
}
