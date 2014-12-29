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
 * Test for CartesianProductIndexer
 * martinelli
 * 25 Feb 2014
 *
 */

#include "../tests.h"

#include <igatools/utils/cartesian_product_indexer.h>

template<int rank>
void
do_test()
{
    out << "========== BEGIN do_test<" << rank << "> ==========" << endl;
    TensorSize<rank> size;
    for (int i = 0 ; i < rank ; ++i)
        size[i] = i+2;

    CartesianProductIndexer<rank> indexer(size);
    out << "Num. indices = " << indexer.get_num_indices() << endl;
    indexer.print_info(out);
    out << endl;
    out << "========== END do_test<" << rank << "> ==========" << endl;
    out << endl;
}


int main()
{
    do_test<1>();
    do_test<2>();
    do_test<3>();

    return 0;
}
