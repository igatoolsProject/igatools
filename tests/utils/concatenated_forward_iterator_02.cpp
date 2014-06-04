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
 *
 * Test for the Concatenated iterator class, using as base iterator the
 * std::vector<int>::iterator.
 *
 * martinelli
 * 30 May 2014
 *
 */

#include "../tests.h"

#include <igatools/basis_functions/bspline_space.h>
#include <igatools/utils/concatenated_forward_iterator.h>

#include <vector>

using std::vector;


template <int dim>
void
do_test()
{
    out << "========== do_test() dim=" << dim << " --- begin ==========" << endl;
    using Grid = CartesianGrid<dim>;
    using RefSpace = BSplineSpace<dim>;

    auto grid_1 = Grid::create(2);
    auto grid_2 = Grid::create(3);


    auto ref_space_1 = RefSpace::create(grid_1,2);
    const int n_dofs_space_1 = ref_space_1->get_num_basis();

    auto ref_space_2 = RefSpace::create(grid_2,3);
    ref_space_2->add_dofs_offset(n_dofs_space_1);

    ref_space_1->print_info(out);
    ref_space_2->print_info(out);

    using DMA = DynamicMultiArray<Index,dim>;
    const DMA &index_space_1 = ref_space_1->get_index_space()(0);
    const DMA &index_space_2 = ref_space_2->get_index_space()(0);

    out << "Index space 1 =" << endl;
    index_space_1.print_info(out);
    out << endl;

    out << "Index space 2 =" << endl;
    index_space_2.print_info(out);
    out << endl;
    using VecIt = typename vector<Index>::const_iterator;
    VecIt index_space_1_begin = index_space_1.get_data().begin();
    VecIt index_space_1_end = index_space_1.get_data().end();

    VecIt index_space_2_begin = index_space_2.get_data().begin();
    VecIt index_space_2_end = index_space_2.get_data().end();

    using PairVecIt = std::pair<VecIt,VecIt>;
    std::vector<PairVecIt> ranges;
    ranges.push_back(PairVecIt(index_space_1_begin,index_space_1_end));
    ranges.push_back(PairVecIt(index_space_2_begin,index_space_2_end));

    ConcatenatedForwardIterator<VecIt> dofs_iterator_begin(ranges,0);
    ConcatenatedForwardIterator<VecIt> dofs_iterator_end(ranges,IteratorState::pass_the_end);

    out << "DOFs = [ ";
    for (; dofs_iterator_begin != dofs_iterator_end ; ++dofs_iterator_begin)
    {
        Index dof_id = *dofs_iterator_begin;
        out << dof_id << " ";
    }
    out << "]" << endl;

    out << "========== do_test() dim=" << dim << " --- end ==========" << endl;
    out << endl;
}

int main()
{
    vector<int> v0 = {1,2,3,4};
    vector<int> v1 = {5,6,7,8,9};
    vector<int> v2 = {10,11,12};
    vector<int> v3 = {13};

    using VecIterator = vector<int>::iterator;

    std::vector<std::pair<VecIterator,VecIterator>> ranges;
    ranges.push_back(std::make_pair<VecIterator,VecIterator>(v0.begin(),v0.end()));
    ranges.push_back(std::make_pair<VecIterator,VecIterator>(v1.begin(),v1.end()));
    ranges.push_back(std::make_pair<VecIterator,VecIterator>(v2.begin(),v2.end()));
    ranges.push_back(std::make_pair<VecIterator,VecIterator>(v3.begin(),v3.end()));


    ConcatenatedForwardIterator<VecIterator> begin(ranges,0);
    ConcatenatedForwardIterator<VecIterator> end(ranges,IteratorState::pass_the_end);



    using std::endl;
    int i = 0;
    for (; begin != end ; ++begin, ++i)
        out << "i = " << i << "     value = " << *begin << std::endl;


    do_test<1>();
//    do_test<2>();
//    do_test<3>();
}
