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
#include <memory>

using std::vector;
using std::shared_ptr;


template <int dim>
void
do_test_1()
{
    out << "========== do_test_1() dim=" << dim << " --- begin ==========" << endl;
    using Grid = CartesianGrid<dim>;
    using RefSpace = BSplineSpace<dim>;
    using DMA = DynamicMultiArray<Index,dim>;
    using VecIt = typename DMA::const_iterator;
    using VecView = ConstView<VecIt>;

    std::vector<VecView> ranges;

    int n_spaces = 3;
    vector<shared_ptr<RefSpace>> ref_spaces(n_spaces);


    int dofs_offset = 0;
    for (int i_sp = 0 ; i_sp < n_spaces ; ++i_sp)
    {
        int n_knots = i_sp + 2;
        auto grid = Grid::create(n_knots);

        int degree = i_sp + 2;

        ref_spaces[i_sp] = RefSpace::create(grid,degree);

        ref_spaces[i_sp]->add_dofs_offset(dofs_offset);

        dofs_offset += ref_spaces[i_sp]->get_num_basis();

        const DMA &index_space = ref_spaces[i_sp]->get_index_space()(0);

        out << "Index space " << i_sp << " =" << endl;
        index_space.print_info(out);
        out << endl;

        ranges.push_back(
            VecView(index_space.begin(),index_space.end()));
    }


    ConcatenatedForwardConstIterator<VecIt> dofs_iterator_begin(ranges,0);
    ConcatenatedForwardConstIterator<VecIt> dofs_iterator_end(ranges,IteratorState::pass_the_end);

    out << "DOFs = [ ";
    for (; dofs_iterator_begin != dofs_iterator_end ; ++dofs_iterator_begin)
    {
        Index dof_id = *dofs_iterator_begin;
        out << dof_id << " ";
    }
    out << "]" << endl;

    out << "========== do_test_1() dim=" << dim << " --- end ==========" << endl;
    out << endl << endl;
}



template <int dim>
void
do_test_2()
{
    out << "========== do_test_2() dim=" << dim << " --- begin ==========" << endl;
    using Grid = CartesianGrid<dim>;
    using RefSpace = BSplineSpace<dim>;
    using DMA = DynamicMultiArray<Index,dim>;
    using VecIt = typename vector<Index>::const_iterator;
    using VecView = ConstView<VecIt>;

    std::vector<VecView> ranges;

    int n_spaces = 3;
    vector<shared_ptr<RefSpace>> ref_spaces(n_spaces);


    int dofs_offset = 0;
    for (int i_sp = 0 ; i_sp < n_spaces ; ++i_sp)
    {
        int n_knots = i_sp + 2;
        auto grid = Grid::create(n_knots);

        int degree = i_sp + 2;

        ref_spaces[i_sp] = RefSpace::create(grid,degree);

        ref_spaces[i_sp]->add_dofs_offset(dofs_offset);

        dofs_offset += ref_spaces[i_sp]->get_num_basis();

        const DMA &index_space = ref_spaces[i_sp]->get_index_space()(0);

        out << "Index space " << i_sp << " =" << endl;
        index_space.print_info(out);
        out << endl;


        ranges.push_back(
            VecView(index_space.get_data().begin(),index_space.get_data().end()));
    }


    ConcatenatedForwardConstIterator<VecIt> dofs_iterator_begin(ranges,0);
    ConcatenatedForwardConstIterator<VecIt> dofs_iterator_end(ranges,IteratorState::pass_the_end);

    out << "DOFs = [ ";
    for (; dofs_iterator_begin != dofs_iterator_end ; ++dofs_iterator_begin)
    {
        Index dof_id = *dofs_iterator_begin;
        out << dof_id << " ";
    }
    out << "]" << endl;

    out << "========== do_test_2() dim=" << dim << " --- end ==========" << endl;
    out << endl << endl;
}

template <int dim>
void do_tests()
{
    do_test_1<dim>();
    do_test_2<dim>();
}

int main()
{
    do_tests<1>();
    do_tests<2>();
    do_tests<3>();
}
