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
 *  Test for the NURBS space iterator
 *
 *  author: pauletti
 *  date: Jun 11, 2014
 *
 */

#include "../tests.h"

#include <igatools/base/quadrature_lib.h>
#include <igatools/basis_functions/nurbs_space.h>
#include <igatools/basis_functions/nurbs_element_accessor.h>

template< int dim, int range, int rank = 1>
void test()
{
    const int r = 2;
    out << "test<" << dim << "," << range << ">" << endl;

    using Space = NURBSSpace< dim, range, rank >;
    using WeightsTable = typename Space::WeightsTable;
    using DegreeTable = typename Space::DegreeTable;
    auto  knots = CartesianGrid<dim>::create();

    auto degree = TensorIndex<dim>(r);
    DegreeTable deg(degree);

    auto  bsp = BSplineSpace<dim, range, rank >::create(deg, knots);
    WeightsTable weights;
    const auto n_basis = bsp->get_num_basis_table();
    for (auto comp : Space::components)
        weights(comp).resize(n_basis(comp),1.0);

    auto space = Space::create(deg, knots, weights);

    const int n_points = 3;
    QGauss<dim> quad(n_points);

    auto elem     = space->begin();
    auto end_element = space->end();

    const auto flag = ValueFlags::value|ValueFlags::gradient|ValueFlags::hessian;
    elem->init_cache(flag, quad);

    for (; elem != end_element; ++elem)
    {
        elem->fill_cache();
        out << "Element: " << elem->get_flat_index()<< endl;

        out << "Values basis functions:" << endl;
        auto values = elem->get_basis_values();
        values.print_info(out);

        out << "Gradients basis functions:" << endl;
        auto gradients = elem->get_basis_gradients();
        gradients.print_info(out);

        out << "Hessians basis functions:" << endl;
        auto hessians = elem->get_basis_hessians();
        hessians.print_info(out);
    }
}


int main()
{
    test<1, 1>();
    test<1, 2>();
    test<1, 3>();
    test<2, 1>();
    test<2, 2>();
    test<2, 3>();
    test<3, 1>();
    test<3, 3>();

    return 0;
}
