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

#include <igatools/base/quadrature_lib.h>
#include <igatools/basis_functions/nurbs_space.h>
#include <igatools/basis_functions/nurbs_element_accessor.h>

/*
 *  Test for the NURBS space iterator using the no cache evaluations
 *
 *  author: pauletti
 *  date: Jun 11, 2014
 *
 */

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


    //----------------------------------------------------------------------------------------------
    // for the basis functions evaluation we need a set of points (with tensor product structure)
    // to do so, we get the points from a Gauss quadrature scheme with 3 points

    const int n_points = 3;
    QGauss<dim> quad(n_points);
    const auto points = quad.get_points().get_flat_cartesian_product();


    auto elem     = space->begin();
    auto end = space->end();

    for (; elem != end; ++elem)
    {
        out << "Element: " << elem->get_flat_index()<< endl;

        out << "Values basis functions:" << endl ;
        const auto values = elem->evaluate_basis_values_at_points(points);
        values.print_info(out) ;

        out << "Gradients basis functions:" << endl ;
        const auto gradients = elem->evaluate_basis_gradients_at_points(points);
        gradients.print_info(out) ;

        out << "Hessians basis functions:" << endl ;
        const auto hessians = elem->evaluate_basis_hessians_at_points(points);
        hessians.print_info(out) ;
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
