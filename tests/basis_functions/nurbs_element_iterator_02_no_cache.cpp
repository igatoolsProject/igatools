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
 *  Test for the NURBS space iterator using the no cache evaluations
 *
 *  author: pauletti
 *  date: Jun 11, 2014
 *
 */

#include "../tests.h"

#include <igatools/base/quadrature_lib.h>
#include <igatools/basis_functions/nurbs_space.h>
#include <igatools/basis_functions/nurbs_element.h>

using namespace EpetraTools;

template<int dim, int range, int rank = 1>
void test()
{
    OUTSTART
    const int r = 2;

    using Space = NURBSSpace< dim, range, rank >;
    auto  knots = CartesianGrid<dim>::create(3);

    auto degree = TensorIndex<dim>(r);
    auto bsp_space = BSplineSpace<dim,range,rank>::create(degree, knots);

    using ScalarSpSpace = BSplineSpace<dim,1,1>;
    auto scalar_bsp_space = ScalarSpSpace::create(degree, knots);

    const auto n_scalar_basis = scalar_bsp_space->get_num_basis_table()[0];

    using WeightFunc = IgFunction<dim,0,1,1>;
    SafeSTLVector<Real> weights_coef(n_scalar_basis.flat_size());
    const int n_entries = weights_coef.size();
    for (int i = 0 ; i < n_entries ; ++i)
        weights_coef[i] = (i+1) * (1.0 / n_entries) ;
    Epetra_SerialComm comm;
    auto map = create_map(*scalar_bsp_space, "active", comm);

    auto w_func = WeightFunc::create(scalar_bsp_space,
                                     std::make_shared<typename EpetraTools::Vector>(Copy, *map, weights_coef.data()));

    auto nrb_space = Space::create(bsp_space,w_func);

    const int n_points = 3;
    QGauss<dim> quad(n_points);

    auto elem     = nrb_space->begin();
    auto end_element = nrb_space->end();

    for (; elem != end_element; ++elem)
    {
        out << "Element: " << elem->get_flat_index()<< endl;

        out << "Values basis functions:" << endl ;
        const auto values = elem->template evaluate_basis_at_points<_Value>(quad,DofProperties::active);
        values.print_info(out) ;

        out << "Gradients basis functions:" << endl ;
        const auto gradients = elem->template evaluate_basis_at_points<_Gradient>(quad,DofProperties::active);
        gradients.print_info(out) ;

        out << "Hessians basis functions:" << endl ;
        const auto hessians = elem->template evaluate_basis_at_points<_Hessian>(quad,DofProperties::active);
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
