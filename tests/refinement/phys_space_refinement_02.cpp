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
 *  Test refinement of a basic PhysicalSpace using the NURBSSpace as reference space
 *  and the IdentityFunction as mapping.
 *
 *  author: pauletti
 *  date: 2013-10-02
 *
 */

#include "../tests.h"

#include <igatools/base/quadrature_lib.h>
#include <igatools/geometry/push_forward.h>
#include <igatools/basis_functions/nurbs_space.h>
#include <igatools/basis_functions/physical_space.h>
#include <igatools/basis_functions/physical_space_element.h>
#include <igatools/functions/identity_function.h>


template <int dim>
void test_evaluate()
{
    auto grid = CartesianGrid<dim>::create();
    grid->refine();
    out << endl;


    const int deg = 2;
    auto bsp_space = BSplineSpace<dim>::create(deg, grid);


    using ScalarSpSpace = BSplineSpace<dim,1,1>;
    auto scalar_bsp_space = ScalarSpSpace::create(deg, grid);

    const auto n_scalar_basis = scalar_bsp_space->get_num_basis();

    IgCoefficients weights_coef;
    for (int i = 0 ; i < n_scalar_basis ;)
    {
        weights_coef[i++] = 1.0;
        weights_coef[i++] = 0.4;
        weights_coef[i++] = 0.65;
        weights_coef[i++] = 1.0;
    }

    using WeightFunc = IgFunction<dim,0,1,1>;
    auto w_func = WeightFunc::create(scalar_bsp_space,weights_coef);

    using RefSpace = ReferenceSpace<dim>;
    using RefSpacePtr = std::shared_ptr<RefSpace>;
    RefSpacePtr ref_space = NURBSSpace<dim>::create(bsp_space,w_func);

    auto phys_space =
        PhysicalSpace<dim,1,1,0,Transformation::h_grad>::create(
            ref_space,IdentityFunction<dim>::create(grid));


    out << endl;
    out << endl;

    out << "===============================================================" << endl;
    out << "O R I G I N A L     S P A C E" << endl;
    phys_space->print_info(out);
    out << "===============================================================" << endl;
    out << endl;

    out << "===============================================================" << endl;
    out << "R E F I N E D     S P A C E" << endl;
    phys_space->refine_h();
    phys_space->print_info(out);
    out << "===============================================================" << endl;
    out << endl;
    //*/
}

int main()
{
    out.depth_console(10);

    test_evaluate<1>();
    out << endl ;

    test_evaluate<2>();
    out << endl ;

    test_evaluate<3>();
    out << endl ;

    return 0;
}
