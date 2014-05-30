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
 *  Test for the boundary projection function.
 *  This test ....
 *  author: pauletti
 *  date: 2013-03-19
 *
 */

#include "../tests.h"

#include <igatools/base/quadrature_lib.h>
#include <igatools/basis_functions/bspline_space.h>

#include <igatools/basis_functions/space_tools.h>
#include <igatools/linear_algebra/dof_tools.h>


template<int dim>
class BoundaryFunction : public Function<dim,1,1>
{

public:
    BoundaryFunction() : Function<dim,1,1>() {}

    void evaluate(const vector< Point<dim> > &points,
                  std::vector<Point<1> > &values) const
    {
        for (int i =0; i<points.size(); i++)
        {
            Point<dim> p = points[i];
            values[i][0] = p[0];
        }
    };

};


template<int dim>
void do_test(const int p)
{
    out << "Dimension: " << dim << endl;
    typedef BSplineSpace<dim> space_ref_t;

    const int num_knots = 2;
    TensorSize<dim> n_knots;
    for (int i = 0; i < dim; ++i)
    {
        n_knots[i] = num_knots+i;
    }
    auto grid = CartesianGrid<dim>::create(n_knots);
    auto space = space_ref_t::create(grid, p) ;

    const int n_qpoints = 4;
    QGauss<dim-1> quad(n_qpoints);

    BoundaryFunction<dim> f;

    const boundary_id dirichlet = 1;
    grid->set_boundary_id(2, dirichlet);
    std::set<boundary_id> face_id;
    face_id.insert(dirichlet);


#if defined(USE_TRILINOS)
    const auto linear_algebra_package = LinearAlgebraPackage::trilinos;
#elif defined(USE_PETSC)
    const auto linear_algebra_package = LinearAlgebraPackage::petsc;
#endif

    std::map<Index,iga::Real> boundary_values;
    space_tools::project_boundary_values<space_ref_t,linear_algebra_package>(
        f, const_pointer_cast<const space_ref_t>(space), quad, face_id,
        boundary_values);

    out << "basis index \t value" << endl;
    for (auto entry : boundary_values)
        out << entry.first << "\t" << entry.second << endl;

}



int main()
{
    out.depth_console(20);

    do_test<2>(2);
//    do_test<2,1,1>(3);
//    do_test<3,1,1>(2);

    return 0;
}

