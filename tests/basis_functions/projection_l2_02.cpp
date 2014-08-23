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
 *  Physical spaces version
 *  author: pauletti
 *  date: 2013-10-10
 *
 */

#include "../tests.h"

#include <igatools/base/quadrature_lib.h>
#include <igatools/basis_functions/bspline_space.h>
#include <igatools/basis_functions/physical_space.h>
#include <igatools/geometry/mapping_lib.h>
#include <igatools/geometry/push_forward.h>
#include <igatools/basis_functions/space_tools.h>
#include <igatools/io/writer.h>

using numbers::PI;

template<int dim>
class BoundaryFunction : public Function<dim,1,1>
{
public:
    BoundaryFunction() : Function<dim,1,1>() {}

    Real value(Points<dim> x) const
    {
        Real f = 1;
        for (int i = 0; i<dim; ++i)
            f = f * cos(2*PI*x[i]);
        return f;
    }

    void evaluate(const ValueVector< Points<dim> > &points, ValueVector<Points<1> > &values) const
    {
        for (int i = 0; i<points.size(); ++i)
        {
            Points<dim> p = points[i];
            values[i][0] = this->value(p);
        }
    }

};


template<int dim , int spacedim>
void do_test(const int p)
{
    const int codim = spacedim-dim;
    using space_ref_t = BSplineSpace<dim,1,1>;
    using pushforward_t = PushForward<Transformation::h_grad, dim, codim>;
    using space_t = PhysicalSpace<space_ref_t, pushforward_t>;

    const int num_knots = 10;
    auto grid = CartesianGrid<dim>::create(num_knots);
    auto ref_space = space_ref_t::create(p, grid);
    Points<spacedim> b;
    Derivatives<dim, spacedim, 1, 1> A;
    for (int i = 0; i < dim; ++i)
    {
        A[i][i] = 1+i;
    }
    auto map = LinearMapping<dim, codim>::create(grid, A, b);
    auto pf = pushforward_t::create(map);
    auto space = space_t::create(ref_space,pf);


    const int n_qpoints = 4;
    QGauss<dim> quad(n_qpoints);

    BoundaryFunction<dim> f;

#if defined(USE_TRILINOS)
    const auto la_pack = LAPack::trilinos;
#elif defined(USE_PETSC)
    const auto la_pack = LAPack::petsc;
#endif

    auto proj_values = space_tools::projection_l2
                       <space_t,la_pack>(
                           f,const_pointer_cast<const space_t>(space),quad);

    proj_values.print(out);

    Writer<dim> output(grid, 4);
    output.add_field(space, proj_values, "projected function");
    string filename = "proj_function-" + to_string(dim) +"d";
    output.save(filename);
}



int main()
{
    out.depth_console(20);
    // do_test<1,1,1>(3);
    do_test<2,2>(3);
    //do_test<3,1,1>(1);

    return 0;
}

