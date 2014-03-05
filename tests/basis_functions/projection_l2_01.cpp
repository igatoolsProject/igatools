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
 *  Test for the l2_projection function.
 *  Bspline spaces case
 *
 *  author: pauletti
 *  date: 2013-10-10
 *  QA: The 0 dim case not returning appropriate value
 */

#include "../tests.h"

#include <igatools/base/quadrature_lib.h>
#include <igatools/basis_functions/bspline_space.h>
#include <igatools/basis_functions/space_tools.h>
#include <igatools/io/writer.h>


using numbers::PI;

template<int dim>
class BoundaryFunction : public Function<dim,1,1>
{
public:
    BoundaryFunction() : Function<dim,1,1>() {}

    iga::Real value(Point<dim> x) const
    {
        iga::Real f = 1;
        for (int i = 0; i<dim; ++i)
            f = f * cos(2*PI*x[i]);
        return f;
    }

    void evaluate(const std::vector< Point<dim> > &points, std::vector<Point<1> > &values) const
    {
        for (int i = 0; i<points.size(); ++i)
        {
            Point<dim> p = points[i];
            values[i][0] = this->value(p);
        }
    }

};


template<int dim , int range ,int rank>
void do_test(const int p)
{
    typedef BSplineSpace<dim,range,rank> space_ref_t;

    const int num_knots = 10;
    auto knots = CartesianGrid<dim>::create(num_knots);
    auto space = space_ref_t::create(knots, p) ;

    const int n_qpoints = 4;
    QGauss<dim> quad(n_qpoints);

    BoundaryFunction<dim> f;

    auto proj_values = space_tools::projection_l2(f,
                                                  const_pointer_cast<const space_ref_t>(space),quad);

    proj_values.print(out);

//    Writer<dim> output(knots, 4);
//    output.add_field(space, proj_values, "projected function");
//    string filename = "proj_function-" + to_string(dim) +"d";
//    output.save(filename);
}



int main()
{
    out.depth_console(20);
    do_test<0,1,1>(1);
    do_test<1,1,1>(3);
    do_test<2,1,1>(3);
    do_test<3,1,1>(1);

    return 0;
}

