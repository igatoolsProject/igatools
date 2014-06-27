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
 *  Test for the integrate_difference function.
 *
 *  author: pauletti
 *  date: 26 Jun 2014
 */

#include "../tests.h"

#include <igatools/base/quadrature_lib.h>
#include <igatools/basis_functions/bspline_space.h>
#include <igatools/basis_functions/space_tools.h>
#include <igatools/io/writer.h>

/**
 * The p-norm of this function on the unit square is
 * (1/(p+1))^(n/p)
 */
template<int dim>
class NormFunction : public Function<dim>
{
public:
    using Base = Function<dim>;
    using typename Base::Point;
    using typename Base::Value;
    using typename Base::Gradient;
    using typename Base::Hessian;

    void evaluate(const std::vector<Point> &points,
                  std::vector<Value> &values) const
    {
        auto pt = points.begin();
        for (auto &val : values)
        {
            val = pt->norm();
            ++pt;
        }
    }

    void evaluate_gradients(
            const std::vector<Point> &points,
            std::vector<Gradient> &gradient) const
    {}

    void evaluate_hessians(
            const std::vector<Point> &points,
            std::vector<Hessian> &hessians) const
    {}
};



template<int dim, int range = 1, int rank = 1>
void test(const int deg,  const int n_knots)
{
    using Space = BSplineSpace<dim, range, rank>;


    auto knots = CartesianGrid<dim>::create(n_knots);
    auto space = Space::create(deg, knots);

    const int n_qpoints = ceil((2*deg + 1)/2.);
    QGauss<dim> quad(n_qpoints);

    NormFunction<dim> f;

#if defined(USE_TRILINOS)
    const auto la_pack = LAPack::trilinos;
#elif defined(USE_PETSC)
    const auto la_pack = LAPack::petsc;
#endif

    auto coeffs = space_tools::projection_l2<Space,la_pack>(f,space, quad);

    vector<Real> elem_err(space->get_grid()->get_num_elements());

    Real err = space_tools::integrate_difference<Space,la_pack>
    (f, space, quad, Norm::L2, coeffs, elem_err);

    out << err << endl;
    //out << elem_err << endl;;

//    Writer<dim> output(knots, 4);
//    output.add_field(space, proj_values, "projected function");
//    string filename = "proj_function-" + to_string(dim) +"d";
//    output.save(filename);
}



int main()
{
    out.depth_console(20);

    for (int n=2; n< std::pow(2,5); n*=2)
    {
        //    test<1,1,1>(1);
        //    test<2,1,1>(1);
        test<3,1,1>(1,n);
    }

    return 0;
}

