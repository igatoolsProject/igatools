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
#include <igatools/basis_functions/new_bspline_space.h>
#include <igatools/base/formula_function.h>
#include <igatools/base/identity_function.h>
#include <igatools/base/function_lib.h>
#include <igatools/basis_functions/space_tools.h>
//#include <igatools/io/writer.h>

/**
 * The p-norm of this function on the unit square is
 * (1/(p+1))^(n/p)
 */
template<int dim, int codim=0, int range = 1, int rank = 1>
class ProductFunction : public FormulaFunction<dim>
{
private:
    using base_t = NewFunction<dim, codim, range, rank>;
    using parent_t = FormulaFunction<dim, codim, range, rank>;
    using self_t = ProductFunction<dim, codim, range, rank>;
public:
    using typename parent_t::Point;
    using typename parent_t::Value;
    using typename parent_t::Gradient;
    using typename parent_t::ElementIterator;
    using typename parent_t::ElementAccessor;
    template <int order>
    using Derivative = typename parent_t::template Derivative<order>;

    using parent_t::FormulaFunction;
//    static std::shared_ptr<base_t>
//    create(std::shared_ptr<const CartesianGrid<dim>> grid,
//           const Value &b);

//protected:
//    ConstantFunction(std::shared_ptr<const CartesianGrid<dim>> grid,
//                     const Value &b);

private:
    void evaluate_0(const ValueVector<Point> &points,
                    ValueVector<Value> &values) const
    {
        auto pt = points.begin();
        auto val = values.begin();

        for (; pt != points.end(); ++pt, ++val)
        {
            *val = 1.;
            for (int i=0; i<dim; ++i)
                (*val) *= (*pt)[i];
        }
    }

    void evaluate_1(const ValueVector<Point> &points,
                    ValueVector<Derivative<1>> &values) const
    {}

    void evaluate_2(const ValueVector<Point> &points,
                    ValueVector<Derivative<2>> &values) const
    {}

};



template<int dim, int range = 1, int rank = 1>
void do_test(const int deg)
{
    using Space = NewBSplineSpace<dim, range, rank>;

    const int n_knots = 10;
    auto grid = CartesianGrid<dim>::create(n_knots);
    auto space = Space::create(deg, grid);

    const int n_qpoints = ceil((2*dim + 1)/2.);
    QGauss<dim> quad(n_qpoints);

    ProductFunction<dim> f(grid, IdentityFunction<dim>::create(grid));
    typename functions::ConstantFunction<dim,0,1>::Value val {0.};
    auto g = functions::ConstantFunction<dim,0,1>::create(grid, IdentityFunction<dim>::create(grid), val);

    vector<Real> elem_err(grid->get_num_active_elems());
    Real err = space_tools::integrate_difference<dim>
               (f, *g, quad, Norm::L2, elem_err);

    const Real p=2;
    out << std::pow(p+1, -dim/p) << "\t" << err << endl;

#if 0
#if defined(USE_TRILINOS)
    const auto la_pack = LAPack::trilinos;
#elif defined(USE_PETSC)
    const auto la_pack = LAPack::petsc;
#endif



    Vector<la_pack> coeffs(space->get_num_basis());
    vector<Real> elem_err(space->get_grid()->get_num_active_elems());

    Real err = space_tools::integrate_difference<Space,la_pack>
               (f, space, quad, Norm::L2, coeffs, elem_err);

    out << std::pow(p+1, -dim/p) << "\t" << err << endl;
    // out << elem_err << endl;;

//    Writer<dim> output(knots, 4);
//    output.add_field(space, proj_values, "projected function");
//    string filename = "proj_function-" + to_string(dim) +"d";
//    output.save(filename);
#endif
}



int main()
{
    out.depth_console(20);
    // do_test<0,1,1>(1);
    do_test<1,1,1>(3);
    do_test<2,1,1>(3);
    do_test<3,1,1>(1);

    return 0;
}

