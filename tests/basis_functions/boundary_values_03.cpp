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
 *  Test for the boundary projection function.
 *  This test ....
 *  author: pauletti
 *  date: 2013-03-19
 *
 */
//TODO: this test should be merge into the other ones

#include "../tests.h"

#include <igatools/base/quadrature_lib.h>
#include <igatools/base/identity_function.h>
#include <igatools/base/formula_function.h>
#include <igatools/basis_functions/bspline_space.h>
#include <igatools/basis_functions/bspline_element.h>

#include <igatools/basis_functions/space_tools.h>
#include <igatools/linear_algebra/dof_tools.h>


template<int dim>
class XProject : public FormulaFunction<dim>
{
private:
    using base_t = Function<dim>;
    using parent_t = FormulaFunction<dim>;
    using self_t = XProject<dim>;
    using typename base_t::GridType;
public:
    using typename parent_t::Point;
    using typename parent_t::Value;
    template <int order>
    using Derivative = typename parent_t::template Derivative<order>;
public:
    XProject(std::shared_ptr<GridType> grid)
        : FormulaFunction<dim>(grid, IdentityFunction<dim>::create(grid))
    {}

    static std::shared_ptr<base_t>
    create(std::shared_ptr<GridType> grid)
    {
        return std::shared_ptr<base_t>(new self_t(grid));
    }

    std::shared_ptr<base_t> clone() const override
    {
        return std::make_shared<self_t>(self_t(*this));
    }

    void evaluate_0(const ValueVector<Point> &points,
                    ValueVector<Value> &values) const override
    {
        for (int i = 0; i<points.size(); ++i)
        {
            Points<dim> p = points[i];
            values[i][0] = p[0];
        }
    }
    void evaluate_1(const ValueVector<Point> &points,
                    ValueVector<Derivative<1>> &values) const override
    {}

    void evaluate_2(const ValueVector<Point> &points,
                    ValueVector<Derivative<2>> &values) const override
    {}
};




template<int dim , int range ,int rank>
void do_test(const int p, TensorSize<dim> n_knots)
{
    const int sub_dim = dim - 1;
    out << "Dimension: " << dim << endl;
    using Space = BSplineSpace<dim, range, rank>;


    auto grid = CartesianGrid<dim>::create(n_knots);
    auto space = Space::create(p, grid) ;
    auto f = XProject<dim>::create(grid);

    const int n_qpoints = 4;
    QGauss<sub_dim> quad(n_qpoints);

    const boundary_id dirichlet = 1;
    grid->set_boundary_id(2, dirichlet);
    std::set<boundary_id> bdry_ids;
    bdry_ids.insert(dirichlet);

    std::map<Index,Real> boundary_values;
    space_tools::project_boundary_values<Space>(
        f, space, quad, bdry_ids,
        boundary_values);

    out << "basis index \t value" << endl;
    for (auto entry : boundary_values)
        out << entry.first << "\t" << entry.second << endl;

}


int main()
{
    {
        const int dim = 2;
        TensorSize<dim> n_knots { arr::sequence<dim>(2)};
        do_test<dim, 1, 1>(2, n_knots);
    }
//    do_test<2,1,1>(3);
//    do_test<3,1,1>(2);

    return 0;
}

