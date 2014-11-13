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
#include <igatools/base/formula_function.h>

#include <igatools/basis_functions/new_bspline_space.h>
#include <igatools/basis_functions/bspline_element.h>

#include <igatools/basis_functions/space_tools.h>

#include <igatools/io/writer.h>


using numbers::PI;

template<int dim>
class BoundaryFunction : public FormulaFunction<dim>
{
private:
   using base_t = NewFunction<dim>;
    using parent_t = FormulaFunction<dim>;
   using self_t = BoundaryFunction<dim>;
    using typename base_t::GridType;
public:
    using typename parent_t::Point;
    using typename parent_t::Value;
    template <int order>
        using Derivative = typename parent_t::template Derivative<order>;
public:
    BoundaryFunction(std::shared_ptr<GridType> grid)
    : FormulaFunction<dim>(grid)
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

    Real value(Points<dim> x) const
    {
        Real f = 1;
        for (int i = 0; i<dim; ++i)
            f = f * cos(2*PI*x[i]);
        return f;
    }

    void evaluate_0(const ValueVector<Point> &points,
                       ValueVector<Value> &values) const override
    {
        for (int i = 0; i<points.size(); ++i)
        {
            Points<dim> p = points[i];
            values[i][0] = this->value(p);
        }
    }
    void evaluate_1(const ValueVector<Point> &points,
                       ValueVector<Derivative<1>> &values) const override
                               {}

       void evaluate_2(const ValueVector<Point> &points,
                       ValueVector<Derivative<2>> &values) const override
                               {}
};


template<int dim , int range ,int rank, LAPack la_pack>
void do_test(const int p, const int num_knots = 10)
{
    using Space =  NewBSplineSpace<dim,range,rank>;

    auto knots = CartesianGrid<dim>::create(num_knots);
    auto space = Space::create(p, knots) ;

    const int n_qpoints = 4;
    QGauss<dim> quad(n_qpoints);

    auto f = BoundaryFunction<dim>::create(knots);




    auto proj_func = space_tools::projection_l2<Space,la_pack>(f, space, quad);
    proj_func->print_info(out);

//    Writer<dim> output(knots, 4);
//    output.add_field(space, proj_values, "projected function");
//    string filename = "proj_function-" + to_string(dim) +"d";
//    output.save(filename);
}



int main()
{
#if defined(USE_TRILINOS)
    const auto la_pack = LAPack::trilinos;
#elif defined(USE_PETSC)
    const auto la_pack = LAPack::petsc;
#endif
    out.depth_console(20);
//    // do_test<0,1,1>(1);
    do_test<1,1,1, la_pack>(3);
    do_test<2,1,1, la_pack>(3);
    do_test<3,1,1, la_pack>(1);

    return 0;
}

