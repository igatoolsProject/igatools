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

#include "../tests.h"

#include <igatools/base/quadrature_lib.h>
#include <igatools/geometry/identity_mapping.h>

#include <igatools/basis_functions/bspline_space.h>
#include <igatools/basis_functions/physical_space.h>

#include <igatools/basis_functions/space_tools.h>
#include <igatools/linear_algebra/dof_tools.h>


template<int dim>
class BoundaryFunction : public Function<dim,1,1>
{

public:
    BoundaryFunction() : Function<dim,1,1>() {}

    Real value(Points<dim> P_) const
    {
        Real  PI = numbers::PI;

        Real f = 1;
        for (int cnt = 0; cnt<dim; cnt++)
        {
            f = f * cos(Real(2*PI*P_[cnt]));
        }

        return f;
    }

    void evaluate(const ValueVector< Points<dim> > &points, ValueVector<Points<1> > &values) const
    {
        for (int i =0; i<points.size(); i++)
        {
            Points<dim> p = points[i];
            values[i][0] = this->value(p);
        }
    };

};


template<int dim, int space_dim, int range, int rank=1>
void do_test(const int p)
{
    const int codim = space_dim - dim;
    using RefSpace = BSplineSpace<dim,range,rank>;

    typedef PushForward<Transformation::h_grad,dim,codim> PushForward ;
    typedef PhysicalSpace<RefSpace, PushForward> space_phys_t ;

    const int num_knots = 10;
    auto knots = CartesianGrid<dim>::create(num_knots);
    auto space = RefSpace::create(p, knots) ;
    auto map = IdentityMapping<dim,codim>::create(knots);
    auto phys_space = space_phys_t::create(space, PushForward::create(map));

    //Quadrature
    const int n_qpoints = 4;
    QGauss<dim-1> quad(n_qpoints);

    BoundaryFunction<space_dim> bc;

    knots->set_boundary_id(0,1);
    std::set<boundary_id> face_id;
    face_id.insert(1);

#if defined(USE_TRILINOS)
    const auto la_pack = LAPack::trilinos;
#elif defined(USE_PETSC)
    const auto la_pack = LAPack::petsc;
#endif

    std::map<Index, Real> boundary_values;
    space_tools::project_boundary_values<space_phys_t,la_pack>(
        bc, phys_space, quad, face_id, boundary_values);

    for (auto entry: boundary_values)
        out << entry.first << "\t" << entry.second << endl;
}



int main()
{
    out.depth_console(20);
    do_test<2,2,1>(3);

    //do_test<3,3,1,1>(3);

    //do_test<2,3,1,0>(3);

    return 0;
}

