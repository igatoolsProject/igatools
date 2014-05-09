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

    iga::Real value(Point<dim> P_) const
    {
        iga::Real  PI = numbers::PI;

        iga::Real f = 1;
        for (int cnt = 0; cnt<dim; cnt++)
        {
            f = f * cos(iga::Real(2*PI*P_[cnt]));
        }

        return f;
    }

    void evaluate(const std::vector< Point<dim> > &points, std::vector<Point<1> > &values) const
    {
        for (int i =0; i<points.size(); i++)
        {
            Point<dim> p = points[i];
            values[i][0] = this->value(p);
        }
    };

};


template<int dim_ref_domain ,int dim_phys_domain,int dim_range ,int rank>
void do_test(const int p)
{
    const int codim = dim_phys_domain - dim_ref_domain;
    typedef BSplineSpace<dim_ref_domain,dim_range,rank> space_ref_t ;
    typedef PushForward<Transformation::h_grad,dim_ref_domain,codim> PushForward ;
    typedef PhysicalSpace<space_ref_t,PushForward> space_phys_t ;

    const int num_knots = 10;
    auto knots = CartesianGrid<dim_ref_domain>::create(num_knots);
    auto space = space_ref_t::create(knots, p) ;
    auto map = IdentityMapping<dim_ref_domain,codim>::create(knots);
    auto phys_space = space_phys_t::create(space, PushForward::create(map));

    //Quadrature
    const int n_qpoints = 4;
    QGauss<dim_ref_domain-1> quad(n_qpoints);

    BoundaryFunction<dim_phys_domain> bc;

    knots->set_boundary_id(0,1);
    std::set<boundary_id> face_id;
    face_id.insert(1);

#if defined(USE_TRILINOS)
    const auto linear_algebra_package = LinearAlgebraPackage::trilinos;
#elif defined(USE_PETSC)
    const auto linear_algebra_package = LinearAlgebraPackage::petsc;
#endif

    std::map<Index,iga::Real> boundary_values;
    space_tools::project_boundary_values<space_phys_t,linear_algebra_package>(
        bc,
        phys_space,
        quad,
        face_id,
        boundary_values);

    for (auto entry: boundary_values)
        out << entry.first << "\t" << entry.second << endl;
}



int main()
{
    out.depth_console(20);
    do_test<2,2,1,1>(3);

    //do_test<3,3,1,1>(3);

    //do_test<2,3,1,0>(3);

    return 0;
}

