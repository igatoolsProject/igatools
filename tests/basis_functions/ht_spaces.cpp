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
#include "../tests.h"

#include <igatools/basis_functions/bspline_space.h>
#include <igatools/basis_functions/bspline_element_accessor.h>

#include <igatools/base/quadrature_lib.h>

template <int dim_ref, int dim_phys, int dim_range, int rank>
void do_test()
{
    typedef BSplineSpace< dim_ref, dim_range, rank > vel_space_ref_t;
    typedef BSplineSpace< dim_ref, 1, 1 > prex_space_ref_t;

    shared_ptr<vel_space_ref_t> vel_space;
    shared_ptr<prex_space_ref_t> prex_space;

    out << "defining space:" << endl;
    const int deg = 1, reg = 0, n_knots = 3;
    const int mult = deg-reg;

    std::array< std::vector< Index >, dim_ref > mult_vector;
    TensorIndex< dim_ref > deg_vector;

    for (int i=0; i<dim_ref; i++)
    {
        mult_vector[i].push_back(deg+1);
        for (int j=1; j<n_knots-1; j++)
            mult_vector[i].push_back(mult);
        mult_vector[i].push_back(deg+1);
        deg_vector[i]=deg;
    }

    auto mesh = CartesianGrid<dim_ref>::create(n_knots);

    prex_space.reset(new prex_space_ref_t(
                         mesh,
                         typename prex_space_ref_t::MultiplicityTable(Multiplicity<dim_ref>(mult_vector)),
                         StaticMultiArray< TensorIndex<dim_ref>,1,1>(deg_vector))) ;

    out << "prex #dofs = " << prex_space->get_num_basis() << endl;

    for (int i = 0; i < mult_vector.size(); i++)
        for (int j = 0; j < mult_vector[i].size(); j++)
            mult_vector[i][j] += 1;
    deg_vector.fill(deg+1);

    vel_space.reset(new vel_space_ref_t(
                        mesh,
                        StaticMultiArray< Multiplicity<dim_ref>,dim_range,rank>(Multiplicity<dim_ref>(mult_vector)),
                        StaticMultiArray< TensorIndex<dim_ref>,dim_range,rank>(deg_vector))) ;

    out << "vel  #dofs = " << vel_space->get_num_basis() << endl;

    int n_basis_vel = vel_space->get_num_basis_per_element();

    out << "vel  #dofs per element = " << n_basis_vel << endl;
    vector<Index> vel_local_dofs(n_basis_vel);

    int n_basis_prex = prex_space->get_num_basis_per_element();
    vector<Index> prex_local_dofs(n_basis_prex);

    out << "prex #dofs per element = " << n_basis_prex << endl;

    auto vel_el = vel_space->begin();
    auto prex_el = prex_space->begin();
    const auto end_e = vel_space->end();

    QGauss< dim_ref > quad(4);
    vel_el->init_values(ValueFlags::value,quad);
    prex_el->init_values(ValueFlags::value,quad);

    for (; vel_el != end_e; ++vel_el, ++prex_el)
    {
        vel_el->fill_values();
        prex_el->fill_values();

        prex_local_dofs = prex_el->get_local_to_global();
        vel_local_dofs = vel_el->get_local_to_global();

        out << "element = " << vel_el->get_flat_index() <<endl;
        for (int i = 0; i < n_basis_vel; i++)
        {
            out << "row = " << vel_local_dofs[i] << ", cln = " ;
            for (int j = 0; j <n_basis_prex; j ++)
                out  << prex_local_dofs[j]<< ", ";
            out  << endl;
        }
    }
}

int main(int argc, char *argv[])
{
    out.depth_console(10);
    do_test< 2, 2, 2, 1 >() ;

    return 0;
}
