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
 *  Test for building a matrix on a space of an igfunction
 *
 *  author: pauletti
 *  date: 2015-03-17
 *
 */

#include "../tests.h"

#include <igatools/basis_functions/bspline_space.h>
#include <igatools/basis_functions/bspline_element.h>
#include <igatools/linear_algebra/distributed_matrix.h>



template<int dim, int range = 1, int rank = 1>
void filtered_dofs(const int deg = 1, const int n_knots = 3)
{
    OUTSTART

    std::string interior = "interior";
    using RefSpace = ReferenceSpace<dim, range, rank>;
    using Space = BSplineSpace<dim, range, rank>;

    auto grid = CartesianGrid<dim>::create(n_knots);
    auto space = Space::create(deg, grid);
    auto dof_dist = space->get_dof_distribution();
    dof_dist->add_dofs_property(interior);
    std::set<Index> int_dofs={4};
    dof_dist->set_dof_property_status(interior, int_dofs,true);

    auto elem = space->begin();
    auto end  = space->end();

    for(;elem != end; ++elem)
    {
        auto loc_to_glob = elem->get_local_to_global("interior");
        loc_to_glob.print_info(out);
        out << endl;
    }

    OUTEND
}




int main()
{
	const int dim = 2;
	filtered_dofs<dim>();

    return 0;
}
