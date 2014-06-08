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
 *  Test for the BsplineSpace class face subspace extraction function.
 *  Here we print the information of the face spaces thus extracted.
 *  author: pauletti
 *  date: Jan 29, 2013
 *  updated: 2013-04-02
 */

#include "../tests.h"
#include <igatools/basis_functions/bspline_space.h>
#include <igatools/basis_functions/space_tools.h>

template< int dim_domain, int dim_range, int rank >
void run_test()
{
    out << "-------------------------" <<endl;
    out << "dim_dimain = " << dim_domain << endl;
    out << "dim_range  = " << dim_range << endl;
    out << "rank       = " << rank << endl;

    const int degree=1;

    auto grid = CartesianGrid<dim_domain>::create(3);
    auto space = BSplineSpace<dim_domain, dim_range, rank>::create(degree, grid);

    vector<Index> dof_map;

    int face_id = 0;

    auto face_space = space->get_face_space(face_id,dof_map);

    for (vector<Index>::iterator it=dof_map.begin() ; it < dof_map.end(); ++it)
        out<< "face_id = "<< face_id << ", dof_id = "<< *it << endl;

    out << "-------------------------" <<endl;
}



int main()
{
    //out.depth_console(10);

    run_test<2,1,1>();
    run_test<2,2,1>();
    run_test<3,1,1>();

    return  0;
}
