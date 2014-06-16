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
 *  date:
 */

#include "../tests.h"
#include <igatools/geometry/unit_element.h>
#include <igatools/basis_functions/spline_space.h>


template< int dim, int range, int rank >
void run_test()
{
	using SplineSpace = SplineSpace<dim, range, rank>;
	auto grid = CartesianGrid<dim>::create({3,4});
	typename SplineSpace::DegreeTable deg{{1,3}};
	SplineSpace space(deg, grid, SplineSpace::InteriorReg::maximum);

	for (auto  f : UnitElement<dim>::faces)
	{
		out << "face: " << f << endl;

		out << "Multiplicity:\n";
		auto f_mult = space.get_face_mult(f);
		for (auto x: *f_mult)
			x.print_info(out);

		out << "Degree:\n";
		auto f_deg = space.get_face_degree(f);
		for (auto x: f_deg)
			out << x;
		out << endl;
	}
}



int main()
{
    out.depth_console(10);

 //   run_test<1,1,1>();
    run_test<2,1,1>();
 //   run_test<3,1,1>();

    return  0;
}
