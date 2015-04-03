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
//TODO: This test should be splitted in: vector test, matrix test and solver test

/*
 * Test for developing epetra minimal and efficient linear algebra interaction
 *  author: pauletti
 *  date: 2015-03-30
 *
 */

#include "../tests.h"


#include <igatools/linear_algebra/epetra.h>

#include <igatools/base/quadrature_lib.h>
#include <igatools/basis_functions/bspline_space.h>
#include <igatools/basis_functions/bspline_element.h>





using namespace Teuchos;


template<int dim>
void matrix_map(const int deg, const int n_knots)
{
	OUTSTART
	auto grid = CartesianGrid<dim>::create(n_knots);
	auto space = BSplineSpace<dim>::create(deg, grid);

	Epetra_SerialComm comm;

	auto map = EpetraTools::create_map(space, "active", comm);
	auto graph = EpetraTools::create_graph(space, "active", space, "active",map, map);

	auto matrix = EpetraTools::create_matrix(graph);
	auto vector = EpetraTools::create_vector(map);
	auto sol = EpetraTools::create_vector(map);
//	matrix->Print(out.get_file_stream());
//	vector->Print(out.get_file_stream());

	auto quad = QGauss<dim>(2);
	auto flag = ValueFlags::value | ValueFlags::w_measure;
	auto elem_handler = space->create_elem_handler();
	elem_handler->reset(flag, quad);

	auto elem = space->begin();
	auto end  = space->end();
	elem_handler->init_element_cache(elem);
	const int n_qp = quad.get_num_points();

	for (; elem != end; ++elem)
	{
		const int n_basis = elem->get_num_basis("active");

		DenseMatrix loc_mat(n_basis, n_basis);
		loc_mat = 0.0;
		DenseVector loc_rhs(n_basis);
		loc_rhs = 0.0;

		elem_handler->fill_element_cache(elem);
		auto phi = elem->template get_values<0, dim>(0, "active");
		auto w_meas = elem->template get_w_measures<dim>(0);

		for (int i = 0; i < n_basis; ++i)
		{
			auto phi_i = phi.get_function_view(i);
			for (int j = 0; j < n_basis; ++j)
			{
				auto phi_j = phi.get_function_view(j);
				for (int qp = 0; qp < n_qp; ++qp)
					loc_mat(i,j) +=
							scalar_product(phi_i[qp], phi_j[qp])
							* w_meas[qp];
			}

			for (int qp=0; qp<n_qp; ++qp)
				loc_rhs(i) += phi_i[qp][0] // f=1
				* w_meas[qp];
		}

		const auto loc_dofs = elem->get_local_to_global("active");
		matrix->add_block(loc_dofs, loc_dofs, loc_mat);
		vector->add_block(loc_dofs, loc_rhs);
	}

	matrix->FillComplete();

	auto solver = EpetraTools::create_solver(matrix, sol, vector);
	auto result = solver->solve();
	AssertThrow(result == Belos::ReturnType::Converged,
			ExcMessage("No convergence."));
	out << solver->getNumIters() << endl;

	matrix->Print(out.get_file_stream());
	vector->Print(out.get_file_stream());
	sol->Print(out.get_file_stream());

	OUTEND
}



int main()
{
    const int deg = 1;
    const int n_knots = 10;
	matrix_map<2>(deg, n_knots);
	return 0;

}
