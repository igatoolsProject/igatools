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
//TODO: This test should be splitted in: SafeSTLVector test, matrix test and solver test

/*
 * Test for developing epetra minimal and efficient linear algebra interaction
 *  author: pauletti
 *  date: 2015-03-30
 *
 */

#include "../tests.h"

#include <igatools/linear_algebra/epetra.h>
#include <igatools/base/quadrature_lib.h>
#include <igatools/basis_functions/bspline.h>
#include <igatools/basis_functions/bspline_element.h>

using namespace Teuchos;

template<int dim>
void matrix_map(const int deg, const int n_knots)
{
  OUTSTART
  auto grid = Grid<dim>::create(n_knots);
  auto space = SplineSpace<dim>::create(deg, grid);
  auto basis = BSpline<dim>::create(space);

  Epetra_SerialComm comm;

//    auto map = EpetraTools::create_map(*basis, "active", comm);
  auto graph = EpetraTools::create_graph(*basis, "active", *basis, "active",comm);

  auto matrix = EpetraTools::create_matrix(*graph);
  auto rhs = EpetraTools::create_vector(matrix->RangeMap());
  auto sol = EpetraTools::create_vector(matrix->DomainMap());

  using Flags = basis_element::Flags;
  auto flag = Flags::value | Flags::w_measure;
  auto elem_handler = basis->create_cache_handler();
  elem_handler->set_element_flags(flag);

  auto elem = basis->begin();
  auto end  = basis->end();

  auto quad = QGauss<dim>::create(2);
  elem_handler->init_element_cache(elem,quad);

  const int n_qp = quad->get_num_points();

  for (; elem != end; ++elem)
  {
    const int n_basis = elem->get_num_basis();

    DenseMatrix loc_mat(n_basis, n_basis);
    loc_mat = 0.0;
    DenseVector loc_rhs(n_basis);
    loc_rhs = 0.0;

    elem_handler->fill_element_cache(elem);
    auto phi = elem->get_element_values();
    auto w_meas = elem->get_element_w_measures();

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
    rhs->add_block(loc_dofs, loc_rhs);
  }

  matrix->FillComplete();

  auto solver = EpetraTools::create_solver(*matrix, *sol, *rhs);
  auto result = solver->solve();
  AssertThrow(result == Belos::ReturnType::Converged,
              ExcMessage("No convergence."));
  out << solver->getNumIters() << endl;

  matrix->print_info(out);
  rhs->print_info(out);
  sol->print_info(out);

  OUTEND
}



int main()
{
  const int deg = 1;
  const int n_knots = 10;
  matrix_map<2>(deg, n_knots);
  return 0;

}
