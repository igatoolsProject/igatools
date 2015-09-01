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
 *  Test for periodic BSplineSpace solving poisson
 *
 *  author: pauletti
 *  date: 2015-03-12
 *
 */

// TODO (pauletti, Mar 12, 2015): this is not a unit test, this tests
// has to be splitted into simpler tests (this one goes onto consistency test)

#include "../tests.h"
#include <igatools/functions/function_lib.h>
#include <igatools/basis_functions/bspline_space.h>
#include <igatools/basis_functions/bspline_element.h>
#include <igatools/base/quadrature_lib.h>

#include <igatools/linear_algebra/epetra_solver.h>
#include <igatools/io/writer.h>
#include <igatools/functions/identity_function.h>
#include <igatools/functions/ig_function.h>

#include <igatools/basis_functions/space_tools.h>
#include <igatools/linear_algebra/dof_tools.h>

using space_tools::project_boundary_values;
using dof_tools::apply_boundary_values;
using functions::ConstantFunction;

using namespace EpetraTools;

template <int dim>
void assemble_matrix(const int n_knots, const int deg)
{
  using Space  = BSplineSpace<dim>;
  using RefSpace  = ReferenceSpace<dim>;

  using Function = Function<dim,0,1,1>;
  using ConstFunction = functions::LinearFunction<dim,0,1>;
  using Value = typename Function::Value;
  using Gradient = typename Function::Gradient;

  typename Space::Degrees degt(deg);
  typename Space::Periodicity periodic(false);
  periodic[0] = true;
  typename Space::EndBehaviour end_b(BasisEndBehaviour::interpolatory);

  end_b[0] = BasisEndBehaviour::periodic;

  auto grid  = CartesianGrid<dim>::create(n_knots);
  grid->set_boundary_id(0, 1);
  grid->set_boundary_id(1, 1);

  auto space = Space::create(degt, grid, InteriorReg::maximum, periodic, end_b);

  space->print_info(out);


  Gradient A;
  Value b = {-5.};
  for (int i = 0; i < dim; ++i)
  {
    A[i]=10*(i+1);
  }

  auto f = ConstFunction::create(grid, IdentityFunction<dim>::create(grid), A, b);

  auto matrix = create_matrix(*space,DofProperties::active,Epetra_SerialComm());
  auto rhs = create_vector(matrix->RangeMap());
  auto solution = create_vector(matrix->DomainMap());


  const QGauss<dim>  elem_quad(deg);
  auto elem_handler = space->create_cache_handler();
  auto flag = ValueFlags::value | ValueFlags::gradient |ValueFlags::w_measure;
  elem_handler->reset(flag, elem_quad);
  f->reset(ValueFlags::value, elem_quad);

  auto elem   = space->begin();
  const auto elem_end = space->end();
  auto f_elem = f->begin();

  elem_handler->init_element_cache(elem);
  f->init_element_cache(f_elem);

  const int n_qp = elem_quad.get_num_points();
  for (; elem != elem_end; ++elem, ++f_elem)
  {
    const int n_basis = elem->get_num_basis(DofProperties::active);
    DenseMatrix loc_mat(n_basis, n_basis);
    loc_mat = 0.0;

    DenseVector loc_rhs(n_basis);
    loc_rhs = 0.0;

    elem_handler->fill_element_cache(elem);
    f->fill_element_cache(f_elem);

    auto phi = elem->template get_basis<_Value, dim>(0,DofProperties::active);
    auto grad_phi  = elem->template get_basis<_Gradient, dim>(0,DofProperties::active);
    auto w_meas = elem->template get_w_measures<dim>(0);

    grad_phi.print_info(out);

    auto f_values = f_elem->template get_values<_Value,dim>(0);
    for (int i = 0; i < n_basis; ++i)
    {
      auto grad_phi_i = grad_phi.get_function_view(i);
      auto phi_i = phi.get_function_view(i);
      for (int j = 0; j < n_basis; ++j)
      {
        auto grad_phi_j = grad_phi.get_function_view(j);
        auto phi_j = phi.get_function_view(j);
        for (int qp = 0; qp < n_qp; ++qp)
          loc_mat(i,j) +=
            (scalar_product(grad_phi_i[qp], grad_phi_j[qp])
             +
             scalar_product(phi_i[qp], phi_j[qp])
            )
            * w_meas[qp];
      }

      for (int qp=0; qp<n_qp; ++qp)
        loc_rhs(i) += scalar_product(phi_i[qp], f_values[qp])
                      * w_meas[qp];
    }

    const auto loc_dofs = elem->get_local_to_global(DofProperties::active);
    matrix->add_block(loc_dofs, loc_dofs,loc_mat);
    rhs->add_block(loc_dofs, loc_rhs);
  }
  matrix->FillComplete();
  auto g = ConstantFunction<dim,0,1,1>::create(grid, IdentityFunction<dim>::create(grid), {0.});

  const QGauss<dim-1>  face_quad(deg);
  const std::set<boundary_id> dir_id {0};
  std::map<Index, Real> values;
  // TODO (pauletti, Mar 9, 2015): parametrize with dimension
  project_boundary_values<RefSpace>(
    const_pointer_cast<const Function>(g),
    space,
    face_quad,
    dir_id,
    values);
  apply_boundary_values(values, *matrix, *rhs, *solution);

  matrix->print_info(out);

  auto solver = create_solver(*matrix, *solution, *rhs);
  solver->solve();

  const int n_plot_points = deg+1;
  auto map1 = IdentityFunction<dim>::create(space->get_ptr_const_grid());
  Writer<dim> writer(map1, n_plot_points);


  using IgFunc = IgFunction<dim,0,1,1>;
  auto solution_function = IgFunc::create(space, solution);

  writer.template add_field<1,1>(solution_function, "solution");
  string filename = "poisson_problem-" + to_string(deg) + "-" + to_string(dim) + "d" ;
  writer.save(filename);

}


int main()
{
  for (int deg = 1; deg<3; ++deg)
  {
    const int n_knots = 5 + deg;
    // assemble_matrix<1>(n_knots, deg);
    assemble_matrix<2>(n_knots, deg);
  }
  return 0;
}
