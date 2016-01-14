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

// [functions]
//#include <igatools/functions/identity_function.h>
#include <igatools/functions/grid_function_lib.h>
// [functions]
// [old includes]
#include <igatools/basis_functions/bspline.h>
#include <igatools/basis_functions/bspline_element.h>
#include <igatools/basis_functions/bspline_element_handler.h>
#include <igatools/base/quadrature_lib.h>

//#include <igatools/linear_algebra/dense_matrix.h>
//#include <igatools/linear_algebra/dense_vector.h>

#include <igatools/io/writer.h>
// [old includes]

// [project to boundary]
#include <igatools/basis_functions/space_tools.h>
// [project to boundary]

// [linear system]
#include <igatools/linear_algebra/epetra_solver.h>
#include <igatools/linear_algebra/dof_tools.h>
// [linear system]

using namespace iga;
using namespace std;

using namespace EpetraTools;


// [class functions]
template<int dim>
class PoissonProblem
{
public:
  PoissonProblem(const int n_knots, const int deg);
  void run();

private:
  void assemble();
  void solve();
  void output();
  // [class functions]

  // [members]
private:
  using RefSpace = ReferenceBasis<dim>;
  using Basis = BSpline<dim>;
  shared_ptr<const Basis> basis;

  shared_ptr<const Quadrature<dim>>   elem_quad;
  shared_ptr<const Quadrature<dim-1>> face_quad;
  // [members]

  // [la members]


  shared_ptr<Matrix> matrix;
  shared_ptr<Vector> rhs;
  shared_ptr<Vector> solution;
};
// [la members]


template<int dim>
PoissonProblem<dim>::
PoissonProblem(const int n_knots, const int deg)
  :
  basis(Basis::const_create(
         SplineSpace<dim>::const_create(
           deg, Grid<dim>::const_create(n_knots)))),
  elem_quad(QGauss<dim>::create(deg+1)),
  face_quad(QGauss<dim-1>::create(deg+1)),
  matrix(create_matrix(*basis,DofProperties::active, Epetra_SerialComm())),
  rhs(create_vector(matrix->RangeMap())),
  solution(create_vector(matrix->DomainMap()))
{}



template<int dim>
void PoissonProblem<dim>::assemble()
{
  auto grid = basis->get_grid();

  using ConstFunction = grid_functions::ConstantGridFunction<dim,1>;

  auto func = ConstFunction::const_create(grid, {5.});

  auto basis_elem_handler = basis->create_cache_handler();

  using BasisFlags = space_element::Flags;
  auto basis_flags = BasisFlags::value | BasisFlags::gradient |
                     BasisFlags::w_measure;
  basis_elem_handler->set_element_flags(basis_flags);

  auto func_elem_handler = func->create_cache_handler();
  using FuncFlags = grid_function_element::Flags;
  func_elem_handler->set_element_flags(FuncFlags::D0);

  auto func_elem = func->begin();

  auto basis_elem   = basis->begin();
  const auto basis_elem_end = basis->end();

  basis_elem_handler->init_element_cache(basis_elem,elem_quad);
  func_elem_handler->init_element_cache(func_elem,elem_quad);

  for (; basis_elem != basis_elem_end; ++basis_elem, ++func_elem)
  {
    basis_elem_handler->fill_element_cache(basis_elem);
    func_elem_handler->fill_element_cache(func_elem);

    const auto &f_values = func_elem->get_element_values_D0();

    const auto loc_mat = basis_elem->integrate_gradu_gradv();
    const auto loc_rhs = basis_elem->integrate_u_func(f_values);

    const auto loc_dofs = basis_elem->get_local_to_global();
    matrix->add_block(loc_dofs, loc_dofs,loc_mat);
    rhs->add_block(loc_dofs, loc_rhs);
  }

  matrix->FillComplete();


  // [dirichlet constraint]
  const auto g = ConstFunction::const_create(grid, {0.});

//  LogStream out;
  using SubGridElemMap = typename Grid<dim>::template SubGridMap<dim-1>;
  using SubFunc = SubGridFunction<dim-1,dim,1>;
  SafeSTLMap<int,std::shared_ptr<const SubFunc>> bndry_funcs;
  for (int face_id = 0 ; face_id < UnitElement<dim>::n_faces ; ++face_id)
  {
    SubGridElemMap sub_grid_elem_map;
    const auto sub_grid = grid->template get_sub_grid<dim-1>(face_id,sub_grid_elem_map);

    bndry_funcs[face_id] = g->template get_sub_function<dim-1>(face_id,sub_grid_elem_map,sub_grid);
  }

  SafeSTLMap<Index, Real> bndry_values;
  space_tools::project_boundary_values<dim,1>(bndry_funcs,*basis,face_quad,bndry_values);

  dof_tools::apply_boundary_values(bndry_values, *matrix, *rhs, *solution);
  // [dirichlet constraint]
}


template<int dim>
void PoissonProblem<dim>::solve()
{
  auto solver = create_solver(*matrix, *solution, *rhs);
  solver->solve();
}


template<int dim>
void PoissonProblem<dim>::output()
{
  const int n_plot_points = 2;
  Writer<dim> writer(basis->get_grid(), n_plot_points);


  using IgFunc = IgGridFunction<dim,1>;
  auto solution_function = IgFunc::const_create(basis, *solution);
  writer.template add_field<1>(*solution_function, "solution");
  string filename = "poisson_problem-" + to_string(dim) + "d" ;
  writer.save(filename, true);
}


template<int dim>
void PoissonProblem<dim>::run()
{
  assemble();
  solve();
  output();
}


int main()
{
  const int n_knots = 10;
  const int deg = 1 ;


  PoissonProblem<1> laplace_1d(n_knots, deg);
  laplace_1d.run();

  PoissonProblem<2> laplace_2d(n_knots, deg);
  laplace_2d.run();

  PoissonProblem<3> laplace_3d(n_knots, deg);
  laplace_3d.run();

  return  0;
}
