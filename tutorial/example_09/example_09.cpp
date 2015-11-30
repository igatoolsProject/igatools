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

// [old includes]
#include <igatools/geometry/grid_function_lib.h>
#include <igatools/functions/function_lib.h>
#include <igatools/base/quadrature_lib.h>
#include <igatools/basis_functions/bspline.h>
#include <igatools/basis_functions/physical_space_basis.h>
#include <igatools/basis_functions/physical_space_element.h>
#include <igatools/basis_functions/space_tools.h>
#include <igatools/linear_algebra/dense_matrix.h>
#include <igatools/linear_algebra/dense_vector.h>
#include <igatools/linear_algebra/epetra_solver.h>
#include <igatools/linear_algebra/dof_tools.h>
#include <igatools/io/writer.h>
// [old includes]

// [unqualified names]
using namespace iga;
using namespace std;
using namespace EpetraTools;
using functions::ConstantFunction;
using space_tools::project_boundary_values;
using dof_tools::apply_boundary_values;
using numbers::PI;
// [unqualified names]

// [Problem class]
template<int dim>
class PoissonProblem
{
public:
  PoissonProblem(const int deg, const TensorSize<dim> &n_knots);
  void run();

private:
  void assemble();
  void solve();
  void output();
  // [Problem class]

  // [type aliases]
private:
  using RefSpace = BSpline<dim>;
  using Basis    = PhysicalSpaceBasis<dim>;
  using Value = typename Function<dim>::Value;
  // [type aliases]


  shared_ptr<const Basis>        basis;
//  shared_ptr<MapFunction<dim, dim>> domain;

  const std::shared_ptr<const Quadrature<dim>>   elem_quad;
  const std::shared_ptr<const Quadrature<dim-1>> face_quad;

  const boundary_id dir_id = 0;

  shared_ptr<Matrix> matrix;
  shared_ptr<Vector> rhs;
  shared_ptr<Vector> solution;

};



template<int dim>
PoissonProblem<dim>::
PoissonProblem(const int deg, const TensorSize<dim> &n_knots)
  :
  elem_quad(QGauss<dim>::create(deg+1)),
  face_quad(QGauss<dim-1>::create(deg+1))
{
  BBox<dim> box;
  box[0] = {{0.5,1}};
  for (int i=1; i<dim; ++i)
    box[i] = {{PI/4,PI/2}};

  auto grid = Grid<dim>::create(box, n_knots);
  auto space = SplineSpace<dim>::create(deg, grid);
  auto ref_basis = BSpline<dim>::create(space);
  auto ball_func = grid_functions::BallGridFunction<dim>::create(grid);
  auto ball_domain = Domain<dim>::create(ball_func);
  basis = Basis::create(ref_basis, ball_domain);

  matrix = create_matrix(*basis,DofProperties::active,Epetra_SerialComm());
  rhs = create_vector(matrix->RangeMap());
  solution=create_vector(matrix->DomainMap());

}



template<int dim>
void PoissonProblem<dim>::assemble()
{
//  auto grid = basis->get_grid();

  const auto domain = basis->get_physical_domain();

  using ConstFunc = functions::ConstantFunction<dim,0,1,1>;
  auto f = ConstFunc::const_create(domain, {5.});

  auto basis_elem_handler = basis->create_cache_handler();
  auto f_elem_handler = f->create_cache_handler();

  using Flags = space_element::Flags;
  auto flag = Flags::value |
              Flags::gradient |
              Flags::w_measure;

  basis_elem_handler->set_element_flags(flag);

  f_elem_handler->set_element_flags(function_element::Flags::D0);

  auto f_elem = f->begin();
  auto basis_elem   = basis->begin();
  const auto elem_end = basis->end();

  basis_elem_handler->init_cache(basis_elem,elem_quad);
  f_elem_handler->init_cache(f_elem,elem_quad);

//  const int n_qp = elem_quad->get_num_points();

  for (; basis_elem != elem_end; ++basis_elem, ++f_elem)
  {
    f_elem_handler->fill_element_cache(f_elem);
    const auto &f_values = f_elem->get_element_values_D0();

    basis_elem_handler->fill_element_cache(basis_elem);
    const auto loc_mat = basis_elem->integrate_gradu_gradv();
    const auto loc_rhs = basis_elem->integrate_u_func(f_values);

    const auto loc_dofs = basis_elem->get_local_to_global(DofProperties::active);
    matrix->add_block(loc_dofs, loc_dofs,loc_mat);
    rhs->add_block(loc_dofs, loc_rhs);
  }

  matrix->FillComplete();

  // [dirichlet constraint]


  using SubGridElemMap = typename Grid<dim>::template SubGridMap<dim-1>;
  using FuncAtBndry = Function<dim-1,1,1>;
  SafeSTLMap<int,std::shared_ptr<const FuncAtBndry>> bndry_funcs;
  for (int face_id = 0 ; face_id < UnitElement<dim>::n_faces ; ++face_id)
  {
    SubGridElemMap sub_grid_elem_map;
    const auto sub_grid = domain->get_grid_function()->get_grid()->template get_sub_grid<dim-1>(face_id,sub_grid_elem_map);

    const auto sub_domain = domain->template get_sub_domain<dim-1>(face_id,sub_grid_elem_map,sub_grid);
    bndry_funcs[face_id] = functions::ConstantFunction<dim-1,1,1>::const_create(sub_domain, {0.});
  }

  SafeSTLMap<Index, Real> bndry_values;
  space_tools::project_boundary_values(bndry_funcs,*basis,face_quad,bndry_values);

  apply_boundary_values(bndry_values, *matrix, *rhs, *solution);
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
  auto domain = basis->get_physical_domain();
  Writer<dim> writer(domain, n_plot_points);

  using IgFunc = IgFunction<dim,0,1,1>;
  auto solution_function = IgFunc::const_create(basis, *solution);
  writer.template add_field<1,1>(*solution_function, "solution");

  string filename = "poisson_problem-" + to_string(dim) + "d" ;
  writer.save(filename);
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
  const int deg     = 1;

  PoissonProblem<1> poisson_1d(deg, {n_knots});
  poisson_1d.run();
  //*/

  PoissonProblem<2> poisson_2d(deg, {n_knots, n_knots});
  poisson_2d.run();
//*/

  PoissonProblem<3> poisson_3d(deg, {n_knots, n_knots, n_knots});
  poisson_3d.run();

  return  0;
}
