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

#include <igatools/basis_functions/bspline.h>
#include <igatools/base/logstream.h>

#include <igatools/basis_functions/nurbs.h>
#include <igatools/functions/grid_function_lib.h>
#include <igatools/functions/ig_function.h>
#include <igatools/io/writer.h>
#include <igatools/io/objects_container_xml_writer.h>

#include <igatools/geometry/grid_element.h>
#include <igatools/basis_functions/bspline_element.h>
#include <igatools/basis_functions/bspline_handler.h>
#include <igatools/base/quadrature_lib.h>
#include <igatools/base/objects_container.h>

// [include]
#include <igatools/linear_algebra/dof_tools.h>
#include <igatools/linear_algebra/epetra_matrix.h>
#include <igatools/linear_algebra/epetra_vector.h>
#include <igatools/linear_algebra/epetra_solver.h>
// [include]

using namespace iga;
using namespace std;
// [using]
using namespace EpetraTools;
// [using]

LogStream out;

template<int dim>
class PoissonProblem
{

private:
  shared_ptr<const Grid<dim>>        grid;
  shared_ptr<const SplineSpace<dim>> space;
  shared_ptr<const BSpline<dim>>     basis;
  shared_ptr<const QGauss<dim>>      quad;
// [system]
  shared_ptr<Matrix> mat;
  shared_ptr<Vector> rhs;
  shared_ptr<Vector> sol;
// [system]

  using IgGridFunc_t = IgGridFunction<dim,1>;

public:
  PoissonProblem(const Size nel, const Index deg)
  {
    grid  = Grid<dim>::const_create(nel+1);
    space = SplineSpace<dim>::const_create(deg,grid);
    basis = BSpline<dim>::const_create(space);
    quad  = QGauss<dim>::const_create(deg+1);
// [sys_create]
    mat   = create_matrix(*basis,DofProperties::active,Epetra_SerialComm());
    rhs   = create_vector(mat->RangeMap());
    sol   = create_vector(mat->DomainMap());
// [sys_create]
  };

  void assemble();
  void solve();
  std::shared_ptr<const IgGridFunc_t> get_solution();
  void save();
  void run();
};

template<int dim>
void PoissonProblem<dim>::assemble()
{

  auto basis_el      = basis->begin();
  auto basis_el_end  = basis->end();
  auto cache_handler = basis->create_cache_handler();
  auto flag = basis_element::Flags::value |
              basis_element::Flags::gradient |
              basis_element::Flags::w_measure;
  cache_handler->set_element_flags(flag);
  cache_handler->init_element_cache(basis_el,quad);

// [loop_init]
  const auto num_quad  = quad->get_num_points();
  const auto num_basis = basis_el->get_num_basis(DofProperties::active);
  const auto source = 1.0;
  for (; basis_el!=basis_el_end; ++basis_el)
  {
    cache_handler->fill_element_cache(basis_el);

    DenseMatrix loc_mat(num_basis,num_basis);
    loc_mat = 0.0;
    DenseVector loc_rhs(num_basis);
    loc_rhs = 0.0;

    auto values = basis_el->get_element_values();
    auto grads  = basis_el->get_element_gradients();
    auto w_meas = basis_el->get_element_w_measures();
// [loop_init]

// [stiffness]
    for (int i=0; i<num_basis; i++)
    {
      const auto &grd_i = grads.get_function_view(i);
      for (int j=0; j<num_basis; j++)
      {
        const auto &grd_j = grads.get_function_view(j);
        for (int q=0; q<num_quad; q++)
        {
          loc_mat(i,j) += scalar_product(grd_i[q], grd_j[q]) * w_meas[q];
        }
      }
    }
// [stiffness]

// [rhs]
    for (int i=0; i<num_basis; i++)
    {
      const auto &vals = values.get_function_view(i);
      for (int q=0; q<num_quad; q++)
      {
        loc_rhs(i) += source * vals[q][0] * w_meas[q];
      }
    }
// [rhs]

// [add_block]
    const auto loc_dofs = basis_el->get_local_to_global();
    mat->add_block(loc_dofs, loc_dofs,loc_mat);
    rhs->add_block(loc_dofs, loc_rhs);
  }
  mat->FillComplete();
// [add_block]

// [boundary]
  auto dof_distribution = basis->get_dof_distribution();
  Topology<dim-1> sub_elem_topology;
  std::map<Index,Real> bdr_vals;

  for (int face=0; face<2*dim; face++)
  {
    auto bdr_dofs = dof_distribution->get_boundary_dofs(face,sub_elem_topology);

    for (set<Index>::iterator it=bdr_dofs.begin(); it!=bdr_dofs.end(); it++)
    {
      bdr_vals.insert(std::pair<Index,Real>(*it,0.0));
    }
  }
  dof_tools::apply_boundary_values(bdr_vals,*mat,*rhs,*sol);
}
// [boundary]

// [solve]
template<int dim>
void PoissonProblem<dim>::solve()
{
  auto solver = create_solver(*mat,*sol,*rhs);
  solver->solve();
}
// [solve]

template<int dim>
auto
PoissonProblem<dim>::get_solution() ->
std::shared_ptr<const IgGridFunc_t>
{
  return IgGridFunc_t::const_create(basis, *sol);
}

template<int dim>
void PoissonProblem<dim>::save()
{

  const auto solution = this->get_solution();
  string filename = "problem_" + to_string(dim) + "d" ;

  const auto solution_non_const = std::const_pointer_cast<IgGridFunc_t>(solution);
  solution_non_const->set_name("solution");

#ifdef SERIALIZATION
  const auto objs_container = ObjectsContainer::create();
  objs_container->insert_const_object<GridFunction<dim, 1>>(solution);
  {
    std::ofstream xml_ostream(filename + ".iga");
    OArchive xml_out(xml_ostream);
    xml_out << *objs_container;
  }
#else
#ifdef XML_IO
  const auto objs_container = ObjectsContainer::create();
  objs_container->insert_const_object<GridFunction<dim, 1>>(solution);
  ObjectsContainerXMLWriter::write(filename + ".iga", objs_container);
#else
  Writer<dim> writer(grid,5);
  writer.add_field(*solution, "solution");
  writer.save(filename);
#endif
#endif
}

template<int dim>
void PoissonProblem<dim>::run()
{
  assemble();
  solve();
  save();
}

int main()
{
  const int nel = 4;
  const int deg = 2;

  PoissonProblem<1> problem_1d(nel,deg);
  problem_1d.run();

  PoissonProblem<2> problem_2d(nel,deg);
  problem_2d.run();

  PoissonProblem<3> problem_3d(nel,deg);
  problem_3d.run();

  return 0;
}
