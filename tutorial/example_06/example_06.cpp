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
#include <igatools/base/identity_function.h>
#include <igatools/base/function_lib.h>
// [functions]
// [old includes]
#include <igatools/basis_functions/bspline_space.h>
#include <igatools/basis_functions/bspline_element.h>
#include <igatools/base/quadrature_lib.h>

#include <igatools/linear_algebra/dense_matrix.h>
#include <igatools/linear_algebra/dense_vector.h>

#include <igatools/io/writer.h>
// [old includes]

// [project to boundary]
#include <igatools/basis_functions/space_tools.h>
// [project to boundary]

// [linear system]
#include <igatools/linear_algebra/distributed_matrix.h>
#include <igatools/linear_algebra/distributed_vector.h>
#include <igatools/linear_algebra/linear_solver.h>
#include <igatools/linear_algebra/sparsity_pattern.h>
#include <igatools/linear_algebra/dof_tools.h>
// [linear system]

using namespace iga;
using namespace std;

// [short names]
using functions::ConstantFunction;
using space_tools::project_boundary_values;
using dof_tools::apply_boundary_values;
// [short names]

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
    using Space = BSplineSpace<dim>;
    shared_ptr<Space> space;

    const QuadratureTensorProduct<dim>   elem_quad;
    const QuadratureTensorProduct<dim-1> face_quad;
    // [members]

    // [la members]
    shared_ptr<Matrix<LAPack::trilinos>> matrix;
    shared_ptr<Vector<LAPack::trilinos>> rhs;
    shared_ptr<Vector<LAPack::trilinos>> solution;
};
// [la members]


template<int dim>
PoissonProblem<dim>::
PoissonProblem(const int n_knots, const int deg)
    :
    space(Space::create(deg, CartesianGrid<dim>::create(n_knots))),
    elem_quad(QGauss<dim>(deg+1)),
    face_quad(QGauss<dim-1>(deg+1))
{
    const auto n_basis = space->get_num_basis();
    matrix   = Matrix<LAPack::trilinos>::create(*space->get_space_manager());
    rhs      = Vector<LAPack::trilinos>::create(n_basis);
    solution = Vector<LAPack::trilinos>::create(n_basis);
}



template<int dim>
void PoissonProblem<dim>::assemble()
{
    auto grid = space->get_grid();

    using Function = ConstantFunction<dim,0,1,1>;
    using Value = typename Function::Value;

    Value b = {5.};
    auto f = ConstantFunction<dim,0,1,1>::create(grid, IdentityFunction<dim>::create(grid), b);

    using ElementHandler = typename Space::ElementHandler;
    auto elem_handler = ElementHandler::create(space);

    auto flag = ValueFlags::value | ValueFlags::gradient |
                ValueFlags::w_measure;

    elem_handler->reset(flag, elem_quad);
    f->reset(ValueFlags::value, elem_quad);

    auto f_elem = f->begin();
    auto elem   = space->begin();
    const auto elem_end = space->end();

    const auto topology = Int<dim>();
    elem_handler->init_cache(elem,topology);
    f->init_cache(f_elem,topology);

    const int n_qp = elem_quad.get_num_points();

    for (; elem != elem_end; ++elem, ++f_elem)
    {
        const int n_basis = elem->get_num_basis();

        DenseMatrix loc_mat(n_basis, n_basis);
        loc_mat = 0.0;

        DenseVector loc_rhs(n_basis);
        loc_rhs = 0.0;

        elem_handler->fill_cache(elem,topology,0);
        auto phi = elem->template get_values<0, dim>(0);
        auto grad_phi  = elem->template get_values<1, dim>(0);
        auto w_meas = elem->template get_w_measures<dim>(0);

        f->fill_cache(f_elem, 0, topology);
        auto f_values = f_elem->template get_values<0,dim>(0);

        for (int i = 0; i < n_basis; ++i)
        {
            auto grad_phi_i = grad_phi.get_function_view(i);
            for (int j = 0; j < n_basis; ++j)
            {
                auto grad_phi_j = grad_phi.get_function_view(j);
                for (int qp = 0; qp < n_qp; ++qp)
                    loc_mat(i,j) +=
                        scalar_product(grad_phi_i[qp], grad_phi_j[qp])
                        * w_meas[qp];
            }
            auto phi_i = phi.get_function_view(i);

            for (int qp=0; qp<n_qp; ++qp)
                loc_rhs(i) += scalar_product(phi_i[qp], f_values[qp])
                              * w_meas[qp];
        }

        const auto loc_dofs = elem->get_local_to_global();
        matrix->add_block(loc_dofs, loc_dofs,loc_mat);
        rhs->add_block(loc_dofs, loc_rhs);
    }

    matrix->fill_complete();

    // [dirichlet constraint]

    auto g = Function::
             create(grid, IdentityFunction<dim>::create(grid), {0.});


    const boundary_id dir_id = 0;
    std::map<Index, Real> values;
    project_boundary_values<Space,LAPack::trilinos>(g, space, face_quad, dir_id, values);
    apply_boundary_values(values, *matrix, *rhs, *solution);
    // [dirichlet constraint]
}


template<int dim>
void PoissonProblem<dim>::solve()
{
    using LinSolver = LinearSolver<LAPack::trilinos>;
    LinSolver solver(LinSolver::SolverType::CG);
    solver.solve(*matrix, *rhs, *solution);
}


template<int dim>
void PoissonProblem<dim>::output()
{
    const int n_plot_points = 2;
    Writer<dim> writer(space->get_grid(), n_plot_points);

    writer.add_field(space, *solution, "solution");
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
    const int deg = 1 ;

    PoissonProblem<1> laplace_1d(n_knots, deg);
    laplace_1d.run();

    PoissonProblem<2> laplace_2d(n_knots, deg);
    laplace_2d.run();

    PoissonProblem<3> laplace_3d(n_knots, deg);
    laplace_3d.run();

    return  0;
}
