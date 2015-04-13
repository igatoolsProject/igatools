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
#include <igatools/base/function_lib.h>
#include <igatools/base/identity_function.h>
#include <igatools/base/quadrature_lib.h>
#include <igatools/basis_functions/bspline_space.h>
#include <igatools/basis_functions/physical_space.h>
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
    using RefSpace = BSplineSpace<dim>;
    using Space    = PhysicalSpace<dim>;
    using Value = typename Function<dim>::Value;
    // [type aliases]


    shared_ptr<Space>        space;
    shared_ptr<MapFunction<dim, dim>> map;

    const Quadrature<dim>   elem_quad;
    const Quadrature<dim-1> face_quad;

    const boundary_id dir_id = 0;

    shared_ptr<Matrix> matrix;
    shared_ptr<Vector> rhs;
    shared_ptr<Vector> solution;

};



template<int dim>
PoissonProblem<dim>::
PoissonProblem(const int deg, const TensorSize<dim> &n_knots)
    :
    elem_quad(QGauss<dim>(deg+1)),
    face_quad(QGauss<dim-1>(deg+1))
{
    BBox<dim> box;
    box[0] = {{0.5,1}};
    for (int i=1; i<dim; ++i)
        box[i] = {{PI/4,PI/2}};

    auto grid = CartesianGrid<dim>::create(box, n_knots);
    auto ref_space = RefSpace::create(deg, grid);
    using Function = functions::BallFunction<dim>;
    map = Function::create(grid, IdentityFunction<dim>::create(grid));
    space = Space::create(ref_space, map);

    matrix = create_matrix(space);
    rhs = create_vector(space);
    solution=create_vector(space);

}



template<int dim>
void PoissonProblem<dim>::assemble()
{
    auto grid = space->get_grid();

    using Function = Function<dim,0,1,1>;
    using ConstFunction = ConstantFunction<dim,0,1,1>;
    using Value = typename Function::Value;

    Value b = {5.};
    auto f = ConstFunction::create(grid, IdentityFunction<dim>::create(grid), b);

    auto elem_handler = space->create_elem_handler();

    auto flag = ValueFlags::value | ValueFlags::gradient |
                ValueFlags::w_measure;

    elem_handler->reset(flag, elem_quad);

    f->reset(ValueFlags::value, elem_quad);

    auto f_elem = f->begin();
    auto elem   = space->begin();
    const auto elem_end = space->end();

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
        auto phi = elem->template get_basis<_Value, dim>(0,DofProperties::active);
        auto grad_phi  = elem->template get_basis<_Gradient, dim>(0,DofProperties::active);
        auto w_meas = elem->template get_w_measures<dim>(0);

        f->fill_element_cache(f_elem);
        auto f_values = f_elem->template get_values<_Value, dim>(0);

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

        const auto loc_dofs = elem->get_local_to_global(DofProperties::active);
        matrix->add_block(loc_dofs, loc_dofs,loc_mat);
        rhs->add_block(loc_dofs, loc_rhs);
    }

    matrix->FillComplete();

    // [dirichlet constraint]

    auto g = ConstFunction::
             create(grid, IdentityFunction<dim>::create(grid), {0.});


    const set<boundary_id> dir_id {0};
    std::map<Index, Real> values;
    // TODO (pauletti, Mar 9, 2015): parametrize with dimension
    project_boundary_values<Space>(
        const_pointer_cast<const Function>(g),
        space,
        face_quad,
        dir_id,
        values);
    apply_boundary_values(values, *matrix, *rhs, *solution);
    // [dirichlet constraint]
}


template<int dim>
void PoissonProblem<dim>::solve()
{
    auto solver = create_solver(matrix, solution, rhs);
    solver->solve();
}



template<int dim>
void PoissonProblem<dim>::output()
{
    const int n_plot_points = 2;
    auto map = space->get_map_func();
    Writer<dim> writer(map, n_plot_points);

    using IgFunc = IgFunction<Space>;
    auto solution_function = IgFunc::create(space, *solution);
    writer.template add_field<1,1>(solution_function, "solution");
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

    PoissonProblem<2> poisson_2d(deg, {n_knots, n_knots});
    poisson_2d.run();

    PoissonProblem<3> poisson_3d(deg, {n_knots, n_knots, n_knots});
    poisson_3d.run();

    return  0;
}
