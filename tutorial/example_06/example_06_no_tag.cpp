
#include <igatools/base/function_lib.h>
#include <igatools/basis_functions/bspline_space.h>
#include <igatools/basis_functions/bspline_element_accessor.h>
#include <igatools/base/quadrature_lib.h>

#include <igatools/linear_algebra/dense_matrix.h>
#include <igatools/linear_algebra/dense_vector.h>

#include <igatools/io/writer.h>

#include <igatools/basis_functions/space_tools.h>

#include <igatools/linear_algebra/distributed_matrix.h>
#include <igatools/linear_algebra/distributed_vector.h>
#include <igatools/linear_algebra/linear_solver.h>
#include <igatools/linear_algebra/dof_tools.h>

using namespace iga;
using namespace std;

using functions::ConstantFunction;
using space_tools::project_boundary_values;
using dof_tools::apply_boundary_values;
using dof_tools::get_sparsity_pattern;

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

private:
    using Space = BSplineSpace<dim>;
    shared_ptr<Space> space;

    const Quadrature<dim>   elem_quad;
    const Quadrature<dim-1> face_quad;

    shared_ptr<Matrix> matrix;
    shared_ptr<Vector> rhs;
    shared_ptr<Vector> solution;
};


template<int dim>
PoissonProblem<dim>::
PoissonProblem(const int n_knots, const int deg)
    :
    space(Space::create(CartesianGrid<dim>::create(n_knots), deg)),
    elem_quad(QGauss<dim>(deg+1)),
    face_quad(QGauss<dim-1>(deg+1))
{
    const auto n_basis = space->get_num_basis();
    matrix   = Matrix::create(get_sparsity_pattern(*space));
    rhs      = Vector::create(n_basis);
    solution = Vector::create(n_basis);
}



template<int dim>
void PoissonProblem<dim>::assemble()
{
    const int n_basis = space->get_num_basis_per_element();
    DenseMatrix loc_mat(n_basis, n_basis);
    DenseVector loc_rhs(n_basis);
    vector<Index> loc_dofs(n_basis);

    ConstantFunction<dim> f({5.});
    using ValueType = typename Function<dim>::ValueType;
    const int n_qp = elem_quad.get_num_points();
    vector<ValueType> f_values(n_qp);

    auto elem = space->begin();
    const auto elem_end = space->end();
    ValueFlags fill_flags = ValueFlags::value | ValueFlags::gradient |
                            ValueFlags::w_measure | ValueFlags::point;
    elem->init_values(fill_flags, elem_quad);

    for (; elem != elem_end; ++elem)
    {
        loc_mat = 0.;
        loc_rhs = 0.;
        elem->fill_values();

        auto points  = elem->get_points();
        auto phi     = elem->get_basis_values();
        auto grd_phi = elem->get_basis_gradients();
        auto w_meas  = elem->get_w_measures();

        f.evaluate(points, f_values);

        for (int i = 0; i < n_basis; ++i)
        {
            auto grd_phi_i = grd_phi.get_function_view(i);
            for (int j = 0; j < n_basis; ++j)
            {
                auto grd_phi_j = grd_phi.get_function_view(j);
                for (int qp = 0; qp < n_qp; ++qp)
                    loc_mat(i,j) +=
                        scalar_product(grd_phi_i[qp], grd_phi_j[qp])
                        * w_meas[qp];
            }
            auto phi_i = phi.get_function_view(i);

            for (int qp=0; qp<n_qp; ++qp)
                loc_rhs(i) += scalar_product(phi_i[qp], f_values[qp])
                              * w_meas[qp];
        }

        loc_dofs = elem->get_local_to_global();
        matrix->add_block(loc_dofs, loc_dofs,loc_mat);
        rhs->add_block(loc_dofs, loc_rhs);
    }

    matrix->fill_complete();

    ConstantFunction<dim> g({0.0});
    const boundary_id dir_id = 0;
    std::map<Index, Real> values;
    project_boundary_values<Space>(g, space, face_quad, dir_id, values);
    apply_boundary_values(values, *matrix, *rhs, *solution);
}


template<int dim>
void PoissonProblem<dim>::solve()
{
    LinearSolver solver(LinearSolver::Type::CG);
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
