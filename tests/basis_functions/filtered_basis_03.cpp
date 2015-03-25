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
 *  Test for solving on a filtered space
 *
 *  author: pauletti
 *  date: 2015-03-17
 *
 */

//TODO (pauletti, Mar 22, 2015): this is a development step, should be
// splitted and modified to became several unit testing afterwards
#include "../tests.h"

#include <igatools/base/quadrature_lib.h>
#include <igatools/basis_functions/bspline_space.h>
#include <igatools/basis_functions/bspline_element.h>
#include <igatools/linear_algebra/distributed_matrix.h>
#include <igatools/linear_algebra/distributed_vector.h>
#include <igatools/linear_algebra/linear_solver.h>
#include <igatools/base/identity_function.h>
#include <igatools/base/ig_function.h>
#include <igatools/io/writer.h>

struct DofProp
{
    static const std::string interior;
    static const std::string dirichlet;
    static const std::string neumman;
};

const std::string DofProp::interior = "interior";
const std::string DofProp::dirichlet  = "dirichlet";
const std::string DofProp::neumman  = "neumman";

static const LAPack la_pack = LAPack::trilinos_epetra;
using Mat = Matrix<la_pack>;
using Vec = Vector<la_pack>;

enum  bc : boundary_id
{
    dir, neu
};


template<int dim, int range = 1, int rank = 1>
void filtered_dofs(const int deg = 1, const int n_knots = 3)
{
    OUTSTART
    using RefSpace = ReferenceSpace<dim, range, rank>;
    using Space = BSplineSpace<dim, range, rank>;

    auto grid = CartesianGrid<dim>::create(n_knots);


    grid->set_boundary_id(3, bc::neu);


    auto space = Space::create(deg, grid);
    auto dof_dist = space->get_dof_distribution();
    dof_dist->add_dofs_property(DofProp::interior);
    dof_dist->add_dofs_property(DofProp::dirichlet);
    dof_dist->add_dofs_property(DofProp::neumman);




    std::set<Index> int_dofs= {4};
    dof_dist->set_dof_property_status(DofProp::interior, int_dofs,true);
    std::set<Index> dir_dofs= {6,3,0, 1, 2, 5, 8};
    dof_dist->set_dof_property_status(DofProp::dirichlet, dir_dofs,true);
    std::set<Index> neu_dofs= {7};
    dof_dist->set_dof_property_status(DofProp::neumman, neu_dofs,true);

    auto elem = space->begin();
    auto end  = space->end();

    auto space_manager =
        build_space_manager_single_patch<RefSpace>(space, DofProp::interior);
    auto matrix   = Mat::create(*space_manager);
    const auto dofs_set = space_manager->get_row_dofs();
//    const vector<Index> dofs_vec(dofs_set.begin(),dofs_set.end());
    auto vec      = Vec::create(dofs_set);
    auto solution     = Vec::create(dofs_set);
    matrix->print_info(out);
    vec->print_info(out);

    auto elem_handler = space->create_elem_handler();
    auto flag = ValueFlags::value | ValueFlags::gradient | ValueFlags::w_measure;
    QGauss<dim> elem_quad(deg+1);
    elem_handler->reset(flag, elem_quad);
    elem_handler->init_element_cache(elem);
    const int n_qp = elem_quad.get_num_points();
    for (; elem != end; ++elem)
    {
        const int n_basis = elem->get_num_basis(DofProp::interior);

        DenseMatrix loc_mat(n_basis, n_basis);
        loc_mat = 0.0;
        DenseVector loc_rhs(n_basis);
        loc_rhs = 0.0;

        elem_handler->fill_element_cache(elem);
        auto phi = elem->template get_values<0, dim>(0,DofProp::interior);
        auto grad_phi  = elem->template get_values<1, dim>(0,DofProp::interior);
        auto w_meas = elem->template get_w_measures<dim>(0);

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
                loc_rhs(i) += phi_i[qp][0] // f=1
                              * w_meas[qp];
        }

        const auto loc_dofs = elem->get_local_to_global(DofProp::interior);
        matrix->add_block(loc_dofs, loc_dofs, loc_mat);
        vec->add_block(loc_dofs, loc_rhs);
    }

    matrix->fill_complete();
    matrix->print_info(out);
    vec->print_info(out);

    using LinSolver = LinearSolverIterative<la_pack>;
    LinSolver solver(LinSolver::SolverType::CG);
    solver.solve(*matrix, *vec, *solution);

    const int n_plot_points = 4;
    auto map = IdentityFunction<dim>::create(space->get_grid());
    Writer<dim> writer(map, n_plot_points);


    using IgFunc = IgFunction<RefSpace>;
    auto solution_function = IgFunc::create(space,solution->as_ig_fun_coefficients(),
                                            DofProp::interior);

    writer.template add_field<1,1>(solution_function, "solution");
    string filename = "poisson_problem-" + to_string(dim) + "d" ;
    writer.save(filename);

    OUTEND
}




int main()
{
    const int dim = 2;
    filtered_dofs<dim>();

    return 0;
}
