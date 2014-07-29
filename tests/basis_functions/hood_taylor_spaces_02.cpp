//-+--------------------------------------------------------------------
// Igatools a general purpose Isogeometric analysis library.
// Copyright (C) 2012-2014  by the igatools authors (see authors.txt).
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

// TODO (pauletti, Jun 11, 2014): put appropriate header

#include "../tests.h"

#include <igatools/basis_functions/bspline_space.h>
#include <igatools/basis_functions/bspline_element_accessor.h>

#include <igatools/base/quadrature_lib.h>
#include <igatools/linear_algebra/dense_matrix.h>
#include <igatools/linear_algebra/distributed_matrix.h>
#include <igatools/linear_algebra/dof_tools.h>

template <int dim>
class StokesProblem
{
public:
    StokesProblem(const int deg, const int n_knots);
    void run();
private:
    void assemble_Bt();

private:
    using PreSpace = BSplineSpace<dim>;
    using VelSpace = BSplineSpace<dim, dim>;

    shared_ptr<PreSpace> pre_space_;
    shared_ptr<VelSpace> vel_space_;

    const Quadrature<dim> elem_quad_;

#if defined(USE_TRILINOS)
    const static auto la_pack = LAPack::trilinos;
#elif defined(USE_PETSC)
    const static auto la_pack = LAPack::petsc;
#endif

    shared_ptr<Matrix<la_pack>> Bt_;
};

template <int dim>
void StokesProblem<dim>::assemble_Bt()
{
    const int vel_n_basis = vel_space_->get_num_basis_per_element();
    const int pre_n_basis = pre_space_->get_num_basis_per_element();
    DenseMatrix loc_mat(vel_n_basis, pre_n_basis);
    vector<Index> vel_loc_dofs(vel_n_basis);
    vector<Index> pre_loc_dofs(pre_n_basis);

    auto vel_el = vel_space_->begin();
    auto pre_el = pre_space_->begin();
    auto end_el = vel_space_->end();

    ValueFlags vel_flag = ValueFlags::divergence | ValueFlags::w_measure;
    ValueFlags pre_flag = ValueFlags::value;
    vel_el->init_cache(vel_flag,elem_quad_);
    pre_el->init_cache(pre_flag,elem_quad_);

    const int n_qp = elem_quad_.get_num_points();
    for (; vel_el != end_el; ++vel_el, ++pre_el)
    {
        loc_mat.clear();
        vel_el->fill_cache();
        pre_el->fill_cache();

        auto q      = pre_el->get_basis_values();
        auto div_v  = vel_el->get_basis_divergences();
        auto w_meas = vel_el->get_w_measures();

        for (int i=0; i<vel_n_basis; ++i)
        {
            auto div_i = div_v.get_function_view(i);
            for (int j=0; j<pre_n_basis; ++j)
            {
                auto q_j = q.get_function_view(j);
                for (int qp=0; qp<n_qp; ++qp)
                    loc_mat(i,j) -=  scalar_product(div_i[qp], q_j[qp])
                                     * w_meas[qp];
            }
        }
        vel_loc_dofs = vel_el->get_local_to_global();
        pre_loc_dofs = pre_el->get_local_to_global();
        Bt_->add_block(vel_loc_dofs, pre_loc_dofs, loc_mat);

        out << loc_mat << endl;
    }
    Bt_->fill_complete();
}

template <int dim>
void StokesProblem<dim>::run()
{
    assemble_Bt();
}


template <int dim>
StokesProblem<dim>::
StokesProblem(const int deg, const int n_knots)
    :
    elem_quad_(QGauss<dim>(deg+3))
{
    const int reg = 0;

    typename PreSpace::Multiplicity pre_mult;
    typename VelSpace::Multiplicity vel_mult;

    vector<Index> mult_p(n_knots-2, deg - reg);
    vector<Index> mult_v(n_knots-2, deg + 1 -reg);
    for (int i = 0; i < dim; ++i)
    {
        pre_mult.copy_data_direction(i, mult_p) ;
        vel_mult.copy_data_direction(i, mult_v) ;
    }


    auto pres_m = make_shared<typename PreSpace::MultiplicityTable>(pre_mult);
    auto vel_m  = make_shared<typename VelSpace::MultiplicityTable>(vel_mult);

    typename PreSpace::DegreeTable pre_deg;//(TensorIndex<dim>(deg));
    pre_deg(0) = TensorIndex<dim>(deg);

    const typename VelSpace::DegreeTable vel_deg(TensorIndex<dim>(deg+1));

    auto grid = CartesianGrid<dim>::create(n_knots);
    pre_space_ = PreSpace::create(pre_deg, grid, pres_m);

    vel_space_ = VelSpace::create(vel_deg, grid,vel_m);

    const SparsityPattern sparsity_pattern(*vel_space_->get_dofs_manager(),*pre_space_->get_dofs_manager());

#if defined(USE_TRILINOS)
    const auto la_pack = LAPack::trilinos;
#elif defined(USE_PETSC)
    const auto la_pack = LAPack::petsc;
#endif

    Bt_ = Matrix<la_pack>::create(sparsity_pattern);
}


int main()
{
    out.depth_console(10);
    StokesProblem<2> stokes(1,3);
    stokes.run();

    return 0;
}
