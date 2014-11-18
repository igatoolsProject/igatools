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

#ifndef SPACE_TOOLS_H_
#define SPACE_TOOLS_H_

#include <igatools/base/new_function.h>
#include <igatools/base/ig_function.h>

#include <igatools/linear_algebra/distributed_matrix.h>
#include <igatools/linear_algebra/linear_solver.h>

#include<set>
#include <igatools/geometry/mapping_slice.h>

#include <igatools/basis_functions/physical_space_element.h>
#include <igatools/basis_functions/space_element_handler.h>

IGA_NAMESPACE_OPEN
namespace space_tools
{
/**
 * Perform an (L2)-Projection the function @p func
 * onto the space @p space using the quadrature rule @p quad.
 *  The projection is a numerical vector (the coefficients of
 *  the projected function)
 */
template<class Space, LAPack la_pack = LAPack::trilinos>
std::shared_ptr<IgFunction<Space> >
projection_l2(const std::shared_ptr<const typename Space::Func> function,
              std::shared_ptr<const Space> space,
              const Quadrature<Space::dim> &quad)
{
    auto func = function->clone();
    const int dim = Space::dim;

    const auto space_manager = space->get_space_manager();
    Matrix<la_pack> matrix(*space_manager);

    const auto space_dofs_set = space_manager->get_row_dofs();
    vector<Index> space_dofs(space_dofs_set.begin(),space_dofs_set.end());
    Vector<la_pack> rhs(space_dofs);
    Vector<la_pack> sol(space_dofs);

    auto func_flag = NewValueFlags::point | NewValueFlags::value;
    func->reset(func_flag, quad);

    typename Space::ElementHandler sp_filler(space);
    auto sp_flag = NewValueFlags::point | NewValueFlags::value| NewValueFlags::w_measure;
    sp_filler.template reset<dim>(sp_flag, quad);

    auto f_elem = func->begin();
    auto elem = space->begin();
    auto end  = space->end();

    func->init_cache(f_elem, Int<dim>());
    sp_filler.template init_cache<dim>(elem);

    const int n_qp = quad.get_num_points();
    const int n_basis = elem->get_num_basis();
    DenseVector loc_rhs(n_basis);
    DenseMatrix loc_mat(n_basis, n_basis);


    for (; elem != end; ++elem, ++f_elem)
    {
        func->fill_cache(f_elem, 0, Int<dim>());
        sp_filler.template fill_cache<dim>(elem, 0);

        loc_mat = 0.;
        loc_rhs = 0.;

        auto f_at_qp = f_elem->template get_values<0,dim>(0);
        auto phi = elem->template get_values<0,dim>(0);



        // computing the upper triangular part of the local matrix
        auto w_meas = elem->template get_w_measures<dim>(0);
        for (int i = 0; i < n_basis; ++i)
        {
            const auto phi_i = phi.get_function_view(i);
            for (int j = i; j < n_basis; ++j)
            {
                const auto phi_j = phi.get_function_view(j);
                for (int q = 0; q < n_qp; ++q)
                    loc_mat(i,j) += scalar_product(phi_i[q], phi_j[q]) * w_meas[q];
            }

            for (int q = 0; q < n_qp; q++)
                loc_rhs(i) += scalar_product(f_at_qp[q], phi_i[q]) * w_meas[q];
        }

        // filling symmetric ;lower part of local matrix
        for (int i = 0; i < n_basis; ++i)
            for (int j = 0; j < i; ++j)
                loc_mat(i, j) = loc_mat(j, i);

        const auto local_dofs = elem->get_local_to_global();
        matrix.add_block(local_dofs,local_dofs,loc_mat);
        rhs.add_block(local_dofs,loc_rhs);
    }
    matrix.fill_complete();

    // TODO (pauletti, Oct 9, 2014): the solver must use a precon
    const Real tolerance = 1.0e-15;
    const int max_num_iter = 1000;
    using LinSolver = LinearSolver<la_pack>;
    LinSolver solver(LinSolver::SolverType::CG,tolerance,max_num_iter);
    solver.solve(matrix, rhs, sol);

    return std::make_shared<IgFunction<Space>>(IgFunction<Space>(space, sol));

}



/**
 * Projects (using the L2 scalar product) a function to the whole or part
 * of the boundary of the domain.
 * The piece of the domain is indicated by the boundary ids and the
 * projection is computed using the provided quadrature rule.
 *
 * The projected function is returned in boundary_values, a map containing all
 * indices of degrees of freedom at the boundary and the computed coefficient value
 * for this degree of freedom.
 *
 */
template<class Space, LAPack la_pack = LAPack::trilinos>
void
project_boundary_values(const std::shared_ptr<const typename Space::Func> function,
                        std::shared_ptr<const Space> space,
                        const Quadrature<Space::dim-1> &quad,
                        const std::set<boundary_id>  &boundary_ids,
                        std::map<Index, Real>  &boundary_values)
{
    const int dim   = Space::dim;
    const int range = Space::range;
    const int rank  = Space::rank;
    const int codim = Space::codim;
    //const int space_dim = Space::space_dim;

    const int sub_dim = dim - 1;
    using GridType = typename Space::GridType;
    using SubSpace = typename Space::template SubSpace<sub_dim>;
    using InterSpaceMap = typename Space::template InterSpaceMap<sub_dim>;
    using SubFunc = SubFunction<sub_dim, dim, codim, range, rank>;


    auto grid = space->get_grid();

    std::set<int> sub_elems;
    auto bdry_begin = boundary_ids.begin();
    auto bdry_end   = boundary_ids.end();
    for (auto &s_id : UnitElement<Space::dim>::template elems_ids<sub_dim>())
    {
        const auto bdry_id = grid->get_boundary_id(s_id);
        if (find(bdry_begin, bdry_end, bdry_id) != bdry_end)
            sub_elems.insert(s_id);
    }

    for (const Index &s_id : sub_elems)
    {
        using  InterGridMap = typename GridType::template InterGridMap<sub_dim>;
        auto elem_map = std::make_shared<InterGridMap>(InterGridMap());

        auto grid = space->get_grid();
        auto sub_grid = grid->template get_sub_grid<sub_dim>(s_id, *elem_map);

        InterSpaceMap  dof_map;
        auto sub_space = space->template get_sub_space<sub_dim>(s_id, dof_map, sub_grid, elem_map);
        auto sub_func = SubFunc::create(sub_grid, function, s_id, *elem_map);

        auto proj = projection_l2<SubSpace,la_pack>(sub_func, sub_space, quad);

        const auto coef = proj->get_coefficients();
        const int face_n_dofs = dof_map.size();
        for (Index i = 0; i< face_n_dofs; ++i)
            boundary_values[dof_map[i]] = coef(i);
    }
}




template<int dim, int codim = 0, int range = 1, int rank = 1>
Real integrate_difference(NewFunction<dim, codim, range, rank> &f,
                          NewFunction<dim, codim, range, rank> &g,
                          const Quadrature<dim> &quad,
                          const Norm &norm_flag,
                          vector<Real> &element_error)
{
    using Func = NewFunction<dim, codim, range, rank>;

    bool is_L2_norm     = contains(norm_flag, Norm::L2);
    bool is_H1_norm     = contains(norm_flag, Norm::H1);
    bool is_H1_seminorm = contains(norm_flag, Norm::H1_semi);

    Assert(is_L2_norm || is_H1_seminorm || is_H1_norm,
           ExcMessage("No active flag for the error norm."));


    Assert(!((is_L2_norm && is_H1_seminorm) ||
             (is_L2_norm && is_H1_norm) ||
             (is_H1_seminorm && is_H1_norm)),
           ExcMessage("Only a single flag for the error norm can be used."));


    if (is_H1_norm)
    {
        is_L2_norm     = true;
        is_H1_seminorm = true;
    }

    auto flag = NewValueFlags::point | NewValueFlags::w_measure;

    if (is_L2_norm)
        flag |= NewValueFlags::value;

    if (is_H1_seminorm)
        flag |= NewValueFlags::gradient;

    f.reset(flag, quad);
    g.reset(flag, quad);
    const int n_points   =  quad.get_num_points();
    //const int n_elements =  element_error.size();
    auto elem_f = f.begin();
    auto elem_g = g.begin();
    auto end = f.end();

    f.init_cache(elem_f, Int<dim>());
    g.init_cache(elem_g, Int<dim>());
    typename Func::Value err;
    for (; elem_f != end; ++elem_f, ++elem_g)
    {
        f.fill_cache(elem_f, 0, Int<dim>());
        g.fill_cache(elem_g, 0, Int<dim>());

        const int elem_id = elem_f->get_flat_index();
        element_error[ elem_id ] = 0.0;

        if (is_L2_norm)
        {
            auto f_val = elem_f->template get_values<0,dim>(0);
            auto g_val = elem_g->template get_values<0,dim>(0);
            auto w_meas = elem_f->template get_w_measures<dim>(0);

            Real element_err_L2_pow2 = 0.0;
            for (int iPt = 0; iPt < n_points; ++iPt)
            {
                err = f_val[iPt] - g_val[iPt];
                element_err_L2_pow2 += err.norm_square() * w_meas[iPt];
            }
            element_error[ elem_id ] += element_err_L2_pow2;
        }
        element_error[ elem_id ] = sqrt(element_error[ elem_id ]);
    }

#if 0
    using Value = typename Space::Func::Value;
    using Gradient = typename Space::Func::Gradient;

    ValueVector<Value>    u(n_points);
    ValueVector<Gradient> grad_u(n_points);

    Value err;
    Gradient grad_err;

    vector<Real>     norm_err_L2_square(n_elements);
    vector<Real> seminorm_err_H1_square(n_elements);

    for (; elem != end; ++elem)
    {



        if (is_H1_seminorm)
        {
            const auto &grad_uh = elem->evaluate_field_gradients(solution_coefs_elem);
            exact_solution.evaluate_gradients(map_at_points, grad_u);

            Real element_err_semiH1_pow2 = 0.0;
            for (int iPt = 0; iPt < n_points; ++iPt)
            {
                grad_err = grad_uh[iPt] - grad_u[iPt];

                element_err_semiH1_pow2 += grad_err.norm_square() * elem->get_w_measures()[iPt];
            }
            element_error[ elem_id ] += element_err_semiH1_pow2;
        }


    }
#endif
    Real err_pow2 = 0.0;
    for (const Real &elem_err : element_error)
        err_pow2 += elem_err * elem_err;

    Real total_error = sqrt(err_pow2);

    return total_error;

}



};




IGA_NAMESPACE_CLOSE

#if 0
#include <igatools/base/config.h>
#include <igatools/base/quadrature.h>
#include <igatools/geometry/cartesian_grid.h>
#include <igatools/linear_algebra/distributed_matrix.h>
#include <igatools/geometry/mapping.h>
#include <igatools/geometry/mapping_slice.h>
#include <igatools/base/function.h>

#include <igatools/basis_functions/physical_space.h>
#include <boost/optional.hpp>

#include <memory>
#include <map>

IGA_NAMESPACE_OPEN

/**
 * This namespace collects functions that work on the
 * physical spaces and fields or functions.
 * Such as:
 * - projections
 * - interpolations
 * - error computations
 */
namespace space_tools
{
/**
 * Determine the knot span index.
 *
 * @return The knot span index of the value @p u in the knot vector @p U.
 * @param[in] p Degree.
 * @param[in] u Knot values for which the span is requested.
 * @param[in] U Knot vector with repeated values.
 *
 * @note The implementation of this function is based on "The NURBS Book" Algorithm A2.1
 */
Index find_span(
    const int p,
    const Real u,
    const vector<Real> &U);


//TODO the order of parameters should be consistent
/**
 * Computes the norm distance between a generic function and
 * and an isogeometric type function.
 * More precisely,
 * res = | D^i(f-g) |_L2 with i=0 or i=1.
 * It also stores the local norm of each grid element in
 * the element error vector.
 *
 * @note mostly use to compute the convergence rates when the exact solution is
 *       known.
 * @todo document a little more
 */
template<class Space, LAPack la_pack = LAPack::trilinos>
Real integrate_difference(const typename Space::Func &exact_solution,
                          std::shared_ptr<const Space> space,
                          const Quadrature< Space::dim > &quad,
                          const Norm &norm_flag,
                          const Vector<la_pack> &solution_coefs,
                          vector<Real> &element_error);




// TODO (pauletti, Jun 18, 2014):use space::Function
// TODO (pauletti, Jun 18, 2014): use a quadrature table, and use a default quad if none provided
/**
 * Perform an (L2)-Projection the function @p func
 * onto the space @p space using the quadrature rule @p quad.
 *  The projection is a numerical vector (the coefficients of
 *  the projected function)
 */
template<class Space, LAPack la_pack = LAPack::trilinos>
Vector<la_pack>
projection_l2(const typename Space::Func &func,
              std::shared_ptr<const Space> space,
              const Quadrature<Space::dim> &quad);

/**
 * Projects (using the L2 scalar product) a function to the whole or part
 * of the boundary of the domain.
 * The piece of the domain is indicated by the boundary ids and the
 * projection is computed using the provided quadrature rule.
 *
 * The projected function is returned in boundary_values, a map containing all
 * indices of degrees of freedom at the boundary and the computed coefficient value
 * for this degree of freedom.
 *
 */
template<class Space, LAPack la_pack = LAPack::trilinos>
void project_boundary_values(
    const typename Space::Func &func,
    std::shared_ptr<const Space> space,
    const Quadrature<Space::dim-1> &quad,
    const std::set<boundary_id>  &boundary_ids,
    std::map<Index, Real>  &boundary_values);

/**
 * See documentation above.
 */
template<class Space, LAPack la_pack = LAPack::trilinos>
void project_boundary_values(
    const typename Space::Func &func,
    std::shared_ptr<const Space> space,
    const Quadrature<Space::dim-1> &quad,
    const boundary_id bdry_id,
    std::map<Index,Real>  &boundary_values) ;

//TODO: who uses the next function? delete
/**
 * Transform a set (CartesianProductArray) of points from the unit reference cube [0,1]^{dim} to the
 * element-based reference domain. It returns also the interval id of the elements that contains the points.
 * @param[in] reference_patch - Reference patch.
 * @param[in] points_ref - Set of points in the unit cube [0,1]^{dim}
 * @param[out] points_element - Set of point in the element-based reference domain.
 * @param knot_interval_id - ID of the intervals owning the points.
 */
template < int dim >
void reference_to_element(
    const CartesianGrid< dim > &reference_patch,
    const CartesianProductArray< Real, dim> &points_ref,
    CartesianProductArray< Real, dim> &points_element,
    CartesianProductArray< int, dim> &knot_interval_id) ;

} ;

IGA_NAMESPACE_CLOSE

#endif // #ifndef SPACE_TOOLS_H_
#endif

