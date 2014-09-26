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

// TODO: should be unified with dof tools?
//TODO Add doxygen description to the content of this file

#include <igatools/basis_functions/space_tools.h>
#include <igatools/basis_functions/physical_space_element_accessor.h>

#include <igatools/linear_algebra/dof_tools.h>
#include <igatools/linear_algebra/linear_solver.h>
#include <igatools/base/quadrature_lib.h>
#include <igatools/geometry/mapping_lib.h>
#include <igatools/geometry/mapping_slice.h>

#include <boost/numeric/ublas/io.hpp>


#include <chrono>

using std::shared_ptr;

using std::array;
using std::set;
using std::map;

using namespace std::chrono;


IGA_NAMESPACE_OPEN
#if 0

namespace
{
template <class RefSpace>
typename RefSpace::RefFaceSpace::MultiplicityTable
get_face_mult(std::shared_ptr<const RefSpace> ref_space, const Index face_id)
{
    const auto &active_dirs = UnitElement<RefSpace::dim>::face_active_directions[face_id];
    auto v_mult = ref_space->get_multiplicities();

    auto v_degree = ref_space->get_degree();
    typename RefSpace::RefFaceSpace::DegreeTable f_degree;
    for (int comp=0; comp<RefSpace::n_components; ++comp)
        for (int j=0; j<RefSpace::dim-1; ++j)
            f_degree(comp)[j] = v_degree(comp)[active_dirs[j]];

    typename RefSpace::RefFaceSpace::MultiplicityTable::parent_t f_mult;
    for (int comp=0; comp<RefSpace::n_components; ++comp)
        for (int j=0; j<RefSpace::dim-1; ++j)
            f_mult(comp).copy_data_direction(j, v_mult(comp).get_data_direction(active_dirs[j]));
    return typename RefSpace::RefFaceSpace::MultiplicityTable(f_mult, f_degree);
}



template <class RefSpace>
typename RefSpace::RefFaceSpace::DegreeTable
get_face_degree(std::shared_ptr<const RefSpace> ref_space, const Index face_id)
{
    const auto &active_dirs = UnitElement<RefSpace::dim>::face_active_directions[face_id];
    auto v_degree = ref_space->get_degree();
    typename RefSpace::RefFaceSpace::DegreeTable f_degree;
    for (int comp=0; comp<RefSpace::n_components; ++comp)
        for (int j=0; j<RefSpace::dim-1; ++j)
            f_degree(comp)[j] = v_degree(comp)[active_dirs[j]];
    return f_degree;
}



template <typename RefSpace>
EnableIf<!RefSpace::has_weights, std::shared_ptr<typename RefSpace::RefFaceSpace> >
create_face_ref_space(std::shared_ptr<const RefSpace> ref_space,
                      const Index face_id,
                      std::map<int,int> &elem_map)
{
    auto face_grid = ref_space->get_grid()->get_face_grid(face_id, elem_map);

    auto f_mult = get_face_mult<RefSpace> (ref_space, face_id);

    auto f_degree = get_face_degree<RefSpace> (ref_space, face_id);
    return RefSpace::RefFaceSpace::create(face_grid, f_mult);
}



/**
 * Given a Spline type reference space, create a spline type
 * space given by the trace of the function over the face.
 *
 * @param ref_space
 * @param face_id
 * @param elem_map
 * @return
 */
template <typename RefSpace>
EnableIf<RefSpace::has_weights, std::shared_ptr<typename RefSpace::RefFaceSpace> >
create_face_ref_space(std::shared_ptr<const RefSpace> ref_space,
                      const Index face_id,
                      std::map<int,int> &elem_map)
{
    auto face_grid = ref_space->get_grid()->get_face_grid(face_id, elem_map);

    auto f_mult = get_face_mult<RefSpace> (ref_space, face_id);
    auto f_degree = get_face_degree<RefSpace> (ref_space, face_id);

    //todo: auto f_weights = get_face_weight<RefSpace>(ref_space, face_id);
    //todo: return RefSpace::ref_face_space_t::create(face_grid, f_mult, f_degree, f_weights);
    return RefSpace::RefFaceSpace::create(face_grid, f_mult, f_degree);
}

}
#endif

namespace space_tools
{
Index find_span(
    const int p,
    const Real u,
    const vector<Real> &U)
{
    const Index m = U.size()-1;
    const Index n = m-p;

    // treat special case
    if (u == U[n+1])
        return n;

    // do binary search
    Index low = p;
    Index high = n+1;
    Index mid = (low+high)/2;
    while (u < U[mid] || u >= U[mid+1])
    {
        if (u < U[mid])
            high = mid;
        else
            low = mid;

        mid = (low+high)/2;
    }
    return mid;
}



template<class Space, LAPack la_pack>
Real
integrate_difference(const typename Space::Func &exact_solution,
                     std::shared_ptr<const Space> space,
                     const Quadrature< Space::dim > &quad,
                     const Norm &norm_flag,
                     const Vector<la_pack> &solution_coefs,
                     vector<Real> &element_error)
{
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

    ValueFlags flag = ValueFlags::point | ValueFlags::w_measure;

    if (is_L2_norm)
        flag |= ValueFlags::value;

    if (is_H1_seminorm)
        flag |= ValueFlags::gradient;

    const int n_points   =  quad.get_num_points();
    const int n_elements =  space->get_grid()->get_num_active_elems();

    Assert((element_error.size() == n_elements) || (element_error.size() == 0),
           ExcMessage("The size of the ouput vector is not correct."));
    if (element_error.size() == 0)
    {
        element_error.resize(n_elements);
    }

    using Value = typename Space::Func::Value;
    using Gradient = typename Space::Func::Gradient;

    ValueVector<Value>    u(n_points);
    ValueVector<Gradient> grad_u(n_points);

    Value err;
    Gradient grad_err;

    auto elem = space->begin();
    const auto end = space->end();
    // TODO (pauletti, Sep 12, 2014): fix next line
    Assert(true, ExcMessage(" fix next line "));
    //elem->init_cache(flag, quad);

    vector<Real>     norm_err_L2_square(n_elements);
    vector<Real> seminorm_err_H1_square(n_elements);

    for (; elem != end; ++elem)
    {
        // TODO (pauletti, Sep 12, 2014): fix next line
        Assert(true, ExcMessage(" fix next line "));
        // elem->fill_cache();
        const int elem_id = elem->get_flat_index();
        element_error[ elem_id ] = 0.0;

        const auto &map_at_points = elem->get_points();

        auto solution_coefs_elem =
            solution_coefs.get_local_coefs(elem->get_local_to_global());

        if (is_L2_norm)
        {
            const auto &uh = elem->evaluate_field(solution_coefs_elem);
            exact_solution.evaluate(map_at_points, u);

            Real element_err_L2_pow2 = 0.0;
            for (int iPt = 0; iPt < n_points; ++iPt)
            {
                err = uh[iPt] - u[iPt];
                element_err_L2_pow2 += err.norm_square() * elem->get_w_measures()[iPt];
            }
            element_error[ elem_id ] += element_err_L2_pow2;
        }


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

        element_error[ elem_id ] = sqrt(element_error[ elem_id ]);
    }

    Real err_pow2 = 0.0;
    for (const Real &elem_err : element_error)
        err_pow2 += elem_err * elem_err;

    Real total_error = sqrt(err_pow2);

    return total_error;
}



template<class Space, LAPack la_pack>
Vector<la_pack>
projection_l2(const typename Space::Func &func,
              shared_ptr<const Space> space,
              const Quadrature<Space::dim> &quad)
{
    const auto space_manager = space->get_space_manager();

//    const SparsityPattern sparsity_pattern(*space_manager);
    Matrix<la_pack> matrix(*space_manager);

    const auto space_dofs_set = space_manager->get_row_dofs();
    vector<Index> space_dofs(space_dofs_set.begin(),space_dofs_set.end());
    Vector<la_pack> rhs(space_dofs);
    Vector<la_pack> sol(space_dofs);



    //ValueFlags flag = ValueFlags::point | ValueFlags::value| ValueFlags::w_measure;
    const int n_qpoints = quad.get_num_points();

    ValueVector< typename Space::Point> eval_points(n_qpoints);
    ValueVector< typename Space::Value> func_at_eval_pts(n_qpoints);

    auto elem = space->begin();
    const auto elem_end = space->end();
    // TODO (pauletti, Sep 12, 2014): fix next line
    Assert(true, ExcMessage(" fix next line "));
    //elem->init_cache(flag, quad);
    const int n_basis = elem->get_num_basis();

    DenseVector local_rhs(n_basis);
    DenseMatrix local_matrix(n_basis,n_basis);

    for (; elem != elem_end; ++elem)
    {
        // TODO (pauletti, Sep 12, 2014): fix next line
        Assert(true, ExcMessage(" fix next line "));
        //elem->fill_cache();

        const auto eval_points = elem->get_points();
        func.evaluate(eval_points, func_at_eval_pts);

        local_matrix.clear();
        local_rhs.clear();

        const auto local_dofs = elem->get_local_to_global();

        // computing the upper triangular part of the local
        auto w_measures = elem->get_w_measures();
        for (int i = 0; i < n_basis; ++i)
        {
            const auto phi_i = elem->get_basis_values(i);

            for (int j = i; j < n_basis; ++j)
            {
                const auto phi_j = elem->get_basis_values(j);

                Real matrix_entry_ij = 0.0;
                for (int q = 0; q < n_qpoints; ++q)
                    matrix_entry_ij += scalar_product(phi_i[q], phi_j[q]) * w_measures[q];

                local_matrix(i,j) = matrix_entry_ij;
            }


            Real rhs_entry = 0.0;
            for (int q = 0; q < n_qpoints; q++)
                rhs_entry += scalar_product(func_at_eval_pts[q], phi_i[q]) * w_measures[q];

            local_rhs(i) = rhs_entry;
        }


        // copying the upper triangular part of the local matrix to the lower triangular part
        for (int i = 0; i < n_basis; ++i)
            for (int j = 0; j < i; ++j)
                local_matrix(i, j) = local_matrix(j, i);


        matrix.add_block(local_dofs,local_dofs,local_matrix);

        rhs.add_block(local_dofs,local_rhs);
    }
    matrix.fill_complete();

    const Real tolerance = 1.0e-15;
    const int max_num_iter = 1000;
    using LinSolver = LinearSolver<la_pack>;
    LinSolver solver(LinSolver::SolverType::CG,tolerance,max_num_iter);
    solver.solve(matrix, rhs, sol);

    return sol;
}




template<class Space, LAPack la_pack>
void
project_boundary_values(const typename Space::Func &func,
                        std::shared_ptr<const Space> space,
                        const Quadrature<Space::dim-1> &quad,
                        const std::set<boundary_id>  &boundary_ids,
                        std::map<Index, Real>  &boundary_values)
{
    const int n_faces = UnitElement<Space::dim>::faces_per_element;

    auto grid = space->get_grid();

    set<int> faces;
    auto bdry_begin = boundary_ids.begin();
    auto bdry_end   = boundary_ids.end();
    for (int face=0; face < n_faces; ++face)
    {
        const auto boundary_id = grid->get_boundary_id(face);
        if (find(bdry_begin, bdry_end, boundary_id) != bdry_end)
            faces.insert(face);
    }

    for (const Index &face_id : faces)
    {
        vector<Index> dof_map;
        auto face_space = space->get_face_space(face_id, dof_map);

        Vector<la_pack> proj_on_face =
            projection_l2<typename Space::FaceSpace,la_pack>(func, face_space, quad);

        const int face_n_dofs = dof_map.size();
        for (Index i = 0; i< face_n_dofs; ++i)
            boundary_values[dof_map[i]] = proj_on_face(i);
    }
}


template<class Space, LAPack la_pack>
void
project_boundary_values(const typename Space::Func &func,
                        std::shared_ptr<const Space> space,
                        const Quadrature<Space::dim-1> &quad,
                        const boundary_id bdry_id,
                        std::map<Index,Real>  &boundary_values)
{
    project_boundary_values<Space,la_pack>(
    func, space, quad,std::set<boundary_id>({{bdry_id}}),boundary_values);
}


#if 0
template < int dim >
void reference_to_element(
    const CartesianGrid< dim > &reference_patch,
    const CartesianProductArray< Real, dim > &points_ref,
    CartesianProductArray< Real, dim > &points_element,
    CartesianProductArray<   int, dim > &knot_interval_id)
{
    //----------------------------------------------------------------------------------------------
    for (int iDim = 0; iDim < dim; iDim++)
    {
        Assert(points_ref[ iDim ].size() == points_element[ iDim ].size(),
               ExcDimensionMismatch(points_ref[ iDim ].size(), points_element[ iDim ].size()));


        // get the point coordinates along the i-th direction
        const vector< Real > pt_coords = points_ref[ iDim ];


        const vector< Real > knot_coords = reference_patch.get_knot_coordinates(iDim);
        const Real knot_min = knot_coords.front();
        const Real knot_max = knot_coords.back();

        for (Real pt : pt_coords)
        {
            // check if the current point coordinate is contained in the reference patch
            AssertThrow(pt >= knot_min && pt <= knot_max, ExcMessage("An evaluation point is not in the current parametric domain."));

            // find the id of the knot interval for which the point coordinate belongs to
            int id;
            if (pt != knot_max)
            {
                // for the points that are not equal to the last knot
                // if u_{i} <= pt < u_{i+1} then the knot interval id is "i"
                id = upper_bound(knot_coords.begin(), knot_coords.end(), pt) - knot_coords.begin() - 1;
            }
            else
            {
                // for the points that are equal to the last knot
                // the last element_id is equal to the last knot_id-1
                // because always num_elements=num_knots-1 and the first id is zero
                id = knot_coords.size() - 2;
            }
            knot_interval_id[ iDim ].push_back(id);
        }

    }
    //----------------------------------------------------------------------------------------------

    array< int, dim > num_points_dim = points_ref.get_size();

    //----------------------------------------------------------------------------------------------
    // scale the points coordinates from the reference domain to the local element
    array< vector< Real >, dim > coords_scaled;
    for (int iDim = 0; iDim < dim; iDim++)
    {
        const vector< Real >   pt_coords = points_ref[ iDim ];
        const vector< Real > knot_coords = reference_patch.get_knot_coordinates(iDim);


        for (int iPt = 0; iPt < num_points_dim[ iDim ]; iPt++)
        {
            const int knot_id = knot_interval_id[ iDim ][ iPt ];

            const Real interval_size = knot_coords[ knot_id + 1 ] - knot_coords[ knot_id ];

            points_element[ iDim ][ iPt ] = (pt_coords[ iPt ] - knot_coords[ knot_id ]) / interval_size;
        }
    }
    //----------------------------------------------------------------------------------------------
}

#endif
};

IGA_NAMESPACE_CLOSE

#include <igatools/basis_functions/space_tools.inst>
