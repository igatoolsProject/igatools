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
#include <igatools/basis_functions/multiplicity.h>
#include <igatools/basis_functions/physical_space_element_accessor.h>

#include <igatools/linear_algebra/dof_tools.h>
#include <igatools/linear_algebra/linear_solver.h>
#include <igatools/base/quadrature_lib.h>
#include <igatools/geometry/mapping_lib.h>
#include <igatools/geometry/mapping_slice.h>

#include <boost/numeric/ublas/io.hpp>

using std::shared_ptr;
using std::vector;
using std::array;
using std::set;
using std::map;

IGA_NAMESPACE_OPEN


namespace
{
template <class RefSpace>
StaticMultiArray< Multiplicity< RefSpace::RefFaceSpace::dim >,
                  RefSpace::RefFaceSpace::range, RefSpace::RefFaceSpace::rank>
                  get_face_mult(std::shared_ptr<const RefSpace> ref_space, const Index face_id)
{
    const auto &active_dirs = UnitElement<RefSpace::dim>::face_active_directions[face_id];
    auto v_mult = ref_space->get_multiplicities();
    StaticMultiArray< Multiplicity < RefSpace::RefFaceSpace::dim >,
                      RefSpace::RefFaceSpace::range, RefSpace::RefFaceSpace::rank>  f_mult;
    for (int comp=0; comp<RefSpace::n_components; ++comp)
        for (int j=0; j<RefSpace::dim-1; ++j)
            f_mult(comp).copy_data_direction(j, v_mult(comp).get_data_direction(active_dirs[j]));
    return f_mult;
}


template <class RefSpace>
StaticMultiArray<TensorIndex<RefSpace::dim-1>, RefSpace::range, RefSpace::rank>
get_face_degree(std::shared_ptr<const RefSpace> ref_space, const Index face_id)
{
    const auto &active_dirs = UnitElement<RefSpace::dim>::face_active_directions[face_id];
    auto v_degree = ref_space->get_degree();
    StaticMultiArray<TensorIndex<RefSpace::dim-1>, RefSpace::range, RefSpace::rank>  f_degree;
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
    return RefSpace::RefFaceSpace::create(face_grid, f_mult, f_degree);
}


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

namespace space_tools
{

Index find_span(
    const int p,
    const Real u,
    const std::vector<Real> &U)
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


template <class Space>
std::shared_ptr< FaceSpace<Space> >
get_face_space(std::shared_ptr<const Space> space,
               const Index face_id,
               std::vector<Index> &face_to_element_dofs)
{
    using RefSpace = typename Space::RefSpace;

    using FPF = typename Space::PushForwardType::FacePushForward;
    auto ref_space = space->get_reference_space();


    auto elem_map = std::make_shared<std::map<int,int> >();

    auto face_ref_sp = create_face_ref_space<RefSpace>(ref_space, face_id, *elem_map);

    auto pf   = space->get_push_forward();
    auto map  = pf->get_mapping();
    auto fmap = MappingSlice<Space::dim-1, Space::codim + 1>::
                create(map, face_id, face_ref_sp->get_grid(), elem_map);

    auto fpf = FPF::create(fmap);

    auto face_space = FaceSpace<Space>::create(face_ref_sp,fpf);

    const auto &active_dirs = UnitElement<Space::dim>::face_active_directions[face_id];
    const auto const_dir = UnitElement<Space::dim>::face_constant_direction[face_id];
    const auto face_side = UnitElement<Space::dim>::face_side[face_id];

    TensorIndex<Space::dim> tensor_index;

    face_to_element_dofs.resize(face_ref_sp->get_num_basis());
    int k=0;
    int offset=0;
    for (int comp = 0; comp < Space::RefFaceSpace::n_components; ++comp)
    {
        const int face_n_basis = face_ref_sp->get_num_basis(comp);
        for (Index i = 0; i < face_n_basis; ++i, ++k)
        {
            const auto f_tensor_idx = face_ref_sp->flat_to_tensor(i,comp);
            const int fixed_idx =
                face_side * (ref_space->get_num_basis(comp,const_dir) - 1);
            for (int j = 0; j < Space::dim-1; ++j)
                tensor_index[active_dirs[j]] =  f_tensor_idx[j];
            tensor_index[const_dir] = fixed_idx;

            const Index dof = ref_space->tensor_to_flat(tensor_index, comp);

            face_to_element_dofs[k] = offset + dof;
        }
        offset += ref_space->get_num_basis(comp);
    }

    return face_space;
}



template<class Space>
Real integrate_difference(std::shared_ptr<const Func<Space> > exact_solution,
                          std::shared_ptr<const Space> space,
                          const Quadrature< Space::dim > &quad,
                          const Norm &norm_flag,
                          const Vector &solution_coefs,
                          std::vector< Real > &element_error)
{
    bool is_L2_norm     = contains(norm_flag, Norm::L2) ;
    bool is_H1_norm     = contains(norm_flag, Norm::H1) ;
    bool is_H1_seminorm = contains(norm_flag, Norm::H1_semi) ;

    Assert(is_L2_norm || is_H1_seminorm || is_H1_norm,
           ExcMessage("No active flag for the error norm.")) ;


    Assert(!((is_L2_norm && is_H1_seminorm) ||
             (is_L2_norm && is_H1_norm) ||
             (is_H1_seminorm && is_H1_norm)),
           ExcMessage("Only a single flag for the error norm can be used.")) ;


    if (is_H1_norm)
    {
        is_L2_norm     = true ;
        is_H1_seminorm = true ;
    }

    ValueFlags flag = ValueFlags::point | ValueFlags::w_measure ;

    if (is_L2_norm)
        flag |= ValueFlags::value ;

    if (is_H1_seminorm)
        flag |= ValueFlags::gradient ;



    const int n_points   =  quad.get_num_points();
    const int n_elements =  space->get_grid()->get_num_elements() ;


    Assert((element_error.size() == n_elements) || (element_error.size() == 0),
           ExcMessage("The size of the ouput vector is not correct.")) ;
    if (element_error.size() == 0)
    {
        element_error.resize(n_elements) ;
    }


    typedef typename Func<Space>::ValueType ValuePhys_t ;
    typedef typename Func<Space>::GradientType GradientPhys_t;

    vector< ValuePhys_t > u(n_points) ;
    vector< GradientPhys_t > grad_u(n_points) ;

    ValuePhys_t err ;
    GradientPhys_t grad_err ;

    auto elem = space->begin() ;
    const auto end = space->end() ;
    elem->init_values(flag, quad);

    vector< Real >     norm_err_L2_square(n_elements) ;
    vector< Real > seminorm_err_H1_square(n_elements) ;

    for (; elem != end ; ++elem)
    {
        elem->fill_values();
        const int elem_id = elem->get_flat_index();
        element_error[ elem_id ] = 0.0 ;

        const auto &map_at_points = elem->get_points() ;

        vector<Real> solution_coefs_elem = dof_tools::get_local_coefs(
                                               solution_coefs,elem->get_local_to_global());

        if (is_L2_norm)
        {
            const auto &uh = elem->evaluate_field(solution_coefs_elem);
            exact_solution->evaluate(map_at_points, u);

            Real element_err_L2_pow2 = 0.0 ;
            for (int iPt = 0 ; iPt < n_points ; ++iPt)
            {
                err = uh[iPt] - u[iPt] ;
                element_err_L2_pow2 += err.norm_square() * elem->get_w_measures()[iPt];
            }
            element_error[ elem_id ] += element_err_L2_pow2;
        }


        if (is_H1_seminorm)
        {
            const auto &grad_uh = elem->evaluate_field_gradients(solution_coefs_elem) ;
            exact_solution->evaluate_gradients(map_at_points, grad_u) ;

            Real element_err_semiH1_pow2 = 0.0 ;
            for (int iPt = 0 ; iPt < n_points ; ++iPt)
            {
                grad_err = grad_uh[iPt] - grad_u[iPt] ;

                element_err_semiH1_pow2 += grad_err.norm_square() * elem->get_w_measures()[iPt] ;
            }
            element_error[ elem_id ] = element_err_semiH1_pow2 ;
        }

        element_error[ elem_id ] = sqrt(element_error[ elem_id ]) ;
    }

    Real err_pow2 = 0.0 ;
    for (const Real & elem_err : element_error)
        err_pow2 += elem_err * elem_err ;

    Real total_error = sqrt(err_pow2) ;

    return (total_error) ;
}





template<class Space>
Vector projection_l2(const Function<Space::space_dim,Space::range,Space::rank> &func,
                     shared_ptr<const Space> space,
                     const Quadrature<Space::dim> &quad)
{
    static const int space_dim = Space::space_dim;
    static const int range = Space::range;
    static const int rank = Space::rank;

    const auto sparsity_pattern = dof_tools::get_sparsity_pattern(space) ;
    Matrix matrix(sparsity_pattern);

    const auto space_dofs = sparsity_pattern.get_row_dofs() ;
    Vector rhs(space_dofs) ;
    Vector sol(space_dofs) ;



    ValueFlags flag = ValueFlags::point | ValueFlags::value| ValueFlags::w_measure ;
    const int n_qpoints = quad.get_num_points();

    vector< Point<space_dim> > eval_points(n_qpoints);
    vector< typename Function<space_dim,range,rank>::ValueType > func_at_eval_pts(n_qpoints);

    auto elem = space->begin() ;
    const auto elem_end = space->end() ;
    elem->init_values(flag, quad);
    const int n_basis = elem->get_num_basis();
    DenseVector local_rhs(n_basis);
    DenseMatrix local_matrix(n_basis,n_basis);

    for (; elem != elem_end ; ++elem)
    {
        elem->fill_values();

        const auto eval_points = elem->get_points() ;
        func.evaluate(eval_points, func_at_eval_pts) ;

        local_matrix.clear();
        local_rhs.clear();

        const auto local_dofs = elem->get_local_to_global();

        // computing the upper triangular part of the local
        auto w_measures = elem->get_w_measures();
        for (int i = 0; i < n_basis; ++i)
        {
            const auto phi_i = elem->get_basis_values(i) ;

            for (int j = i ; j < n_basis ; ++j)
            {
                const auto phi_j = elem->get_basis_values(j) ;

                Real matrix_entry_ij = 0.0 ;
                for (int q = 0; q < n_qpoints; ++q)
                    matrix_entry_ij += scalar_product(phi_i[q], phi_j[q]) * w_measures[q] ;

                local_matrix(i,j) = matrix_entry_ij ;
            }


            Real rhs_entry = 0.0 ;
            for (int q = 0; q < n_qpoints; q++)
                rhs_entry += scalar_product(func_at_eval_pts[q], phi_i[q]) * w_measures[q];

            local_rhs(i) = rhs_entry ;
        }


        // copying the upper triangular part of the local matrix to the lower triangular part
        for (int i = 0 ; i < n_basis ; ++i)
            for (int j = 0 ; j < i ; ++j)
                local_matrix(i, j) = local_matrix(j, i) ;


        matrix.add_block(local_dofs,local_dofs,local_matrix) ;

        rhs.add_block(local_dofs,local_rhs) ;
    }
    matrix.fill_complete();

    const Real tolerance = 1.0e-15;
    const int max_num_iter = 1000;
    using lin_solver_t = LinearSolver<LinearAlgebraPackage::trilinos>;
    lin_solver_t solver(lin_solver_t::Type::CG,tolerance,max_num_iter) ;
    solver.solve(matrix, rhs, sol);

    return sol;
}





template<class Space>
void
project_boundary_values(const Function<Space::space_dim,Space::range,Space::rank> &func,
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

    for (const Index& face_id : faces)
    {
        vector<Index> dof_map;
        auto face_space = get_face_space(space, face_id, dof_map);

        Vector proj_on_face =
            projection_l2<FaceSpace<Space> >(func, face_space, quad);

        const int face_n_dofs = dof_map.size() ;
        for (Index i = 0 ; i< face_n_dofs ; ++i)
            boundary_values[dof_map[i]] = proj_on_face(i);
    }
}


template<class Space>
void
project_boundary_values(const Func<Space> &func,
                        std::shared_ptr<const Space> space,
                        const Quadrature<Space::dim-1> &quad,
                        const boundary_id bdry_id,
                        std::map<Index,Real>  &boundary_values)
{
    project_boundary_values(func, space, quad,
    std::set<boundary_id>({{bdry_id}}),
    boundary_values);
}



template < int dim >
void reference_to_element(
    const CartesianGrid< dim > &reference_patch,
    const CartesianProductArray< Real, dim > &points_ref,
    CartesianProductArray< Real, dim > &points_element,
    CartesianProductArray<   int, dim > &knot_interval_id)
{
    //----------------------------------------------------------------------------------------------
    for (int iDim = 0 ; iDim < dim ; iDim++)
    {
        Assert(points_ref[ iDim ].size() == points_element[ iDim ].size(),
               ExcDimensionMismatch(points_ref[ iDim ].size(), points_element[ iDim ].size())) ;


        // get the point coordinates along the i-th direction
        const vector< Real > pt_coords = points_ref[ iDim ] ;


        const vector< Real > knot_coords = reference_patch.get_knot_coordinates(iDim) ;
        const Real knot_min = knot_coords.front() ;
        const Real knot_max = knot_coords.back() ;

        for (Real pt : pt_coords)
        {
            // check if the current point coordinate is contained in the reference patch
            AssertThrow(pt >= knot_min && pt <= knot_max, ExcMessage("An evaluation point is not in the current parametric domain.")) ;

            // find the id of the knot interval for which the point coordinate belongs to
            int id ;
            if (pt != knot_max)
            {
                // for the points that are not equal to the last knot
                // if u_{i} <= pt < u_{i+1} then the knot interval id is "i"
                id = upper_bound(knot_coords.begin(), knot_coords.end(), pt) - knot_coords.begin() - 1 ;
            }
            else
            {
                // for the points that are equal to the last knot
                // the last element_id is equal to the last knot_id-1
                // because always num_elements=num_knots-1 and the first id is zero
                id = knot_coords.size() - 2 ;
            }
            knot_interval_id[ iDim ].push_back(id) ;
        }

    }
    //----------------------------------------------------------------------------------------------

    array< int, dim > num_points_dim = points_ref.get_size() ;

    //----------------------------------------------------------------------------------------------
    // scale the points coordinates from the reference domain to the local element
    array< vector< Real >, dim > coords_scaled ;
    for (int iDim = 0 ; iDim < dim ; iDim++)
    {
        const vector< Real >   pt_coords = points_ref[ iDim ] ;
        const vector< Real > knot_coords = reference_patch.get_knot_coordinates(iDim) ;


        for (int iPt = 0 ; iPt < num_points_dim[ iDim ] ; iPt++)
        {
            const int knot_id = knot_interval_id[ iDim ][ iPt ] ;

            const Real interval_size = knot_coords[ knot_id + 1 ] - knot_coords[ knot_id ] ;

            points_element[ iDim ][ iPt ] = (pt_coords[ iPt ] - knot_coords[ knot_id ]) / interval_size ;
        }
    }
    //----------------------------------------------------------------------------------------------
}


};

IGA_NAMESPACE_CLOSE

#include <igatools/basis_functions/space_tools.inst>
