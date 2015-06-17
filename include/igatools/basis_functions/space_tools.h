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

#ifndef SPACE_TOOLS_H_
#define SPACE_TOOLS_H_

#include <igatools/functions/ig_function.h>
#include <igatools/geometry/grid_tools.h>
#include <igatools/functions/sub_function.h>

#include <igatools/basis_functions/physical_space_element.h>
#include <igatools/basis_functions/phys_space_element_handler.h>

#include <igatools/linear_algebra/epetra_solver.h>

#include<set>

IGA_NAMESPACE_OPEN
namespace space_tools
{


/**
 * Perform an (L2)-Projection the function @p func
 * onto the space @p space using the quadrature rule @p quad.
 *  The projection is a numerical vector (the coefficients of
 *  the projected function)
 */
template<class Space, LAPack la_pack = LAPack::trilinos_epetra>
std::shared_ptr<IgFunction<Space::dim,Space::codim,Space::range,Space::rank> >
projection_l2(const std::shared_ptr<const Function<Space::dim,Space::codim,Space::range,Space::rank>> function,
              std::shared_ptr<const Space> space,
              const Quadrature<Space::dim> &quad,
              const std::string &dofs_property = DofProperties::active)
{
    using ProjFunc = IgFunction<Space::dim,Space::codim,Space::range,Space::rank>;
//    std::shared_ptr<ProjFunc> projection;

    Epetra_SerialComm comm;

    auto map = EpetraTools::create_map(*space, "active", comm);
    auto graph = EpetraTools::create_graph(*space,"active",*space,"active",*map,*map);

    auto matrix = EpetraTools::create_matrix(graph);
    auto rhs = EpetraTools::create_vector(map);
    auto sol = EpetraTools::create_vector(map);

    const auto space_grid =    space->get_grid();
    const auto func_grid = function->get_grid();

    if (space_grid == func_grid)
    {
        auto func = function->clone();
        const int dim = Space::dim;

        auto func_flag = ValueFlags::point | ValueFlags::value;
        func->reset(func_flag, quad);

//        using ElementHandler = typename Space::ElementHandler;
//        auto sp_filler = ElementHandler::create(space);
        auto sp_filler = space->get_elem_handler();
        auto sp_flag = ValueFlags::point | ValueFlags::value |
                       ValueFlags::w_measure;
        sp_filler->reset(sp_flag, quad);

        auto f_elem = func->begin();
        auto elem = space->begin();
        auto end  = space->end();

        func->init_element_cache(f_elem);
        sp_filler->init_element_cache(elem);

        const int n_qp = quad.get_num_points();

        for (; elem != end; ++elem, ++f_elem)
        {
            const int n_basis = elem->get_num_basis(dofs_property);
            DenseVector loc_rhs(n_basis);
            DenseMatrix loc_mat(n_basis, n_basis);

            func->fill_element_cache(f_elem);
            sp_filler->fill_element_cache(elem);

            loc_mat = 0.;
            loc_rhs = 0.;

            auto f_at_qp = f_elem->template get_values<_Value,dim>(0);
            auto phi = elem->template get_basis_element<_Value>(dofs_property);

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

            const auto elem_dofs = elem->get_local_to_global(dofs_property);
            matrix->add_block(elem_dofs,elem_dofs,loc_mat);
            rhs->add_block(elem_dofs,loc_rhs);
        }
        matrix->FillComplete();
    }
    else
    {
        Assert(space_grid->same_knots_or_refinement_of(*func_grid),
               ExcMessage("The space grid is not a refinement of the function grid."));

//        Assert(false,ExcNotImplemented());
        auto func = function->clone();
        const int dim = Space::dim;


        auto func_flag = ValueFlags::point | ValueFlags::value;
        func->reset(func_flag, quad);

        auto sp_filler = space->get_elem_handler();
        auto sp_flag = ValueFlags::point | ValueFlags::value |
                       ValueFlags::w_measure;
        sp_filler->reset(sp_flag, quad);

        auto f_elem = func->begin();
        auto elem = space->begin();
        auto end  = space->end();

        auto map_elems_id_fine_coarse =
            grid_tools::build_map_elements_id_between_cartesian_grids(
                *space->get_grid(),*func->get_grid());

        func->init_element_cache(f_elem);
        sp_filler->init_element_cache(elem);

        const int n_qp = quad.get_num_points();

        for (const auto &elems_id_pair : map_elems_id_fine_coarse)
        {
            elem->move_to(elems_id_pair.first);
            f_elem->move_to(elems_id_pair.second);

            const int n_basis = elem->get_num_basis(dofs_property);
            DenseVector loc_rhs(n_basis);
            DenseMatrix loc_mat(n_basis, n_basis);

            func->fill_element_cache(f_elem);
            sp_filler->fill_element_cache(elem);

            loc_mat = 0.;
            loc_rhs = 0.;


            //---------------------------------------------------------------------------
            // the function is supposed to be defined on the same grid of the space or coarser
            const auto &elem_grid_accessor = elem->as_cartesian_grid_element_accessor();
            auto quad_in_func_elem = quad;
            quad_in_func_elem.dilate_translate(
                elem_grid_accessor.
                template get_coordinate_lengths<dim>(0),
                elem_grid_accessor.vertex(0));

            auto one_div_f_elem_size = f_elem->template get_coordinate_lengths<dim>(0);
            for (int dir : UnitElement<dim>::active_directions)
                one_div_f_elem_size[dir] = 1.0/one_div_f_elem_size[dir];

            auto f_elem_vertex = -f_elem->vertex(0);
            quad_in_func_elem.translate(f_elem_vertex);
            quad_in_func_elem.dilate(one_div_f_elem_size);


            auto f_at_qp = f_elem->template evaluate_at_points<_Value>(quad_in_func_elem);
            //---------------------------------------------------------------------------


            auto phi = elem->template get_basis_element<_Value>(dofs_property);

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

            const auto elem_dofs = elem->get_local_to_global(dofs_property);
            matrix->add_block(elem_dofs,elem_dofs,loc_mat);
            rhs->add_block(elem_dofs,loc_rhs);
        }
        matrix->FillComplete();
    }

    auto solver = EpetraTools::create_solver(matrix, sol, rhs);
    auto result = solver->solve();
    AssertThrow(result == Belos::ReturnType::Converged,
                ExcMessage("No convergence."));

    return ProjFunc::create(std::const_pointer_cast<Space>(space),
                            sol, dofs_property);
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
template<class Space>
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

        auto proj = projection_l2<SubSpace>(sub_func, sub_space, quad);

        const auto &coef = proj->get_coefficients();
        const int face_n_dofs = dof_map.size();
        for (Index i = 0; i< face_n_dofs; ++i)
            boundary_values[dof_map[i]] = coef[i];
    }
}



/**
 * Returns the list of global ids of the non zero basis functions
 * on the faces with the given boundary ids.
 */
template<class Space>
std::set<Index>
get_boundary_dofs(std::shared_ptr<const Space> space,
                  const std::set<boundary_id>  &boundary_ids)
{
    const int dim   = Space::dim;
    std::set<Index> dofs;
    const int sub_dim = dim - 1;

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
        auto s_dofs = space->template get_boundary_dofs<sub_dim>(s_id);
        dofs.insert(s_dofs.begin(), s_dofs.end());
    }

    return dofs;
}



// TODO (pauletti, Mar 18, 2015): this could be given a more general use
static const SafeSTLArray<ValueFlags, 3> order_to_flag =
{ValueFlags::value,ValueFlags::gradient,ValueFlags::hessian};

template<int order, int dim, int codim = 0, int range = 1, int rank = 1>
void norm_difference(Function<dim, codim, range, rank> &f,
                     Function<dim, codim, range, rank> &g,
                     const Quadrature<dim> &quad,
                     const Real p,
                     SafeSTLVector<Real> &element_error)
{
    const bool is_inf = p==std::numeric_limits<Real>::infinity()? true : false;
    auto flag = ValueFlags::point | ValueFlags::w_measure | order_to_flag[order];

    f.reset(flag, quad);
    g.reset(flag, quad);
    const int n_points = quad.get_num_points();

    auto elem_f = f.begin();
    auto elem_g = g.begin();
    auto end = f.end();

    const auto topology = Topology<dim>();

    f.init_cache(elem_f, topology);
    g.init_cache(elem_g, topology);

    for (; elem_f != end; ++elem_f, ++elem_g)
    {
        f.fill_cache(elem_f, topology, 0);
        g.fill_cache(elem_g, topology, 0);

        const int elem_id = elem_f->get_flat_index();

        auto f_val = elem_f->template get_values<ValueType<order>,dim>(0);
        auto g_val = elem_g->template get_values<ValueType<order>,dim>(0);
        auto w_meas = elem_f->template get_w_measures<dim>(0);

        Real elem_diff_pow_p = 0.0;
        Real val;
        for (int iPt = 0; iPt < n_points; ++iPt)
        {
            const auto err = f_val[iPt] - g_val[iPt];
            val = err.norm_square();
            if (is_inf)
                elem_diff_pow_p = std::max(elem_diff_pow_p, fabs(sqrt(val)));
            else
                elem_diff_pow_p += std::pow(val,p/2.) * w_meas[iPt];
        }
        element_error[ elem_id ] += elem_diff_pow_p;
    }
}



template<int dim, int codim = 0, int range = 1, int rank = 1>
Real l2_norm_difference(Function<dim, codim, range, rank> &f,
                        Function<dim, codim, range, rank> &g,
                        const Quadrature<dim> &quad,
                        SafeSTLVector<Real> &elem_error)
{
    const Real p=2.;
    const Real one_p = 1./p;
    const int order=0;

    space_tools::norm_difference<order,dim, codim, range, rank>(f, g, quad, p, elem_error);

    Real err = 0;
    for (Real &loc_err : elem_error)
    {
        err += loc_err;
        loc_err = std::pow(loc_err,one_p);

    }

    return std::pow(err,one_p);
}



template<int dim, int codim = 0, int range = 1, int rank = 1>
Real h1_norm_difference(Function<dim, codim, range, rank> &f,
                        Function<dim, codim, range, rank> &g,
                        const Quadrature<dim> &quad,
                        SafeSTLVector<Real> &elem_error)
{
    const Real p=2.;
    const Real one_p = 1./p;

    space_tools::norm_difference<0,dim, codim, range, rank>(f, g, quad, p, elem_error);
    space_tools::norm_difference<1,dim, codim, range, rank>(f, g, quad, p, elem_error);

    Real err = 0;
    for (Real &loc_err : elem_error)
    {
        err += loc_err;
        loc_err = std::pow(loc_err,one_p);
    }

    return std::pow(err,one_p);
}



template<int dim, int codim = 0, int range = 1, int rank = 1>
Real inf_norm_difference(Function<dim, codim, range, rank> &f,
                         Function<dim, codim, range, rank> &g,
                         const Quadrature<dim> &quad,
                         SafeSTLVector<Real> &elem_error)
{
    const Real p=std::numeric_limits<Real>::infinity();
    space_tools::norm_difference<0, dim, codim, range, rank>(f, g, quad, p, elem_error);
    Real err = 0;
    for (const Real &loc_err : elem_error)
        err = std::max(err,loc_err);

    return err;
}

};

IGA_NAMESPACE_CLOSE

#endif
