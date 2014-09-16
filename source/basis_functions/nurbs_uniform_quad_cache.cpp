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

#include <igatools/basis_functions/nurbs_uniform_quad_cache.h>
#include <igatools/utils/multi_array_utils.h>

#include <algorithm>
using std::shared_ptr;

IGA_NAMESPACE_OPEN

namespace
{
ValueFlags
nurbs_to_bspline_flags(const ValueFlags nurbs_flags)
{
    int max_der_order = -1;
    ValueFlags bspline_flags = nurbs_flags;
    if (contains(nurbs_flags, ValueFlags::value))
    {
        max_der_order=std::max(max_der_order,0);
        bspline_flags |= ValueFlags::value;
    }

    if (contains(nurbs_flags, ValueFlags::face_value))
    {
        max_der_order=std::max(max_der_order,0);
        bspline_flags |= ValueFlags::face_value;
    }

    if (contains(nurbs_flags, ValueFlags::gradient))
    {
        max_der_order=std::max(max_der_order,1);
        bspline_flags |= ValueFlags::value |
                         ValueFlags::gradient;
    }

    if (contains(nurbs_flags, ValueFlags::face_gradient))
    {
        max_der_order=std::max(max_der_order,1);
        bspline_flags |= ValueFlags::face_value |
                         ValueFlags::face_gradient;
    }

    if (contains(nurbs_flags, ValueFlags::hessian))
    {
        max_der_order=std::max(max_der_order,2);
        bspline_flags |= ValueFlags::value |
                         ValueFlags::gradient |
                         ValueFlags::hessian;
    }

    if (contains(nurbs_flags, ValueFlags::face_hessian))
    {
        max_der_order=std::max(max_der_order,2);
        bspline_flags |= ValueFlags::face_value |
                         ValueFlags::face_gradient |
                         ValueFlags::face_hessian;
    }
    Assert(max_der_order>=0, ExcMessage("Not a right ValueFlag"));

    return bspline_flags;
}

};



template<int dim_, int range_ , int rank_>
NURBSUniformQuadCache<dim_, range_, rank_>::
NURBSUniformQuadCache(shared_ptr<const Space> space,
                      const ValueFlags flag,
                      const Quadrature<dim> &quad)
    :
    GridUniformQuadCache<dim_>(space->get_grid(), flag, quad),
    space_(space),
//    n_basis_(space_->get_num_all_element_basis()),
    flags_(flag),
    face_flags_(flag),
    quad_(quad),
    bspline_uniform_quad_cache_(
        space->get_spline_space(),
        nurbs_to_bspline_flags(flag),
        quad)
{
    Assert(contains(flag, ValueFlags::none),
           ExcMessage("Nothing to reset"));
}

#if 0
template<int dim_, int range_ , int rank_>
BSplineUniformQuadCache<dim_, range_, rank_>::
BSplineUniformQuadCache(shared_ptr<const Space> space,
                        const ValueFlags flag,
                        const Quadrature<dim> &quad)
    :
    GridUniformQuadCache<dim_>(space->get_grid(), flag, quad),
    space_(space),
    n_basis_(space_->get_num_all_element_basis()),
    flags_(flag),
    face_flags_(flag),
    quad_(quad),
    splines1d_(space->get_grid()->get_num_intervals(),
               BasisValues(space->get_components_map()))
{


    // Compute the component offsets
    comp_offset_[0] = 0;
    for (int j = 1; j < Space::n_components; ++j)
        comp_offset_[j] = comp_offset_[j-1] + n_basis_.comp_dimension[j-1];


    const auto &n_inter = space->get_grid()->get_num_intervals();
    const auto &n_points = quad.get_num_points_direction();

    // Allocate space for the BasisValues1D
    for (int dir = 0 ; dir < dim ; ++dir)
    {
        const auto &n_pts = n_points[dir];
        for (int j = 0 ; j < n_inter[dir] ; ++j)
        {
            auto &splines1d = splines1d_.entry(dir, j);
            for (auto comp : splines1d.get_active_components_id())
                splines1d[comp].resize(n_derivatives, n_basis_[comp][dir], n_pts);
        }
    }

    /*
     * For each direction, interval and component we compute the 1D bspline
     * basis evaluate at the 1D component of the tensor product quadrature
     */
    const auto &degree      = space->get_degree();
    const auto &bezier_op   = space_->operators_;
    const auto &points      = quad_.get_points();
    const auto &lengths = this->lengths_;

    BasisValues bernstein_values(n_basis_.get_comp_map());

    for (int dir = 0 ; dir < dim ; ++dir)
    {
        // fill values and derivatives of the Bernstein's polynomials at
        // quad points in [0,1]
        for (auto comp : bernstein_values.get_active_components_id())
        {
            const int deg = degree[comp][dir];
            bernstein_values[comp].resize(n_derivatives, deg+1, n_points[dir]);
            const auto &pt_coords = points.get_data_direction(dir);
            for (int order = 0; order < n_derivatives; ++order)
                bernstein_values[comp].get_derivative(order) =
                    BernsteinBasis::derivative(order, deg, pt_coords);
        }

        const auto &inter_lengths = lengths.get_data_direction(dir);
        for (int j = 0 ; j < n_inter[dir] ; ++j)
        {
            auto &splines1d = splines1d_.entry(dir, j);
            for (auto comp : splines1d.get_active_components_id())
            {
                const auto &berns_values = bernstein_values[comp];
                auto &basis = splines1d[comp];
                const auto &oper = bezier_op.get_operator(comp,dir)[j];
                const Real one_div_size = 1.0 / inter_lengths[j];
                for (int order = 0; order < n_derivatives; ++order)
                {
                    const Real scale = std::pow(one_div_size, order);
                    const auto &b_values = berns_values.get_derivative(order);
                    basis.get_derivative(order) =
                        scale * prec_prod(oper, b_values);
                }
            }
        }

    }

//        for (auto comp : splines1d_.get_active_components_id())
//    {
//        auto &splines1d = splines1d_[comp];
//
//        {
//            const int num_intervals = n_intervals[jDim];
//            const int deg = degree[comp][jDim];
//            BasisValues1d bernstein_values(n_derivatives, deg+1, n_points[jDim]);
//
//            const auto &pt_coords = points.get_data_direction(jDim);
//
//            // fill values and derivatives of the Bernstein's polynomials at
//            // quad points in [0,1]
//            for (int order = 0; order < n_derivatives; ++order)
//                bernstein_values.get_derivative(order) =
//                    BernsteinBasis::derivative(order, deg, pt_coords);
//
//            const auto &bez_iComp_jDim = bezier_op.get_operator(comp,jDim);
//            const auto &lengths_jDim = lengths.get_data_direction(jDim);
//
//            // compute the one dimensional B-splines at quad point on the reference interval
//            for (int i = 0 ; i < num_intervals ; ++i)
//            {
//                const auto &M = bez_iComp_jDim[i];
//                const Real one_div_size = 1.0 / lengths_jDim[i];
//                BasisValues1d &basis = splines1d.entry(jDim,i);
//
//                for (int order = 0; order < n_derivatives; ++order)
//                {
//                    const Real scaling_coef = std::pow(one_div_size, order);
//                    basis.get_derivative(order) = scaling_coef * prec_prod(M, bernstein_values.get_derivative(order));
//                } //end loop order
//
//            } // end loop interval
//        } // end loop jDim
//    } // end loop comp
}
#endif


template<int dim_, int range_ , int rank_>
void
NURBSUniformQuadCache<dim_, range_, rank_>::
init_element_cache(ElementAccessor &elem)
{
    base_t::init_element_cache(elem);

    // init the element values for the cache of the BSplineElementAccessor
    bspline_uniform_quad_cache_.init_element_cache(elem.bspline_element_accessor_);


    //-----------------------------------------------------------------------------
    auto &cache = elem.local_cache_;
    if (cache == nullptr)
    {
        using Cache = typename ElementAccessor::LocalCache;
        cache = shared_ptr<Cache>(new Cache);
    }


    auto n_basis = space_->get_num_all_element_basis();
    auto &elem_cache = cache->elem_values_;
    elem_cache.resize(flags_, quad_, n_basis);

    auto &face_cache = cache->face_values_;
    for (auto f: base_t::faces)
        face_cache[f].resize(face_flags_, quad_, n_basis, f);
    //-----------------------------------------------------------------------------
}



template<int dim_, int range_ , int rank_>
void
NURBSUniformQuadCache<dim_, range_, rank_>::
init_element_cache(ElementIterator &elem)
{
    init_element_cache(elem.get_accessor());
}

#if 0
template <int dim, int range, int rank>
void
BSplineUniformQuadCache<dim, range, rank>::
copy_to_inactive_components_values(const vector<Index> &inactive_comp,
                                   const std::array<Index, n_components> &active_map,
                                   ValueTable<Value> &D_phi) const
{
    const Size num_points = quad_.get_num_points_direction().flat_size();
    for (int comp : inactive_comp)
    {
        const auto act_comp = active_map[comp];
        const auto n_basis = n_basis_.comp_dimension[comp];
        const Size act_offset = comp_offset_[act_comp];
        const Size offset     = comp_offset_[comp];
        for (Size basis_i = 0; basis_i < n_basis;  ++basis_i)
        {
            const auto act_D_phi = D_phi.get_function_view(act_offset+basis_i);
            auto inact_D_phi = D_phi.get_function_view(offset+basis_i);
            for (int qp = 0; qp < num_points; ++qp)
                inact_D_phi[qp](comp) = act_D_phi[qp](act_comp);
        }
    }
}



template <int dim, int range, int rank>
template <int order>
void
BSplineUniformQuadCache<dim, range, rank>::
copy_to_inactive_components(const vector<Index> &inactive_comp,
                            const std::array<Index, n_components> &active_map,
                            ValueTable<Derivative<order>> &D_phi) const
{
    const Size num_points = quad_.get_num_points_direction().flat_size();
    const Size n_ders = Derivative<order>::size;
    for (int comp : inactive_comp)
    {
        const auto act_comp = active_map[comp];
        const auto n_basis = n_basis_.comp_dimension[comp];
        const Size act_offset = comp_offset_[act_comp];
        const Size offset     = comp_offset_[comp];
        for (Size basis_i = 0; basis_i < n_basis;  ++basis_i)
        {
            const auto act_D_phi = D_phi.get_function_view(act_offset+basis_i);
            auto     inact_D_phi = D_phi.get_function_view(offset+basis_i);
            for (int qp = 0; qp < num_points; ++qp)
                for (int der = 0; der < n_ders; ++der)
                    inact_D_phi[qp](der)(comp) = act_D_phi[qp](der)(act_comp);
        }
    }
}




template <int dim, int range, int rank>
void
BSplineUniformQuadCache<dim, range, rank>::
evaluate_bspline_values(
    const ComponentContainer<TensorProductFunctionEvaluator<dim>> &elem_values,
    ValueTable<Value> &D_phi) const
{
    const auto n_points_direction = quad_.get_num_points_direction();
    const Size num_points = n_points_direction.flat_size();
    const TensorIndex<dim> der_tensor_id; // [0,0,..,0] tensor index
    for (int comp : elem_values.get_active_components_id())
    {
        auto &values = elem_values[comp];
        const int total_n_basis = n_basis_.comp_dimension[comp];
        const Size offset = comp_offset_[comp];

        for (int func_id = 0; func_id < total_n_basis; ++func_id)
        {
            auto D_phi_i = D_phi.get_function_view(offset + func_id);
            auto const &func = values.func_flat_to_tensor(func_id);
            for (int point_id = 0; point_id < num_points; ++point_id)
            {
                auto const &pts  = values.points_flat_to_tensor(point_id);
                D_phi_i[point_id](comp) = values.evaluate(der_tensor_id, func, pts);
            }
        } // end func_id loop
    } // end comp loop

    copy_to_inactive_components_values(elem_values.get_inactive_components_id(),
                                       elem_values.get_comp_map(), D_phi);
}



template <int dim, int range, int rank>
template <int order>
void
BSplineUniformQuadCache<dim, range, rank>::
evaluate_bspline_derivatives(
    const ComponentContainer<TensorProductFunctionEvaluator<dim>> &elem_values,
    ValueTable<Derivative<order>> &D_phi) const
{
    /*
     * This code computes any order of derivatives for a multivariate
     * B-spline on the current element
     * We use the formula
     * \partial_(\alpha_1,...,\alpha_n) B(qp) = \Pi d^{\alpha_i} B_i(qp_i)
     */

    Assert(D_phi.size() > 0, ExcEmptyObject());
    //  Assert(D_phi.get_num_functions() == this->get_num_basis(),
    //           ExcDimensionMismatch(D_phi.get_num_functions(),this->get_num_basis()));

    const auto n_points_direction = quad_.get_num_points_direction();
    const Size num_points = n_points_direction.flat_size();
    Assert(D_phi.get_num_points() == num_points,
           ExcDimensionMismatch(D_phi.get_num_points(),num_points));


    TensorFunctionDerivativesSymmetry<dim,order> sym;
    const auto n_der =  TensorFunctionDerivativesSymmetry<dim,order>::num_entries_eval;

    const auto &univariate_order = sym.univariate_order ;
    const auto &copy_indices = sym.copy_indices;

    for (int comp : elem_values.get_active_components_id())
    {
        auto &values = elem_values[comp];
        const int total_n_basis = n_basis_.comp_dimension[comp];
        const Size offset = comp_offset_[comp];

        for (int func_id = 0; func_id < total_n_basis; ++func_id)
        {

            auto D_phi_i = D_phi.get_function_view(offset + func_id);
            auto const &func = values.func_flat_to_tensor(func_id);
            for (int der_id = 0; der_id<n_der; ++der_id)
            {
                auto const &der_tensor_id = univariate_order[der_id];
                for (int point_id = 0; point_id < num_points; ++point_id)
                {
                    auto const &pts  = values.points_flat_to_tensor(point_id);
                    auto &der = D_phi_i[point_id];
                    der(copy_indices[der_id][0])(comp) = values.evaluate(der_tensor_id, func, pts);
                    for (int k=1; k<copy_indices[der_id].size(); ++k)
                        der(copy_indices[der_id][k])(comp) = der(copy_indices[der_id][0])(comp);
                }
            }

        }
    } // end comp loop

    copy_to_inactive_components<order>(elem_values.get_inactive_components_id(),
                                       elem_values.get_comp_map(), D_phi);
}
#endif



template<int dim_, int range_ , int rank_>
void
NURBSUniformQuadCache<dim_, range_, rank_>::
fill_element_cache(ElementIterator &elem)
{
    fill_element_cache(elem.get_accessor());
}



template<int dim_, int range_ , int rank_>
void
NURBSUniformQuadCache<dim_, range_, rank_>::
fill_element_cache(ElementAccessor &elem)
{
    Assert(false,ExcNotImplemented());
    AssertThrow(false,ExcNotImplemented());
#if 0
    base_t::fill_element_cache(elem);
    auto &cache = elem.get_elem_cache();

    const auto &elem_id = elem.get_tensor_index();
    auto values = splines1d_.get_element_values(elem_id);


    //--------------------------------------------------------------------------
    if (cache.flags_handler_.fill_values())
    {
        evaluate_bspline_values(values, cache.phi_);
        cache.flags_handler_.set_values_filled(true);
    }
    if (cache.flags_handler_.fill_gradients())
    {
        evaluate_bspline_derivatives<1>(values, cache.D1phi_);
        cache.flags_handler_.set_gradients_filled(true);
    }

    if (cache.flags_handler_.fill_hessians())
    {
        evaluate_bspline_derivatives<2>(values, cache.D2phi_);
        cache.flags_handler_.set_hessians_filled(true);
    }

    if (cache.flags_handler_.fill_divergences())
    {
        //TODO(pauletti, Sep 7, 2014): create a specialize exception
        Assert(cache.flags_handler_.gradients_filled(),
               ExcMessage("Divergence requires gradient to be filled."));

        auto D1  = cache.D1phi_.begin();
        auto div = cache.div_phi_.begin();
        auto end = cache.D1phi_.end();
        for (; D1 != end; ++D1, ++div)
            *div = trace(*D1);

        cache.flags_handler_.set_divergences_filled(true);
    }

    //--------------------------------------------------------------------------

    cache.set_filled(true);
#endif
}



template<int dim_, int range_ , int rank_>
void
NURBSUniformQuadCache<dim_, range_, rank_>::
print_info(LogStream &out) const
{
    Assert(false,ExcNotImplemented());
    AssertThrow(false,ExcNotImplemented());
#if 0

    out.begin_item("Grid Cache:");
    base_t::print_info(out);
    out.end_item();


    out.begin_item("One dimensional splines cache:");
    splines1d_.print_info(out);

    out.end_item();
#endif
}

IGA_NAMESPACE_CLOSE

#include <igatools/basis_functions/nurbs_uniform_quad_cache.inst>
