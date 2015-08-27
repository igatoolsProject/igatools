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

#include <igatools/basis_functions/bspline_element.h>
#include <igatools/basis_functions/bspline_element_handler.h>
//#include <igatools/basis_functions/bernstein_basis.h>


#include <algorithm>
using std::shared_ptr;

using std::set;

IGA_NAMESPACE_OPEN

namespace
{
//TODO(pauletti, Sep 9, 2014): this class seems to have a more general use put in own
// file under group tensor product values utilities
/**
 * @class DerivativeEvaluationSymmetryManager
 *
 * This class is used to manage the symmetries in the evaluation of a k-th order
 *  derivative of a scalar function \f$ f \mathbb{R}^d -to \mathbb{R} \f$
 *  according with the Schwarz theorem (equality of mixed partials).
 *
 * The total number of derivatives is given by \f$ d^k \f$ where the number of
 * different values
 * is \f$ \binorm{d-1+k}{d-1} = \binorm{d-1+k}{k} \f$.
 *
 */
template<int size>
SafeSTLVector<TensorIndex<size> >
partition(const int n)
{
    SafeSTLVector<TensorIndex<size>> v;
    TensorIndex<size> arr(0);

    arr[0] = n;
    v.push_back(arr);

    for (int j=1; j<n+1; ++j)
    {
        auto w = partition<size-1>(j);
        for (auto a : w)
        {
            arr[0] = n-j;
            std::copy(a.begin(), a.end(), arr.begin()+1);
            v.push_back(arr);
        }
    }
    return v;
}

template<>
SafeSTLVector<TensorIndex<1> >
partition<1>(const int n)
{
    TensorIndex<1> arr(n);
    return SafeSTLVector<TensorIndex<1>>(1,arr);
}



template<>
SafeSTLVector<TensorIndex<0> >
partition<0>(const int n)
{
    return SafeSTLVector<TensorIndex<0>>();
}



template<int dim, int order>
class TensorFunctionDerivativesSymmetry
{
public:
//    static const int num_entries_total = pow(dim,order);
    static const int num_entries_eval = constexpr_binomial_coefficient(dim-1+order,order);

    TensorFunctionDerivativesSymmetry()
    {
        auto uni_indices = partition<dim>(order);
        std::copy(uni_indices.begin(), uni_indices.end(), univariate_order.begin());



        for (int j=0; j<num_entries_eval; ++j)
        {
            auto &der_ind = eval_indices[j];
            int s=0;
            for (int dir=0; dir<dim; ++dir)
            {
                for (int l=0; l<uni_indices[j][dir]; ++l)
                    der_ind[s+l] = dir;
                s += uni_indices[j][dir];
            }

            auto ind = sequence<order>();
            SafeSTLVector<TensorIndex<order>> v;
            do
            {
                TensorIndex<order> ti;
                for (int i=0; i<order; ++i)
                    ti[i] = eval_indices[j][ind[i]];
                v.push_back(ti);
            }
            while (std::next_permutation(ind.begin(),ind.end()));

            auto it = std::unique(v.begin(), v.end());
            v.resize(std::distance(v.begin(),it));

            copy_indices[j] = v;
        }
    }

    void print_info(LogStream &out) const
    {
        out.begin_item("univariate derivative orders:");
        univariate_order.print_info(out);
        out.end_item();

        out.begin_item("Assigment indices:");
        eval_indices.print_info(out);
        out.end_item();

        out.begin_item("all equal indices indices:");
        copy_indices.print_info(out);
        out.end_item();
    }
    SafeSTLArray<TensorIndex<dim>, num_entries_eval> univariate_order;

    SafeSTLArray<TensorIndex<order>, num_entries_eval> eval_indices;

    SafeSTLArray<SafeSTLVector<TensorIndex<order>>, num_entries_eval> copy_indices;

};


}; // of the namespace




template<int dim_, int range_ , int rank_>
BSplineElementHandler<dim_, range_, rank_>::
BSplineElementHandler(shared_ptr<const Space> space)
    :
    base_t(space)
{
    Assert(space != nullptr, ExcNullPtr());
}


template<int dim_, int range_ , int rank_>
auto
BSplineElementHandler<dim_, range_, rank_>::
create(std::shared_ptr<const Space> space) -> std::shared_ptr<self_t>
{
    return std::shared_ptr<self_t>(new self_t(space));
}


template<int dim_, int range_ , int rank_>
template<int sdim>
void
BSplineElementHandler<dim_, range_, rank_>::
SetFlagDispatcher::
operator()(const Topology<sdim> &topology)
{
    using GridFlags = grid_element::Flags;
    //TODO (martinelli, Aug 27, 2015): select the proper grid flags depending on the BSpline element flags
    const auto grid_flags = GridFlags::point |
                            GridFlags::w_measure;
    grid_handler_.template set_flags<sdim>(grid_flags);

    flags_[sdim] = flag_in_;
}




template<int dim_, int range_ , int rank_>
template<int sdim>
void
BSplineElementHandler<dim_, range_, rank_>::
InitCacheDispatcher::
operator()(const std::shared_ptr<const Quadrature<sdim>> &quad)
{
    grid_handler_.template init_cache<sdim>(elem_.get_grid_element(),quad);

    auto &cache = elem_.get_all_sub_elems_cache();
    if (cache == nullptr)
    {
        using VCache = typename BSplineElement<dim_,range_,rank_>::parent_t::Cache;

        using Cache = AllSubElementsCache<VCache>;
        cache = std::make_shared<Cache>();
    }

    const auto n_basis = elem_.get_max_num_basis();
    const auto n_points = quad->get_num_points();
    const auto flag = flags_[sdim];

    for (auto &s_id: UnitElement<dim_>::template elems_ids<sdim>())
    {
        auto &s_cache = cache->template get_sub_elem_cache<sdim>(s_id);
        s_cache.resize(flag, n_points, n_basis);
    }
}





template <int dim, int range, int rank>
void
BSplineElementHandler<dim, range, rank>::
FillCacheDispatcherNoGlobalCache::
copy_to_inactive_components_values(const SafeSTLVector<Index> &inactive_comp,
                                   const SafeSTLArray<Index, n_components> &active_map,
                                   ValueTable<Value> &D_phi) const
{
    Assert(D_phi.get_num_functions() == elem_.get_max_num_basis(),
           ExcDimensionMismatch(D_phi.get_num_functions(),
                                elem_.get_max_num_basis()));

    const auto &bsp_elem = dynamic_cast<BSplineElement<dim,range,rank> &>(elem_);
    const auto comp_offset = bsp_elem.get_basis_offset();

    const Size n_points = D_phi.get_num_points();
    for (int comp : inactive_comp)
    {
        const auto act_comp = active_map[comp];
        const auto n_basis_comp = bsp_elem.get_num_basis_comp(comp);
        const Size act_offset = comp_offset[act_comp];
        const Size offset     = comp_offset[comp];
        for (Size basis_i = 0; basis_i < n_basis_comp;  ++basis_i)
        {
            const auto act_D_phi = D_phi.get_function_view(act_offset+basis_i);
            auto inact_D_phi = D_phi.get_function_view(offset+basis_i);
            for (int qp = 0; qp < n_points; ++qp)
                inact_D_phi[qp](comp) = act_D_phi[qp](act_comp);
        }
    }
}



template <int dim, int range, int rank>
template <int order>
void
BSplineElementHandler<dim, range, rank>::
FillCacheDispatcherNoGlobalCache::
copy_to_inactive_components(const SafeSTLVector<Index> &inactive_comp,
                            const SafeSTLArray<Index, n_components> &active_map,
                            ValueTable<Derivative<order>> &D_phi) const
{
    Assert(D_phi.get_num_functions() == elem_.get_max_num_basis(),
           ExcDimensionMismatch(D_phi.get_num_functions(),
                                elem_.get_max_num_basis()));

    const auto &bsp_elem = dynamic_cast<BSplineElement<dim,range,rank> &>(elem_);
    const auto comp_offset = bsp_elem.get_basis_offset();

    const Size n_points = D_phi.get_num_points();
    const Size n_ders = Derivative<order>::size;
    for (int comp : inactive_comp)
    {
        const auto act_comp = active_map[comp];
        const auto n_basis_comp = bsp_elem.get_num_basis_comp(comp);
        const Size act_offset = comp_offset[act_comp];
        const Size offset     = comp_offset[comp];
        for (Size basis_i = 0; basis_i < n_basis_comp;  ++basis_i)
        {
            const auto act_D_phi = D_phi.get_function_view(act_offset+basis_i);
            auto     inact_D_phi = D_phi.get_function_view(offset+basis_i);
            for (int qp = 0; qp < n_points; ++qp)
                for (int der = 0; der < n_ders; ++der)
                    inact_D_phi[qp](der)(comp) = act_D_phi[qp](der)(act_comp);
        }
    }
}







template <int dim, int range, int rank>
void
BSplineElementHandler<dim, range, rank>::
FillCacheDispatcherNoGlobalCache::
evaluate_bspline_values(
    const ComponentContainer<std::unique_ptr<const TensorProductFunctionEvaluator<dim>>> &elem_values,
    ValueTable<Value> &D_phi) const
{
    Assert(D_phi.get_num_functions() == elem_.get_max_num_basis(),
           ExcDimensionMismatch(D_phi.get_num_functions(),
                                elem_.get_max_num_basis()));

    const auto &bsp_elem = dynamic_cast<BSplineElement<dim,range,rank> &>(elem_);
    const auto comp_offset = bsp_elem.get_basis_offset();


    const Size n_points = D_phi.get_num_points();
    const TensorIndex<dim> der_tensor_id; // [0,0,..,0] tensor index
    for (int comp : elem_values.get_active_components_id())
    {
        const auto &values = *elem_values[comp];
        const int n_basis_comp = bsp_elem.get_num_basis_comp(comp);
        const Size offset = comp_offset[comp];

        for (int func_id = 0; func_id < n_basis_comp; ++func_id)
        {
            auto D_phi_i = D_phi.get_function_view(offset + func_id);
            auto const &func_tensor_id = values.func_flat_to_tensor(func_id);
            for (int point_id = 0; point_id < n_points; ++point_id)
                D_phi_i[point_id](comp) = values.evaluate(der_tensor_id, func_tensor_id, point_id);
        } // end func_id loop
    } // end comp loop

    copy_to_inactive_components_values(elem_values.get_inactive_components_id(),
                                       elem_values.get_comp_map(), D_phi);
}




template <int dim, int range, int rank>
template <int order>
void
BSplineElementHandler<dim, range, rank>::
FillCacheDispatcherNoGlobalCache::
evaluate_bspline_derivatives(
    const ComponentContainer<std::unique_ptr<const TensorProductFunctionEvaluator<dim>>> &elem_values,
    ValueTable<Derivative<order>> &D_phi) const
{
    /*
     * This code computes any order of derivatives for a multivariate
     * B-spline on the current element
     * We use the formula
     * \partial_(\alpha_1,...,\alpha_n) B(qp) = \Pi d^{\alpha_i} B_i(qp_i)
     */

    Assert(D_phi.size() > 0, ExcEmptyObject());
    const Size n_points = D_phi.get_num_points();


    Assert(D_phi.get_num_functions() == elem_.get_max_num_basis(),
           ExcDimensionMismatch(D_phi.get_num_functions(),
                                elem_.get_max_num_basis()));

    const auto &bsp_elem = dynamic_cast<BSplineElement<dim,range,rank> &>(elem_);
    const auto comp_offset = bsp_elem.get_basis_offset();

    TensorFunctionDerivativesSymmetry<dim,order> sym;
    const auto n_der = TensorFunctionDerivativesSymmetry<dim,order>::num_entries_eval;

    const auto &univariate_order = sym.univariate_order ;
    const auto &copy_indices = sym.copy_indices;

    for (int comp : elem_values.get_active_components_id())
    {
        const auto &values = *elem_values[comp];
        const int n_basis_comp = bsp_elem.get_num_basis_comp(comp);
        const Size offset = comp_offset[comp];

        for (int func_id = 0; func_id < n_basis_comp; ++func_id)
        {
            auto D_phi_i = D_phi.get_function_view(offset + func_id);
            auto const &func_tensor_id = values.func_flat_to_tensor(func_id);
            for (int der_id = 0; der_id<n_der; ++der_id)
            {
                const auto &copy_indices_der = copy_indices[der_id];
                const auto copy_indices_der_size = copy_indices_der.size();

                auto const &der_tensor_id = univariate_order[der_id];
                for (int point_id = 0; point_id < n_points; ++point_id)
                {
                    auto &der = D_phi_i[point_id];
                    der(copy_indices_der[0])(comp) = values.evaluate(der_tensor_id, func_tensor_id, point_id);
                    for (int k=1; k<copy_indices_der_size; ++k)
                        der(copy_indices_der[k])(comp) = der(copy_indices_der[0])(comp);
                }
            }

        }
    } // end comp loop

    copy_to_inactive_components<order>(elem_values.get_inactive_components_id(),
                                       elem_values.get_comp_map(), D_phi);
}



template<int dim_, int range_ , int rank_>
template<int sdim>
void
BSplineElementHandler<dim_, range_, rank_>::
FillCacheDispatcherNoGlobalCache::
operator()(const Topology<sdim> &topology)
{
    auto &grid_elem = elem_.get_grid_element();
    grid_handler_.template fill_cache<sdim>(grid_elem,s_id_);



    //--------------------------------------------------------------------------------------
    // filling the 1D cache --- begin

    const auto &grid = *grid_elem.get_grid();
    const auto n_inter = grid.get_num_intervals();

    const auto elem_size = grid_elem.template get_side_lengths<dim>(0);
    const auto elem_tensor_id = grid_elem.get_index();

    const auto &bsp_space = dynamic_cast<const Space &>(*elem_.get_space());

    const auto &space_data = *bsp_space.space_data_;

    const auto &degree = bsp_space.get_degree_table();

    const auto &active_components_id = space_data.get_active_components_id();

    const auto quad_in_cache = grid_elem.template get_quadrature<sdim>();

    const auto quad = extend_sub_elem_quad<sdim,dim>(*quad_in_cache, s_id_);

    using BasisValues1dTable = ComponentContainer<SafeSTLArray<BasisValues1d,dim>>;
    BasisValues1dTable splines_derivatives_1D_table(space_data.get_components_map());

    const auto &n_coords = quad.get_num_coords_direction();

    /*
     * For each direction, interval and component we compute the 1D bspline
     * basis evaluate at the 1D component of the tensor product quadrature
     */
    const auto &bezier_op   = bsp_space.operators_;
    const auto &end_interval = bsp_space.end_interval_;

    using BasisValues = ComponentContainer<BasisValues1d>;
    const auto &deg_comp_map = degree.get_comp_map();
    BasisValues bernstein_values(deg_comp_map);

    for (const int dir : UnitElement<dim>::active_directions)
    {
        const auto &pt_coords_internal = quad.get_coords_direction(dir);

        const auto len = elem_size[dir];

        const auto interval_id = elem_tensor_id[dir];

        Real alpha;

        const int n_pts_1D = n_coords[dir];

        SafeSTLVector<Real> pt_coords_boundary(n_pts_1D);

        const SafeSTLVector<Real> *pt_coords_ptr = nullptr;

        for (auto comp : active_components_id)
        {
            const int deg = degree[comp][dir];

            auto &splines_derivatives_1D = splines_derivatives_1D_table[comp][dir];
            splines_derivatives_1D.resize(MAX_NUM_DERIVATIVES,deg+1,n_pts_1D);

            if (interval_id == 0) // processing the leftmost interval
            {
                // first interval (i.e. left-most interval)

                alpha = end_interval[comp][dir].first;
                const Real one_minus_alpha = 1. - alpha;

                for (int ipt = 0 ; ipt < n_pts_1D ; ++ipt)
                    pt_coords_boundary[ipt] = one_minus_alpha +
                                              pt_coords_internal[ipt] * alpha;

                pt_coords_ptr = &pt_coords_boundary;
            } // end process_interval_left
            else if (interval_id == n_inter[dir]-1) // processing the rightmost interval
            {
                // last interval (i.e. right-most interval)

                alpha = end_interval[comp][dir].second;

                for (int ipt = 0 ; ipt < n_pts_1D ; ++ipt)
                    pt_coords_boundary[ipt] = pt_coords_internal[ipt] *
                                              alpha;

                pt_coords_ptr = &pt_coords_boundary;
            } // end process_interval_right
            else
            {
                // internal interval

                alpha = 1.0;

                pt_coords_ptr = &pt_coords_internal;
            } // end process_interval_internal


            const Real alpha_div_interval_length = alpha / len;

            const auto &oper = bezier_op.get_operator(dir,interval_id,comp);

            //------------------------------------------------------------
            //resize_and_fill_bernstein_values
            bernstein_values[comp].resize(MAX_NUM_DERIVATIVES,deg+1,n_pts_1D);
            for (int order = 0; order < MAX_NUM_DERIVATIVES; ++order)
            {
                auto &berns = bernstein_values[comp].get_derivative(order);
                berns = BernsteinBasis::derivative(order, deg,*pt_coords_ptr);

                auto &splines = splines_derivatives_1D.get_derivative(order);
                splines = oper.scale_action(
                              std::pow(alpha_div_interval_length, order),
                              berns);
            }
            //------------------------------------------------------------

        } // end loop comp

    } // end loop dir
    //
    // filling the 1D cache --- end
    //-------------------------------------------------------------------------------



    //-------------------------------------------------------------------------------
    // Multi-variate spline evaluation from 1D values --- begin
    using TPFE = const TensorProductFunctionEvaluator<dim>;
    ComponentContainer<std::unique_ptr<TPFE>> val_1d(splines_derivatives_1D_table.get_comp_map());

    SafeSTLArray<BasisValues1dConstView, dim> values_1D;
    for (auto c : val_1d.get_active_components_id())
    {
        const auto &value = splines_derivatives_1D_table[c];
        for (int i = 0 ; i < dim_ ; ++i)
            values_1D[i] = BasisValues1dConstView(value[i]);

        val_1d[c] = std::make_unique<TPFE>(quad,values_1D);
    }
    // Multi-variate spline evaluation from 1D values --- end
    //-------------------------------------------------------------------------------



    //-------------------------------------------------------------------------------
    auto &all_sub_elems_cache = elem_.get_all_sub_elems_cache();
    Assert(all_sub_elems_cache != nullptr, ExcNullPtr());
    auto &sub_elem_cache = all_sub_elems_cache->template get_sub_elem_cache<sdim>(s_id_);


    using Elem = SpaceElement<dim_,0,range_,rank_,Transformation::h_grad>;
    using _Value      = typename Elem::_Value;
    using _Gradient   = typename Elem::_Gradient;
    using _Hessian    = typename Elem::_Hessian;
    using _Divergence = typename Elem::_Divergence;

    if (sub_elem_cache.template status_fill<_Value>())
    {
        auto &values = sub_elem_cache.template get_data<_Value>();
        evaluate_bspline_values(val_1d, values);
        sub_elem_cache.template set_status_filled<_Value>(true);
    }
    if (sub_elem_cache.template status_fill<_Gradient>())
    {
        auto &values = sub_elem_cache.template get_data<_Gradient>();
        evaluate_bspline_derivatives<1>(val_1d, values);
        sub_elem_cache.template set_status_filled<_Gradient>(true);
    }
    if (sub_elem_cache.template status_fill<_Hessian>())
    {
        auto &values = sub_elem_cache.template get_data<_Hessian>();
        evaluate_bspline_derivatives<2>(val_1d, values);
        sub_elem_cache.template set_status_filled<_Hessian>(true);
    }
    if (sub_elem_cache.template status_fill<_Divergence>())
    {
        eval_divergences_from_gradients(
            sub_elem_cache.template get_data<_Gradient>(),
            sub_elem_cache.template get_data<_Divergence>());
        sub_elem_cache.template set_status_filled<_Divergence>(true);
    }

    sub_elem_cache.set_filled(true);
    //-------------------------------------------------------------------------------
}


template<int dim_, int range_ , int rank_>
auto
BSplineElementHandler<dim_, range_, rank_>::
get_bspline_space() const -> std::shared_ptr<const Space>
{
    auto bsp_space = std::dynamic_pointer_cast<const Space>(this->get_space());
    Assert(bsp_space != nullptr,ExcNullPtr());
    return bsp_space;
}





template<int dim_, int range_ , int rank_>
void
BSplineElementHandler<dim_, range_, rank_>::
print_info(LogStream &out) const
{
    /*
    out.begin_item("Grid Cache:");
    this->grid_handler_.print_info(out);
    out.end_item();
    //*/

#if 0
    out.begin_item("Splines 1D Cache:");
    cacheutils::print_caches(splines1d_, out);
//    splines1d_.print_info(out);
    out.end_item();
#endif
}



#if 0
template<int dim_, int range_ , int rank_>
BSplineElementHandler<dim_, range_, rank_>::
GlobalCache::
GlobalCache(const std::shared_ptr<const Quadrature<dim>> &quad, const ComponentMap &component_map)
    :
    quad_(quad),
    basis_values_1d_table_(BasisValues1dTable(component_map))
{}

template<int dim_, int range_ , int rank_>
auto
BSplineElementHandler<dim_, range_, rank_>::
GlobalCache::
entry(const int comp, const int dir, const Index interval_id) -> BasisValues1d &
{
    return basis_values_1d_table_[comp][dir][interval_id];
}


template<int dim_, int range_ , int rank_>
void
BSplineElementHandler<dim_, range_, rank_>::
GlobalCache::
print_info(LogStream &out) const
{
    using std::to_string;
    for (const auto comp : basis_values_1d_table_.get_active_components_id())
    {
        out.begin_item("Active Component ID: " + to_string(comp));

        for (const int dir : UnitElement<dim_>::active_directions)
        {
            out.begin_item("Direction : " + to_string(dir));

            for (const auto &interv_id_and_basis : basis_values_1d_table_[comp][dir])
            {
                const auto interval_id = interv_id_and_basis.first;
                const auto &basis = interv_id_and_basis.second;

                out.begin_item("Interval ID: " + to_string(interval_id));
                basis.print_info(out);
                out.end_item();
            }
            out.end_item();
        } // end loop dir
        out.end_item();
    } // end loop comp
}


template<int dim_, int range_ , int rank_>
auto
BSplineElementHandler<dim_, range_, rank_>::
GlobalCache::
get_element_values(const TensorIndex<dim> &elem_tensor_id) const
-> ComponentContainer<std::unique_ptr<const TensorProductFunctionEvaluator<dim>>>
{
    ComponentContainer<std::unique_ptr<const TensorProductFunctionEvaluator<dim>> >
    result(basis_values_1d_table_.get_comp_map());

    SafeSTLArray<BasisValues1dConstView, dim> values_1D;
    for (auto c : result.get_active_components_id())
    {
        const auto &value = basis_values_1d_table_[c];
        for (int i = 0 ; i < dim_ ; ++i)
            values_1D[i] = BasisValues1dConstView(value[i].at(elem_tensor_id[i]));

        result[c] = std::make_unique<const TensorProductFunctionEvaluator<dim>>(*this->quad_,values_1D);
    }
    return result;
}
#endif

#ifdef SERIALIZATION
template<int dim_, int range_ , int rank_>
template<class Archive>
void
BSplineElementHandler<dim_, range_, rank_>::
serialize(Archive &ar, const unsigned int version)
{
    ar &boost::serialization::make_nvp("BSplineElementHandler_base_t",
                                       boost::serialization::base_object<base_t>(*this));

    ar &boost::serialization::make_nvp("flags_",flags_);
#if 0
    ar &boost::serialization::make_nvp("splines1d_",splines1d_);
#endif
}
#endif


IGA_NAMESPACE_CLOSE

#include <igatools/basis_functions/bspline_element_handler.inst>
