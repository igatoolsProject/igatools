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
vector<TensorIndex<size> >
partition(const int n)
{
    vector<TensorIndex<size>> v;
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
vector<TensorIndex<1> >
partition<1>(const int n)
{
    TensorIndex<1> arr(n);
    return vector<TensorIndex<1>>(1,arr);
}



template<>
vector<TensorIndex<0> >
partition<0>(const int n)
{
    return vector<TensorIndex<0>>();
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
            vector<TensorIndex<order>> v;
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
    special_array<TensorIndex<dim>, num_entries_eval> univariate_order;

    special_array<TensorIndex<order>, num_entries_eval> eval_indices;

    special_array<vector<TensorIndex<order>>, num_entries_eval> copy_indices;

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
template<class T>
void
BSplineElementHandler<dim_, range_, rank_>::
ResetDispatcher::
operator()(const T &quad1)
{
    Assert(grid_handler_ != nullptr,ExcNullPtr());
    grid_handler_->reset(FunctionFlags::to_grid_flags(flag_),quad1);


    const auto k = T::dim;

    Assert(flags_ != nullptr,ExcNullPtr());
    (*flags_)[k] = flag_;


    Assert(splines1d_ != nullptr,ExcNullPtr());


    Assert(space_ != nullptr,ExcNullPtr());

    const auto space_data = space_->space_data_;

//    const auto n_basis = space_->get_num_all_element_basis();
    const auto &degree = space_->get_degree();

    const auto &active_components_id = space_data->get_active_components_id();

    // number of intervals in the cartesian grid
    const auto n_inter = space_->get_grid()->get_num_intervals();

#ifndef NDEBUG
    for (const auto &intervals_id : intervals_id_directions_)
        Assert(!intervals_id.empty(),ExcEmptyObject());
#endif

    const auto &lengths = grid_handler_->get_grid()->get_element_lengths();

    for (auto &s_id: Topology::template elems_ids<k>())
    {

//        auto &g_cache = std::get<k>(*splines1d_)[s_id];
        auto &g_cache = boost::fusion::at_c<k>(*splines1d_).second[s_id];
        g_cache = GlobalCache(space_data->get_components_map());

        const auto quad = extend_sub_elem_quad<k,dim>(quad1, s_id);
        const auto &n_coords = quad.get_num_coords_direction();

        // Allocate space for the BasisValues1D
        for (auto comp : active_components_id)
        {
            for (const int dir : Topology::active_directions)
            {
                const auto &intervals_id = intervals_id_directions_[dir];

                const auto &n_pts = n_coords[dir];

                const auto n_basis_comp_dir = degree[comp][dir]+1;

                for (const int &interv_id : intervals_id)
                {
                    auto &splines1d = g_cache.entry(comp, dir, interv_id);
                    splines1d.resize(max_der, n_basis_comp_dir, n_pts);
                } // end loop interv_id
            } // end loop dir
        } // end loop comp

        /*
         * For each direction, interval and component we compute the 1D bspline
         * basis evaluate at the 1D component of the tensor product quadrature
         */
        const auto &bezier_op   = space_->operators_;
        const auto &end_interval = space_->end_interval_;

        using BasisValues = ComponentContainer<BasisValues1d>;
        const auto &deg_comp_map = degree.get_comp_map();
        BasisValues bernstein_values_internal(deg_comp_map);
        BasisValues bernstein_values_left(deg_comp_map);
        BasisValues bernstein_values_right(deg_comp_map);

        using LengthCompContainer = ComponentContainer<Points<dim>>;

        LengthCompContainer len_left(end_interval.get_comp_map());
        LengthCompContainer len_right(end_interval.get_comp_map());

        // First/last interval treatment
        for (const int dir : Topology::active_directions)
        {
            const auto &intervals_id = intervals_id_directions_[dir];

            const int id_interval_left  = 0;
            const int id_interval_right = n_inter[dir]-1;

            const auto &pt_coords_internal = quad.get_coords_direction(dir);
            const auto &len_internal = lengths.get_data_direction(dir);


            if (intervals_id.front() == id_interval_left) // processing the leftmost interval
            {
                vector<Real> pt_coords_left(n_coords[dir]);

                for (auto comp : bernstein_values_left.get_active_components_id())
                {
                    const Real alpha = end_interval[comp][dir].first;
                    const Real one_alpha = 1. - alpha;
                    len_left[comp][dir] = lengths.get_data_direction(dir)[id_interval_left]/alpha;

                    for (int ipt = 0 ; ipt < n_coords[dir] ; ++ipt)
                        pt_coords_left[ipt] = one_alpha + pt_coords_internal[ipt] * alpha;

                    resize_and_fill_bernstein_values(degree[comp][dir],pt_coords_left,bernstein_values_left[comp]);
                } // end loop comp
            } // end process_interval_left

            if (intervals_id.back() == id_interval_right) // processing the rightmost interval
            {
                vector<Real> pt_coords_right(n_coords[dir]);

                for (auto comp : bernstein_values_right.get_active_components_id())
                {
                    const Real alpha = end_interval[comp][dir].second;
                    len_right[comp][dir] = lengths.get_data_direction(dir)[id_interval_right]/alpha;

                    for (int ipt = 0 ; ipt < n_coords[dir] ; ++ipt)
                        pt_coords_right[ipt] = pt_coords_internal[ipt] * alpha;

                    resize_and_fill_bernstein_values(degree[comp][dir],pt_coords_right,bernstein_values_right[comp]);
                } // end loop comp
            } // end process_interval_right

            if (std::any_of(intervals_id.begin(),intervals_id.end(),
                            [&id_interval_left,&id_interval_right](int i)
        {
            return (i > id_interval_left) && (i < id_interval_right);
            }))
            {
                // processing the internal intervals
                for (auto comp : bernstein_values_internal.get_active_components_id())
                    resize_and_fill_bernstein_values(degree[comp][dir],pt_coords_internal,bernstein_values_internal[comp]);
            } // end process_interval_internal


            for (auto &inter : intervals_id)
            {
                for (auto comp : active_components_id)
                {
                    Real one_div_interval_length;
                    const BasisValues1d *berns_values_ptr = nullptr;
                    if (inter != id_interval_left && inter != id_interval_right)
                    {
                        // internal intervals
                        one_div_interval_length = 1.0 / len_internal[inter];
                        berns_values_ptr = &bernstein_values_internal[comp];
                    }
                    else if (inter == id_interval_left)
                    {
                        // first interval (i.e. left-most interval)
                        one_div_interval_length = 1.0 / len_left[comp][dir];
                        berns_values_ptr = &bernstein_values_left[comp];
                    }
                    else if (inter == id_interval_right)
                    {
                        // last interval (i.e. right-most interval)
                        one_div_interval_length = 1.0 / len_right[comp][dir];
                        berns_values_ptr = &bernstein_values_right[comp];
                    }

                    auto &basis = g_cache.entry(comp, dir, inter);
                    const auto &oper = bezier_op.get_operator(dir, inter, comp);

                    fill_interval_values(one_div_interval_length, oper, *berns_values_ptr, basis);
                } // end loop comp
            } // end loop inter

        } // end loop dir

    } // end loop s_id
}


template<int dim_, int range_ , int rank_>
void
BSplineElementHandler<dim_, range_, rank_>::
reset_selected_elements(
    const ValueFlags &flag,
    const eval_pts_variant &eval_points,
    const vector<int> elements_flat_id)
{
    reset_impl_.grid_handler_ = &(this->grid_handler_);
    reset_impl_.flag_ = flag;
    reset_impl_.flags_ = &flags_;
    reset_impl_.splines1d_ = &splines1d_;
    reset_impl_.space_ = this->get_bspline_space().get();


    //-------------------------------------------------
    // here we get the interval indices from the element indices
    Assert(!elements_flat_id.empty(),ExcEmptyObject());

    const auto grid = reset_impl_.space_->get_grid();
    array<set<int>,dim> intervals_id_unique;
    for (const auto elem_id : elements_flat_id)
    {
        const auto elem_tensor_id = grid->flat_to_tensor(elem_id);

        for (int dir = 0 ; dir < dim ; ++dir)
            intervals_id_unique[dir].insert(elem_tensor_id[dir]);
    }

    for (const int dir : Topology::active_directions)
    {
        reset_impl_.intervals_id_directions_[dir].assign(
            intervals_id_unique[dir].begin(),intervals_id_unique[dir].end());
    }
    //-------------------------------------------------

    boost::apply_visitor(reset_impl_, eval_points);
}







template<int dim_, int range_ , int rank_>
template<class T>
void
BSplineElementHandler<dim_, range_, rank_>::
InitCacheDispatcher::
operator()(const T &quad)
{
    Assert(grid_handler_ != nullptr,ExcNullPtr());
    Assert(elem_ != nullptr,ExcNullPtr());
    grid_handler_->template init_cache<T::k>(elem_->as_cartesian_grid_element_accessor());


    auto &cache = elem_->get_local_cache();
    if (cache == nullptr)
    {
        using VCache = typename BSplineElement<dim_,range_,rank_>::parent_t::Cache;

        using Cache = LocalCache<VCache>;

        cache = shared_ptr<Cache>(new Cache);
    }

    const auto n_basis = elem_->get_num_basis();//elem_->get_num_basis(DofProperties::active);
    const auto n_points = grid_handler_->template get_num_points<T::k>();
    const auto flag = (*flags_)[T::k];

    for (auto &s_id: Topology::template elems_ids<T::k>())
    {
        auto &s_cache = cache->template get_value_cache<T::k>(s_id);
        s_cache.resize(flag, n_points, n_basis);
    }
}

template<int dim_, int range_ , int rank_>
void
BSplineElementHandler<dim_, range_, rank_>::
init_cache(RefElementAccessor &elem, const topology_variant &topology)
{
    init_cache_impl_.grid_handler_ = &(this->grid_handler_);

    Assert(this->get_space() == elem.get_space(),
           ExcMessage("The element accessor and the element handler cannot have different spaces."));

    Assert(elem.get_space()->is_bspline(),ExcMessage("Not a BSplineElement."));
    init_cache_impl_.elem_ = &elem;

    init_cache_impl_.flags_ = &flags_;

    boost::apply_visitor(init_cache_impl_,topology);
}




template <int dim, int range, int rank>
void
BSplineElementHandler<dim, range, rank>::
FillCacheDispatcher::
copy_to_inactive_components_values(const vector<Index> &inactive_comp,
                                   const std::array<Index, n_components> &active_map,
                                   ValueTable<Value> &D_phi) const
{
    Assert(D_phi.get_num_functions() == elem_->get_num_basis(),
           ExcDimensionMismatch(D_phi.get_num_functions(),
                                elem_->get_num_basis()));

    const auto comp_offset = elem_->get_basis_offset();

    const Size n_points = D_phi.get_num_points();
    for (int comp : inactive_comp)
    {
        const auto act_comp = active_map[comp];
        const auto n_basis_comp = elem_->get_num_basis_comp(comp);
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
FillCacheDispatcher::
copy_to_inactive_components(const vector<Index> &inactive_comp,
                            const std::array<Index, n_components> &active_map,
                            ValueTable<Derivative<order>> &D_phi) const
{
    Assert(D_phi.get_num_functions() == elem_->get_num_basis(),
           ExcDimensionMismatch(D_phi.get_num_functions(),
                                elem_->get_num_basis()));

    const auto comp_offset = elem_->get_basis_offset();

    const Size n_points = D_phi.get_num_points();
    const Size n_ders = Derivative<order>::size;
    for (int comp : inactive_comp)
    {
        const auto act_comp = active_map[comp];
        const auto n_basis_comp = elem_->get_num_basis_comp(comp);
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
FillCacheDispatcher::
evaluate_bspline_values(
    const ComponentContainer<TensorProductFunctionEvaluator<dim>> &elem_values,
    ValueTable<Value> &D_phi) const
{
    Assert(D_phi.get_num_functions() == elem_->get_num_basis(),
           ExcDimensionMismatch(D_phi.get_num_functions(),
                                elem_->get_num_basis()));

    const auto comp_offset = elem_->get_basis_offset();

    const Size n_points = D_phi.get_num_points();
    const TensorIndex<dim> der_tensor_id; // [0,0,..,0] tensor index
    for (int comp : elem_values.get_active_components_id())
    {
        auto &values = elem_values[comp];
        const int n_basis_comp = elem_->get_num_basis_comp(comp);
        const Size offset = comp_offset[comp];

        for (int func_id = 0; func_id < n_basis_comp; ++func_id)
        {
            auto D_phi_i = D_phi.get_function_view(offset + func_id);
            auto const &func = values.func_flat_to_tensor(func_id);
            for (int point_id = 0; point_id < n_points; ++point_id)
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
BSplineElementHandler<dim, range, rank>::
FillCacheDispatcher::
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
    const Size n_points = D_phi.get_num_points();


    Assert(D_phi.get_num_functions() == elem_->get_num_basis(),
           ExcDimensionMismatch(D_phi.get_num_functions(),
                                elem_->get_num_basis()));

    const auto comp_offset = elem_->get_basis_offset();

    TensorFunctionDerivativesSymmetry<dim,order> sym;
    const auto n_der =  TensorFunctionDerivativesSymmetry<dim,order>::num_entries_eval;

    const auto &univariate_order = sym.univariate_order ;
    const auto &copy_indices = sym.copy_indices;

    for (int comp : elem_values.get_active_components_id())
    {
        auto &values = elem_values[comp];
        const int n_basis_comp = elem_->get_num_basis_comp(comp);
        const Size offset = comp_offset[comp];

        for (int func_id = 0; func_id < n_basis_comp; ++func_id)
        {
            auto D_phi_i = D_phi.get_function_view(offset + func_id);
            auto const &func = values.func_flat_to_tensor(func_id);
            for (int der_id = 0; der_id<n_der; ++der_id)
            {
                const auto &copy_indices_der = copy_indices[der_id];
                const auto copy_indices_der_size = copy_indices_der.size();

                auto const &der_tensor_id = univariate_order[der_id];
                for (int point_id = 0; point_id < n_points; ++point_id)
                {
                    auto const &pts  = values.points_flat_to_tensor(point_id);
                    auto &der = D_phi_i[point_id];
                    der(copy_indices_der[0])(comp) = values.evaluate(der_tensor_id, func, pts);
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
template<class T>
void
BSplineElementHandler<dim_, range_, rank_>::
FillCacheDispatcher::
operator()(const T &topology)
{
    Assert(grid_handler_ != nullptr,ExcNullPtr());
    grid_handler_->template fill_cache<T::k>(elem_->as_cartesian_grid_element_accessor(),j_);

    Assert(splines1d_ != nullptr, ExcNullPtr());
//    const auto &g_cache = std::get<T::k>(*splines1d_)[j_];
    const auto &g_cache = boost::fusion::at_c<T::k>(*splines1d_).second[j_];


    Assert(elem_ != nullptr, ExcNullPtr());
    Assert(elem_->get_local_cache() != nullptr, ExcNullPtr());
    auto &cache = elem_->get_local_cache()->template get_value_cache<T::k>(j_);

    const auto &elem_t_index = elem_->get_tensor_index();

    auto &flags = cache.flags_handler_;
    auto val_1d = g_cache.get_element_values(elem_t_index);
    if (flags.template fill<_Value>())
    {
        auto &values = cache.template get_der<_Value>();
        evaluate_bspline_values(val_1d, values);
        flags.template set_filled<_Value>(true);
    }
    if (flags.template fill<_Gradient>())
    {
        auto &values = cache.template get_der<_Gradient>();
        evaluate_bspline_derivatives<1>(val_1d, values);
        flags.template set_filled<_Gradient>(true);
    }
    if (flags.template fill<_Hessian>())
    {
        auto &values = cache.template get_der<_Hessian>();
        evaluate_bspline_derivatives<2>(val_1d, values);
        flags.template set_filled<_Hessian>(true);
    }
    if (flags.template fill<_Divergence>())
    {
        eval_divergences_from_gradients(
            cache.template get_der<_Gradient>(),
            cache.template get_der<_Divergence>());
        flags.template set_filled<_Divergence>(true);
    }


    cache.set_filled(true);
}



template<int dim_, int range_ , int rank_>
void
BSplineElementHandler<dim_, range_, rank_>::
fill_cache(RefElementAccessor &elem, const topology_variant &topology, const int j)
{
    fill_cache_impl_.grid_handler_ = &(this->grid_handler_);

    fill_cache_impl_.j_ = j;
    fill_cache_impl_.splines1d_ = &splines1d_;

    Assert(this->get_space() == elem.get_space(),
           ExcMessage("The element accessor and the element handler cannot have different spaces."));

    Assert(elem.get_space()->is_bspline(),ExcMessage("Not a BSplineElement."));
    fill_cache_impl_.elem_ = &elem;

    boost::apply_visitor(fill_cache_impl_,topology);
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

    out.begin_item("Splines 1D Cache:");
    cacheutils::print_caches(splines1d_, out);
//    splines1d_.print_info(out);
    out.end_item();
}

IGA_NAMESPACE_CLOSE

#include <igatools/basis_functions/bspline_element_handler.inst>
