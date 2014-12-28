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
#include <igatools/basis_functions/bernstein_basis.h>
#include <igatools/utils/multi_array_utils.h>

#include <algorithm>
using std::shared_ptr;

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
    base_t(space->get_grid()),
    space_(space),
    n_basis_(space_->get_num_all_element_basis())
{

    // Compute the component offsets
    comp_offset_[0] = 0;
    for (int j = 1; j < Space::n_components; ++j)
        comp_offset_[j] = comp_offset_[j-1] + n_basis_.get_component_size(j-1);
}



template<int dim_, int range_ , int rank_>
template<int k>
void
BSplineElementHandler<dim_, range_, rank_>::
reset(const ValueFlags flag,
      const Quadrature<k> &quad1)
{
    base_t::template reset<k>(FunctionFlags::to_grid_flags(flag), quad1);
    flags_[k] = flag;

    for (auto &s_id: UnitElement<dim>::template elems_ids<k>())
    {
        auto &g_cache = std::get<k>(splines1d_)[s_id];
        g_cache.clear();
        g_cache.resize(space_->get_grid()->get_num_intervals(),
                       BasisValues(space_->get_components_map()));
        const auto &n_inter = space_->get_grid()->get_num_intervals();
        const auto quad = extend_sub_elem_quad<k,dim>(quad1, s_id);
        const auto &n_points = quad.get_num_points_direction();

        // Allocate space for the BasisValues1D
        for (int dir = 0 ; dir < dim ; ++dir)
        {
            const auto &n_pts = n_points[dir];
            for (int j = 0 ; j < n_inter[dir] ; ++j)
            {
                auto &splines1d = g_cache.entry(dir, j);
                for (auto comp : splines1d.get_active_components_id())
                    splines1d[comp].resize(max_der, n_basis_[comp][dir], n_pts);
            }
        }

        /*
         * For each direction, interval and component we compute the 1D bspline
         * basis evaluate at the 1D component of the tensor product quadrature
         */
        const auto &degree      = space_->get_degree();
        const auto &bezier_op   = space_->operators_;
        const auto &end_interval = space_->end_interval_;
        const auto &points      = quad.get_points();
        const auto &lengths = this->lengths_;

        BasisValues bernstein_values(n_basis_.get_comp_map());


        ComponentContainer<Points<dim>> len_left(end_interval.get_comp_map());
        ComponentContainer<TensorProductArray<dim>>
		points_left(end_interval.get_comp_map(), points);
        LogStream out1;
        for (auto comp : points_left.get_active_components_id())
        {
        	Points<dim> dilate;
        	Points<dim> translate;
        	for (int dir=0; dir<dim; ++dir)
        	{
        		const Real alpha = end_interval[comp][dir].first;
        		const Real one_alpha = 1. - alpha;
        		dilate[dir] = alpha;
        		translate[dir] = one_alpha;
        		len_left[comp][dir] = lengths.get_data_direction(dir)[0]*alpha;
        	}
        	out1 << dilate << std::endl;
        	out1 << len_left[comp] << std::endl;
//        	out1 << translate << std::endl;
//        	points_left.print_info(out1);
        	points_left[comp].dilate_translate(dilate, translate);
//        	points_left.print_info(out1);
        }

        ComponentContainer<Points<dim>> len_right(end_interval.get_comp_map());
        ComponentContainer<TensorProductArray<dim>>
		points_right(end_interval.get_comp_map(), points);
        for (auto comp : points_right.get_active_components_id())
        {
        	Points<dim> dilate;
        	for (int dir=0; dir<dim; ++dir)
        	{
        		const Real alpha = end_interval[comp][dir].second;
        		dilate[dir] = alpha;
        		len_right[comp][dir] = lengths.get_data_direction(dir)[n_inter[dir]-1]*alpha;
        	}
        	points_right[comp].dilate(dilate);
        }


// Left interval treatment
        for (int dir = 0 ; dir < dim ; ++dir)
        {
        	for (auto comp : bernstein_values.get_active_components_id())
        	{
        		const int deg = degree[comp][dir];
        		bernstein_values[comp].resize(max_der, deg+1, n_points[dir]);
        		const auto &pt_coords = points_left[comp].get_data_direction(dir);
        		for (int order = 0; order < max_der; ++order)
        			bernstein_values[comp].get_derivative(order) =
        					BernsteinBasis::derivative(order, deg, pt_coords);
        	}

        	//const auto &inter_lengths = lengths.get_data_direction(dir);
        	const int inter = 0;
        	{
        		auto &splines1d = g_cache.entry(dir, inter);
        		for (auto comp : splines1d.get_active_components_id())
        		{
        			const auto &berns_values = bernstein_values[comp];
        			auto &basis = splines1d[comp];
        			const auto &oper = bezier_op.get_operator(dir, inter, comp);

        			const Real one_div_size = 1.0 / len_left[comp][dir];
        			LogStream out2;
        			out2 << one_div_size;
        			fill_interval_values(one_div_size, oper, berns_values, basis);
        			berns_values.print_info(out2);
        			basis.print_info(out2);
        		}
        	}
        }


        // Right interval treatment
        for (int dir = 0 ; dir < dim ; ++dir)
        {
        	for (auto comp : bernstein_values.get_active_components_id())
        	{
        		const int deg = degree[comp][dir];
        		bernstein_values[comp].resize(max_der, deg+1, n_points[dir]);
        		const auto &pt_coords = points_right[comp].get_data_direction(dir);
        		for (int order = 0; order < max_der; ++order)
        			bernstein_values[comp].get_derivative(order) =
        					BernsteinBasis::derivative(order, deg, pt_coords);
        	}

        	//const auto &inter_lengths = lengths.get_data_direction(dir);
        	const int inter = n_inter[dir]-1;
        	{
        		auto &splines1d = g_cache.entry(dir, inter);
        		for (auto comp : splines1d.get_active_components_id())
        		{
        			const auto &berns_values = bernstein_values[comp];
        			auto &basis = splines1d[comp];
        			const auto &oper = bezier_op.get_operator(dir, inter, comp);

        			const Real one_div_size = 1.0 / len_right[comp][dir];
        			fill_interval_values(one_div_size, oper, berns_values, basis);
        		}
        	}
        }



        for (int dir = 0 ; dir < dim ; ++dir)
        {
        	const int inter_begin = 1;
        	const int inter_end = n_inter[dir];//-1;

        	for (auto comp : bernstein_values.get_active_components_id())
            {
                const int deg = degree[comp][dir];
                bernstein_values[comp].resize(max_der, deg+1, n_points[dir]);
                const auto &pt_coords = points.get_data_direction(dir);
                for (int order = 0; order < max_der; ++order)
                    bernstein_values[comp].get_derivative(order) =
                        BernsteinBasis::derivative(order, deg, pt_coords);
            }

            const auto &inter_lengths = lengths.get_data_direction(dir);
            for (int inter = inter_begin ; inter < inter_end ; ++inter)
            {
            	auto &splines1d = g_cache.entry(dir, inter);
                for (auto comp : splines1d.get_active_components_id())
                {
                    const auto &berns_values = bernstein_values[comp];
                    auto &basis = splines1d[comp];
                    const auto &oper = bezier_op.get_operator(dir, inter, comp);

                    const Real one_div_size = 1.0 / inter_lengths[inter];
                    fill_interval_values(one_div_size, oper, berns_values, basis);
                }
            }

        }
    }
}



template<int dim_, int range_ , int rank_>
template<int k>
void
BSplineElementHandler<dim_, range_, rank_>::
init_cache(ElementAccessor &elem)
{
    base_t::template init_cache<k>(elem);

    auto &cache = elem.local_cache_;
    if (cache == nullptr)
    {
        using Cache = typename ElementAccessor::LocalCache;
        cache = shared_ptr<Cache>(new Cache);
    }

    const auto n_basis = elem.get_num_basis();

    for (auto &s_id: UnitElement<dim>::template elems_ids<k>())
    {
        auto &s_cache = cache->template get_value_cache<k>(s_id);
        const auto n_points = this->template get_num_points<k>();
        s_cache.resize(flags_[k], n_points, n_basis);
    }
}



//template<int dim_, int range_ , int rank_>
//template<int k>
//void
//BSplineElementHandler<dim_, range_, rank_>::
//init_all_caches(ElementAccessor &elem)
//{
//    auto &cache = elem.local_cache_;
//    if (cache == nullptr)
//    {
//        using Cache = typename ElementAccessor::LocalCache;
//        cache = shared_ptr<Cache>(new Cache);
//    }
//    init_unif_caches(flags_[dim], std::get<dim>(quad_), cache->values_);
//}



template<int dim_, int range_ , int rank_>
void
BSplineElementHandler<dim_, range_, rank_>::
init_element_cache(ElementIterator &elem)
{
    init_cache<dim>(*elem);
}


template <int dim, int range, int rank>
void
BSplineElementHandler<dim, range, rank>::
copy_to_inactive_components_values(const vector<Index> &inactive_comp,
                                   const std::array<Index, n_components> &active_map,
                                   ValueTable<Value> &D_phi) const
{
    const Size n_points = D_phi.get_num_points();
    for (int comp : inactive_comp)
    {
        const auto act_comp = active_map[comp];
        const auto n_basis = n_basis_.get_component_size(comp);
        const Size act_offset = comp_offset_[act_comp];
        const Size offset     = comp_offset_[comp];
        for (Size basis_i = 0; basis_i < n_basis;  ++basis_i)
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
copy_to_inactive_components(const vector<Index> &inactive_comp,
                            const std::array<Index, n_components> &active_map,
                            ValueTable<Derivative<order>> &D_phi) const
{
    const Size n_points = D_phi.get_num_points();
    const Size n_ders = Derivative<order>::size;
    for (int comp : inactive_comp)
    {
        const auto act_comp = active_map[comp];
        const auto n_basis = n_basis_.get_component_size(comp);
        const Size act_offset = comp_offset_[act_comp];
        const Size offset     = comp_offset_[comp];
        for (Size basis_i = 0; basis_i < n_basis;  ++basis_i)
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
evaluate_bspline_values(
    const ComponentContainer<TensorProductFunctionEvaluator<dim>> &elem_values,
    ValueTable<Value> &D_phi) const
{
    const Size n_points = D_phi.get_num_points();
    const TensorIndex<dim> der_tensor_id; // [0,0,..,0] tensor index
    for (int comp : elem_values.get_active_components_id())
    {
        auto &values = elem_values[comp];
        const int total_n_basis = n_basis_.get_component_size(comp);
        const Size offset = comp_offset_[comp];

        for (int func_id = 0; func_id < total_n_basis; ++func_id)
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
    const Size n_points = D_phi.get_num_points();


    TensorFunctionDerivativesSymmetry<dim,order> sym;
    const auto n_der =  TensorFunctionDerivativesSymmetry<dim,order>::num_entries_eval;

    const auto &univariate_order = sym.univariate_order ;
    const auto &copy_indices = sym.copy_indices;

    for (int comp : elem_values.get_active_components_id())
    {
        auto &values = elem_values[comp];
        const int total_n_basis = n_basis_.get_component_size(comp);
        const Size offset = comp_offset_[comp];

        for (int func_id = 0; func_id < total_n_basis; ++func_id)
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
template <int k>
void
BSplineElementHandler<dim_, range_, rank_>::
fill_cache(ElementAccessor &elem, const int j)
{
    base_t::template fill_cache<k> (elem, j);

    auto &g_cache = std::get<k>(splines1d_)[j];

    Assert(elem.local_cache_ != nullptr, ExcNullPtr());
    auto &cache = elem.local_cache_->template get_value_cache<k>(j);

    const auto &index = elem.get_tensor_index();
    //const TensorIndex<k> active(UnitElement<dim>::template get_elem<k>(j).active_directions);

    auto &flags = cache.flags_handler_;
    auto val_1d = g_cache.get_element_values(index);
    if (flags.fill_values())
    {
        auto &values = cache.template get_der<0>();
        evaluate_bspline_values(val_1d, values);
        flags.set_values_filled(true);
    }
    if (flags.fill_gradients())
    {
        auto &values = cache.template get_der<1>();
        evaluate_bspline_derivatives<1>(val_1d, values);
        flags.set_gradients_filled(true);
    }
    if (flags.fill_hessians())
    {
        auto &values = cache.template get_der<2>();
        evaluate_bspline_derivatives<2>(val_1d, values);
        flags.set_hessians_filled(true);
    }

    cache.set_filled(true);
}



template<int dim_, int range_ , int rank_>
void
BSplineElementHandler<dim_, range_, rank_>::
fill_element_cache(ElementIterator &elem)
{
    fill_cache<dim>(*elem, 0);
}



//template <int dim_, int range_ , int rank_>
//template <int k>
//void
//BSplineElementHandler<dim_, range_, rank_>::
//fill_element_cache_(ElementAccessor &elem, const int j)
//{
//    base_t::template fill_cache<dim-k> (elem, j);
//
//    Assert(elem.local_cache_ != nullptr, ExcNullPtr());
//    auto &cache = elem.local_cache_->template get_value_cache<k>(j);
//
//    const auto &index = elem.get_tensor_index();
//    auto val_1d = splines1d_.get_element_values(index);
//    if (cache.flags_handler_.fill_values())
//    {
//        auto &values = cache.template get_der<0>();
//        evaluate_bspline_values(val_1d, values);
//        cache.flags_handler_.set_values_filled(true);
//    }
//    if (cache.flags_handler_.fill_gradients())
//    {
//        auto &values = cache.template get_der<1>();
//        evaluate_bspline_derivatives<1>(val_1d, values);
//        cache.flags_handler_.set_gradients_filled(true);
//    }
//    if (cache.flags_handler_.fill_hessians())
//    {
//        auto &values = cache.template get_der<2>();
//        evaluate_bspline_derivatives<2>(val_1d, values);
//        cache.flags_handler_.set_hessians_filled(true);
//    }
//
////    if (cache.flags_handler_.fill_divergences())
////    {
////        //TODO(pauletti, Sep 7, 2014): create a specialize exception
////        Assert(cache.flags_handler_.gradients_filled(),
////               ExcMessage("Divergence requires gradient to be filled."));
////
////        auto D1  = cache.D1phi_.begin();
////        auto div = cache.div_phi_.begin();
////        auto end = cache.D1phi_.end();
////        for (; D1 != end; ++D1, ++div)
////            *div = trace(*D1);
////
////        cache.flags_handler_.set_divergences_filled(true);
////    }
//
//    cache.set_filled(true);
//}


//template<int dim_, int range_ , int rank_>
//void
//BSplineElementHandler<dim_, range_, rank_>::
//fill_element_cache(ElementAccessor &elem)
//{
//    base_t::template fill_cache<dim>(elem, 0);
//}



template<int dim_, int range_ , int rank_>
void
BSplineElementHandler<dim_, range_, rank_>::
print_info(LogStream &out) const
{
    out.begin_item("Grid Cache:");
    base_t::print_info(out);
    out.end_item();

    cacheutils::print_caches(splines1d_, out);
    //out.begin_item("One dimensional splines cache:");
    //splines1d_.print_info(out);
    //out.end_item();
}

IGA_NAMESPACE_CLOSE

#include <igatools/basis_functions/bspline_element_handler.inst>
