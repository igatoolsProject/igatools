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

};

template<int dim_, int range_ , int rank_>
BSplineElementHandler<dim_, range_, rank_>::
BSplineElementHandler(shared_ptr<const Space> space,
                      const NewValueFlags flag,
                      const Quadrature<dim> &quad)
    :
    base_t(space->get_grid(), FunctionFlags::to_grid_flags(flag), quad),
    space_(space),
    n_basis_(space_->get_num_all_element_basis()),
    flags_ {flag, flag},
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
}



template<int dim_, int range_ , int rank_>
void
BSplineElementHandler<dim_, range_, rank_>::
init_element_cache(ElementAccessor &elem)
{
    base_t::init_element_cache(elem);

    auto &cache = elem.local_cache_;
    if (cache == nullptr)
    {
        using Cache = typename ElementAccessor::LocalCache;
        cache = shared_ptr<Cache>(new Cache);
    }

    auto n_basis = space_->get_num_all_element_basis();

    auto &elem_cache = cache->template get_value_cache<0>(0);
    elem_cache.resize(std::get<0>(flags_), quad_, n_basis);

    for (auto &f: base_t::faces)
    {
        auto &face_cache = cache->template get_value_cache<1>(f);
        face_cache.resize(std::get<1>(flags_), quad_.collapse_to_face(f), n_basis);
    }

}



template<int dim_, int range_ , int rank_>
void
BSplineElementHandler<dim_, range_, rank_>::
init_element_cache(ElementIterator &elem)
{
    init_element_cache(elem.get_accessor());
}


template <int dim, int range, int rank>
void
BSplineElementHandler<dim, range, rank>::
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
BSplineElementHandler<dim, range, rank>::
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
BSplineElementHandler<dim, range, rank>::
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
                const auto &copy_indices_der = copy_indices[der_id];
                const auto copy_indices_der_size = copy_indices_der.size();

                auto const &der_tensor_id = univariate_order[der_id];
                for (int point_id = 0; point_id < num_points; ++point_id)
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
void
BSplineElementHandler<dim_, range_, rank_>::
fill_element_cache(ElementIterator &elem)
{
    fill_element_cache(elem.get_accessor());
}



template <int dim_, int range_ , int rank_>
template <int k>
void
BSplineElementHandler<dim_, range_, rank_>::
fill_element_cache_(ElementAccessor &elem, const int j)
{
    base_t::template fill_element_cache_<k> (elem, j);

    Assert(elem.local_cache_ != nullptr, ExcNullPtr());
    auto &cache = elem.local_cache_->template get_value_cache<k>(j);

    const auto &index = elem.get_tensor_index();
    auto val_1d = splines1d_.get_element_values(index);
    if (cache.flags_handler_.fill_values())
    {
        auto &values = cache.template get_der<0>();
        evaluate_bspline_values(val_1d, values);
        cache.flags_handler_.set_values_filled(true);
    }
    if (cache.flags_handler_.fill_gradients())
    {
        auto &values = cache.template get_der<1>();
        evaluate_bspline_derivatives<1>(val_1d, values);
        cache.flags_handler_.set_gradients_filled(true);
    }
    if (cache.flags_handler_.fill_hessians())
    {
        auto &values = cache.template get_der<2>();
        evaluate_bspline_derivatives<2>(val_1d, values);
        cache.flags_handler_.set_hessians_filled(true);
    }

//    if (cache.flags_handler_.fill_divergences())
//    {
//        //TODO(pauletti, Sep 7, 2014): create a specialize exception
//        Assert(cache.flags_handler_.gradients_filled(),
//               ExcMessage("Divergence requires gradient to be filled."));
//
//        auto D1  = cache.D1phi_.begin();
//        auto div = cache.div_phi_.begin();
//        auto end = cache.D1phi_.end();
//        for (; D1 != end; ++D1, ++div)
//            *div = trace(*D1);
//
//        cache.flags_handler_.set_divergences_filled(true);
//    }

    cache.set_filled(true);
}


template<int dim_, int range_ , int rank_>
void
BSplineElementHandler<dim_, range_, rank_>::
fill_element_cache(ElementAccessor &elem)
{
    this->template fill_element_cache_<0>(elem, 0);
}



template<int dim_, int range_ , int rank_>
void
BSplineElementHandler<dim_, range_, rank_>::
print_info(LogStream &out) const
{
    out.begin_item("Grid Cache:");
    base_t::print_info(out);
    out.end_item();


    out.begin_item("One dimensional splines cache:");
    splines1d_.print_info(out);
    out.end_item();
}

IGA_NAMESPACE_CLOSE

#include <igatools/basis_functions/bspline_element_handler.inst>
